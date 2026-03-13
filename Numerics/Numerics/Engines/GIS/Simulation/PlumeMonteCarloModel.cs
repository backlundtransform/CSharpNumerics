using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Enums;
using CSharpNumerics.Statistics.MonteCarlo;
using CSharpNumerics.Statistics.Random;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Simulation
{
    /// <summary>
    /// Monte Carlo adapter for plume simulations. Samples physics parameters
    /// from a <see cref="ScenarioVariation"/> and runs a <see cref="PlumeSimulator"/>
    /// per trial.
    /// <para>
    /// Provides two modes of operation:
    /// <list type="bullet">
    ///   <item><see cref="RunBatch"/> — runs N scenarios and returns the full
    ///   scenario matrix (rows = iterations, columns = cells × time steps).</item>
    ///   <item><see cref="IMonteCarloModel.Evaluate"/> — returns a scalar summary
    ///   (peak ground-level concentration) compatible with <see cref="MonteCarloSimulator"/>.</item>
    /// </list>
    /// </para>
    /// </summary>
    public class PlumeMonteCarloModel : IMonteCarloModel
    {
        // ── Baseline parameters ──────────────────────────────────────

        private readonly double _emissionRate;
        private readonly double _windSpeed;
        private readonly Vector _windDirection;
        private readonly double _stackHeight;
        private readonly Vector _sourcePosition;
        private readonly StabilityClass _stability;
        private readonly PlumeMode _mode;
        private readonly double _releaseSeconds;

        // ── Variation and grid ───────────────────────────────────────

        private readonly ScenarioVariation _variation;
        private readonly GeoGrid _grid;
        private readonly TimeFrame _timeFrame;

        /// <summary>
        /// Creates a Monte Carlo plume model.
        /// </summary>
        /// <param name="emissionRate">Baseline emission rate Q (kg/s).</param>
        /// <param name="windSpeed">Baseline wind speed u (m/s).</param>
        /// <param name="windDirection">Baseline wind direction.</param>
        /// <param name="stackHeight">Effective stack height H (m).</param>
        /// <param name="sourcePosition">Source position.</param>
        /// <param name="grid">Spatial grid to evaluate on.</param>
        /// <param name="timeFrame">Time range.</param>
        /// <param name="variation">Parameter variation ranges (may be null for deterministic).</param>
        /// <param name="stability">Baseline stability class.</param>
        /// <param name="mode">Physics mode.</param>
        /// <param name="releaseSeconds">Release duration for transient mode.</param>
        public PlumeMonteCarloModel(
            double emissionRate,
            double windSpeed,
            Vector windDirection,
            double stackHeight,
            Vector sourcePosition,
            GeoGrid grid,
            TimeFrame timeFrame,
            ScenarioVariation variation = null,
            StabilityClass stability = StabilityClass.D,
            PlumeMode mode = PlumeMode.SteadyState,
            double releaseSeconds = 0)
        {
            _emissionRate = emissionRate;
            _windSpeed = windSpeed;
            _windDirection = windDirection;
            _stackHeight = stackHeight;
            _sourcePosition = sourcePosition;
            _grid = grid;
            _timeFrame = timeFrame;
            _variation = variation;
            _stability = stability;
            _mode = mode;
            _releaseSeconds = releaseSeconds;
        }

        // ═══════════════════════════════════════════════════════════════
        //  IMonteCarloModel — scalar summary (peak concentration)
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Runs one stochastic trial and returns the peak ground-level
        /// concentration across all cells and time steps as a scalar.
        /// Compatible with <see cref="MonteCarloSimulator"/>.
        /// </summary>
        public double Evaluate(RandomGenerator rng)
        {
            var snapshots = RunOneScenario(rng);
            double max = 0;
            foreach (var snap in snapshots)
            {
                double snapMax = snap.Max();
                if (snapMax > max) max = snapMax;
            }
            return max;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Batch runner — full scenario matrix
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Runs N stochastic scenarios and returns a <see cref="MonteCarloScenarioResult"/>
        /// containing the full scenario matrix plus individual snapshots.
        /// </summary>
        /// <param name="iterations">Number of Monte Carlo iterations.</param>
        /// <param name="seed">Optional seed for reproducibility.</param>
        public MonteCarloScenarioResult RunBatch(int iterations, int? seed = null)
        {
            if (iterations <= 0) throw new ArgumentException("iterations must be positive.");

            var rng = seed.HasValue ? new RandomGenerator(seed.Value) : new RandomGenerator();
            int cols = _grid.CellCount * _timeFrame.Count;

            var scenarioMatrix = new double[iterations, cols];
            var allSnapshots = new List<List<GridSnapshot>>(iterations);

            for (int i = 0; i < iterations; i++)
            {
                var snapshots = RunOneScenario(rng);
                allSnapshots.Add(snapshots);

                // Flatten snapshots into one row
                int offset = 0;
                foreach (var snap in snapshots)
                {
                    var values = snap.GetValues();
                    for (int j = 0; j < values.Length; j++)
                        scenarioMatrix[i, offset + j] = values[j];
                    offset += values.Length;
                }
            }

            return new MonteCarloScenarioResult(
                new Matrix(scenarioMatrix),
                allSnapshots,
                _grid,
                _timeFrame);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Internal — single scenario with sampled parameters
        // ═══════════════════════════════════════════════════════════════

        private List<GridSnapshot> RunOneScenario(RandomGenerator rng)
        {
            double emissionRate = _emissionRate;
            double windSpeed = _windSpeed;
            Vector windDir = _windDirection;
            StabilityClass stability = _stability;

            if (_variation != null && _variation.HasVariation)
            {
                // Sample wind speed
                if (_variation.WindSpeedMin.HasValue)
                    windSpeed = rng.NextUniform(_variation.WindSpeedMin.Value, _variation.WindSpeedMax.Value);

                // Sample emission rate
                if (_variation.EmissionRateMin.HasValue)
                    emissionRate = rng.NextUniform(_variation.EmissionRateMin.Value, _variation.EmissionRateMax.Value);

                // Sample wind direction jitter (rotation around z-axis)
                if (_variation.WindDirectionJitterDeg.HasValue && _variation.WindDirectionJitterDeg.Value > 0)
                {
                    double jitterRad = rng.NextGaussian(0, _variation.WindDirectionJitterDeg.Value) * Math.PI / 180.0;
                    double cos = Math.Cos(jitterRad);
                    double sin = Math.Sin(jitterRad);
                    windDir = new Vector(
                        windDir.x * cos - windDir.y * sin,
                        windDir.x * sin + windDir.y * cos,
                        0);
                }

                // Sample stability class
                if (_variation.StabilityWeights != null)
                    stability = SampleStability(rng, _variation.StabilityWeights);
            }

            // Clamp wind speed to avoid invalid simulator input
            if (windSpeed <= 0) windSpeed = 0.1;
            if (emissionRate <= 0) emissionRate = 0.01;

            var sim = new PlumeSimulator(
                emissionRate, windSpeed, windDir, _stackHeight,
                _sourcePosition, stability, _mode);
            sim.ReleaseSeconds = _releaseSeconds;

            return sim.Run(_grid, _timeFrame);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Sampling helpers
        // ═══════════════════════════════════════════════════════════════

        private static StabilityClass SampleStability(
            RandomGenerator rng, Dictionary<StabilityClass, double> weights)
        {
            double total = 0;
            foreach (var w in weights.Values)
                total += w;

            if (total <= 0) return StabilityClass.D;

            double u = rng.NextDouble() * total;
            double cumulative = 0;
            foreach (var kv in weights)
            {
                cumulative += kv.Value;
                if (u <= cumulative)
                    return kv.Key;
            }

            // Fallback (shouldn't reach)
            return StabilityClass.D;
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Result container
    // ═══════════════════════════════════════════════════════════════════

    /// <summary>
    /// Holds the output of a batch Monte Carlo plume simulation:
    /// the full scenario matrix, per-iteration snapshots, and grid/time metadata.
    /// </summary>
    public class MonteCarloScenarioResult
    {
        /// <summary>
        /// Scenario matrix of shape (iterations, cells × timeSteps).
        /// Each row is one scenario's flattened concentration values.
        /// </summary>
        public Matrix ScenarioMatrix { get; }

        /// <summary>
        /// All grid snapshots per iteration: outer list = iterations, inner = time steps.
        /// </summary>
        public List<List<GridSnapshot>> Snapshots { get; }

        /// <summary>The spatial grid used.</summary>
        public GeoGrid Grid { get; }

        /// <summary>The time frame used.</summary>
        public TimeFrame TimeFrame { get; }

        /// <summary>Number of iterations (rows in the matrix).</summary>
        public int Iterations => ScenarioMatrix.rowLength;

        /// <summary>Number of columns (cells × time steps).</summary>
        public int FeatureCount => ScenarioMatrix.columnLength;

        public MonteCarloScenarioResult(
            Matrix scenarioMatrix,
            List<List<GridSnapshot>> snapshots,
            GeoGrid grid,
            TimeFrame timeFrame)
        {
            ScenarioMatrix = scenarioMatrix;
            Snapshots = snapshots;
            Grid = grid;
            TimeFrame = timeFrame;
        }

        /// <summary>
        /// Extracts the concentration values for a single cell across all
        /// iterations at a specific time step. Useful for building per-cell
        /// probability distributions.
        /// </summary>
        /// <param name="cellIndex">Flat grid cell index.</param>
        /// <param name="timeIndex">Zero-based time step index.</param>
        public double[] GetCellDistribution(int cellIndex, int timeIndex)
        {
            int col = timeIndex * Grid.CellCount + cellIndex;
            var values = new double[Iterations];
            for (int i = 0; i < Iterations; i++)
                values[i] = ScenarioMatrix.values[i, col];
            return values;
        }

        /// <summary>
        /// Extracts the full flattened scenario vector for iteration <paramref name="iterationIndex"/>.
        /// </summary>
        public double[] GetScenarioVector(int iterationIndex)
        {
            var row = new double[FeatureCount];
            for (int j = 0; j < FeatureCount; j++)
                row[j] = ScenarioMatrix.values[iterationIndex, j];
            return row;
        }
    }
}
