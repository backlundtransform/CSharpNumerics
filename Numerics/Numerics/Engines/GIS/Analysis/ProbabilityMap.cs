using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Simulation;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Analysis
{
    /// <summary>
    /// Per-cell exceedance probability map: for each cell in a <see cref="GeoGrid"/>,
    /// the fraction of Monte Carlo iterations where the concentration exceeds a
    /// threshold at a given time step.
    /// <para>
    /// Optionally filtered to a subset of iterations (e.g. a single cluster from
    /// <see cref="ClusterAnalysisResult.GetClusterIterations"/>).
    /// </para>
    /// </summary>
    public class ProbabilityMap
    {
        private readonly double[] _probabilities;

        /// <summary>The spatial grid.</summary>
        public GeoGrid Grid { get; }

        /// <summary>Simulation time in seconds for this map.</summary>
        public double Time { get; }

        /// <summary>Zero-based time-step index.</summary>
        public int TimeIndex { get; }

        /// <summary>Concentration threshold used to compute exceedance.</summary>
        public double Threshold { get; }

        /// <summary>Number of scenarios considered (all or filtered).</summary>
        public int ScenarioCount { get; }

        /// <summary>Number of cells.</summary>
        public int CellCount => _probabilities.Length;

        private ProbabilityMap(GeoGrid grid, double[] probabilities,
            double time, int timeIndex, double threshold, int scenarioCount)
        {
            Grid = grid;
            _probabilities = probabilities;
            Time = time;
            TimeIndex = timeIndex;
            Threshold = threshold;
            ScenarioCount = scenarioCount;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Value access
        // ═══════════════════════════════════════════════════════════════

        /// <summary>Exceedance probability at a flat cell index.</summary>
        public double this[int flatIndex] => _probabilities[flatIndex];

        /// <summary>Exceedance probability at 3-D cell indices.</summary>
        public double this[int ix, int iy, int iz] => _probabilities[Grid.Index(ix, iy, iz)];

        /// <summary>Exceedance probability at the cell nearest to <paramref name="position"/>.</summary>
        public double At(Vector position) => _probabilities[Grid.NearestFlatIndex(position)];

        /// <summary>Returns a copy of the probability array.</summary>
        public double[] GetValues() => (double[])_probabilities.Clone();

        // ═══════════════════════════════════════════════════════════════
        //  Queries
        // ═══════════════════════════════════════════════════════════════

        /// <summary>Maximum probability in the map.</summary>
        public double Max()
        {
            double max = 0;
            for (int i = 0; i < _probabilities.Length; i++)
                if (_probabilities[i] > max) max = _probabilities[i];
            return max;
        }

        /// <summary>
        /// Returns cells where the exceedance probability is above <paramref name="pThreshold"/>.
        /// </summary>
        public List<GeoCell> CellsAbove(double pThreshold)
        {
            var result = new List<GeoCell>();
            for (int i = 0; i < _probabilities.Length; i++)
            {
                if (_probabilities[i] > pThreshold)
                    result.Add(new GeoCell(Grid.CellCentre(i), _probabilities[i], TimeIndex, i));
            }
            return result;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Factory
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Builds a probability map for one time step.
        /// For each cell, counts the fraction of iterations where the
        /// concentration exceeds <paramref name="threshold"/>.
        /// </summary>
        /// <param name="scenarios">Monte Carlo result to analyse.</param>
        /// <param name="timeIndex">Zero-based time-step index.</param>
        /// <param name="threshold">Concentration threshold (kg/m³).</param>
        /// <param name="iterationFilter">
        /// Optional subset of iteration indices to include (e.g. from
        /// <see cref="ClusterAnalysisResult.GetClusterIterations"/>).
        /// When null, all iterations are used.
        /// </param>
        public static ProbabilityMap Build(
            MonteCarloScenarioResult scenarios,
            int timeIndex,
            double threshold,
            int[] iterationFilter = null)
        {
            if (scenarios == null) throw new ArgumentNullException(nameof(scenarios));
            if (timeIndex < 0 || timeIndex >= scenarios.TimeFrame.Count)
                throw new ArgumentOutOfRangeException(nameof(timeIndex));

            var grid = scenarios.Grid;
            int cellCount = grid.CellCount;

            // Determine iterations to scan
            int iterCount;
            if (iterationFilter != null)
            {
                iterCount = iterationFilter.Length;
            }
            else
            {
                iterCount = scenarios.Iterations;
                iterationFilter = new int[iterCount];
                for (int i = 0; i < iterCount; i++)
                    iterationFilter[i] = i;
            }

            var probs = new double[cellCount];

            if (iterCount == 0)
                return new ProbabilityMap(grid, probs,
                    scenarios.TimeFrame.TimeAt(timeIndex), timeIndex, threshold, 0);

            int colOffset = timeIndex * cellCount;
            for (int cell = 0; cell < cellCount; cell++)
            {
                int count = 0;
                int col = colOffset + cell;
                for (int k = 0; k < iterCount; k++)
                {
                    if (scenarios.ScenarioMatrix.values[iterationFilter[k], col] > threshold)
                        count++;
                }
                probs[cell] = (double)count / iterCount;
            }

            return new ProbabilityMap(grid, probs,
                scenarios.TimeFrame.TimeAt(timeIndex), timeIndex, threshold, iterCount);
        }
    }
}
