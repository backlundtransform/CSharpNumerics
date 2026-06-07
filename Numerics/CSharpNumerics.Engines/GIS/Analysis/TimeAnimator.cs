using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Simulation;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Analysis
{
    /// <summary>
    /// A time-stepped sequence of <see cref="ProbabilityMap"/> objects — one per
    /// time step. Answers questions such as
    /// <em>"what is the probability that gas has spread to position X within 30 minutes?"</em>
    /// </summary>
    public class TimeAnimator
    {
        private readonly ProbabilityMap[] _maps;
        private readonly MonteCarloScenarioResult _scenarios;
        private readonly int[] _iterationFilter;

        /// <summary>The spatial grid.</summary>
        public GeoGrid Grid { get; }

        /// <summary>The simulation time frame.</summary>
        public TimeFrame TimeFrame { get; }

        /// <summary>Concentration threshold used for exceedance.</summary>
        public double Threshold { get; }

        /// <summary>Number of time steps (and probability maps).</summary>
        public int Count => _maps.Length;

        /// <summary>Returns the <see cref="ProbabilityMap"/> at the given time-step index.</summary>
        public ProbabilityMap this[int timeIndex] => _maps[timeIndex];

        private TimeAnimator(
            ProbabilityMap[] maps,
            MonteCarloScenarioResult scenarios,
            int[] iterationFilter,
            double threshold)
        {
            _maps = maps;
            _scenarios = scenarios;
            _iterationFilter = iterationFilter;
            Grid = scenarios.Grid;
            TimeFrame = scenarios.TimeFrame;
            Threshold = threshold;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Point queries
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Instantaneous exceedance probability at a position and time.
        /// </summary>
        public double ProbabilityAt(Vector position, double timeSeconds)
        {
            int tIndex = TimeFrame.NearestIndex(timeSeconds);
            return _maps[tIndex].At(position);
        }

        /// <summary>
        /// Cumulative exceedance probability: fraction of scenarios where
        /// concentration exceeded the threshold at the given cell in
        /// <em>any</em> time step up to and including <paramref name="timeSeconds"/>.
        /// </summary>
        public double CumulativeProbabilityAt(Vector position, double timeSeconds)
        {
            int tIndex = TimeFrame.NearestIndex(timeSeconds);
            int cellIndex = Grid.NearestFlatIndex(position);
            int cellCount = Grid.CellCount;

            int hitCount = 0;
            int iterCount = _iterationFilter.Length;

            for (int k = 0; k < iterCount; k++)
            {
                int iter = _iterationFilter[k];
                bool exceeded = false;
                for (int t = 0; t <= tIndex; t++)
                {
                    int col = t * cellCount + cellIndex;
                    if (_scenarios.ScenarioMatrix.values[iter, col] > Threshold)
                    {
                        exceeded = true;
                        break;
                    }
                }
                if (exceeded) hitCount++;
            }

            return iterCount > 0 ? (double)hitCount / iterCount : 0;
        }

        /// <summary>
        /// Returns all <see cref="ProbabilityMap"/> objects as an array.
        /// </summary>
        public ProbabilityMap[] ToArray() => (ProbabilityMap[])_maps.Clone();

        // ═══════════════════════════════════════════════════════════════
        //  Factory
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Builds a time-animated probability sequence.
        /// </summary>
        /// <param name="scenarios">Monte Carlo result to analyse.</param>
        /// <param name="threshold">Concentration threshold (kg/m³).</param>
        /// <param name="iterationFilter">
        /// Optional subset of iteration indices (e.g. from
        /// <see cref="ClusterAnalysisResult.GetClusterIterations"/>).
        /// When null, all iterations are used.
        /// </param>
        public static TimeAnimator Build(
            MonteCarloScenarioResult scenarios,
            double threshold,
            int[] iterationFilter = null)
        {
            if (scenarios == null) throw new ArgumentNullException(nameof(scenarios));

            // Materialise the iteration list for cumulative queries
            if (iterationFilter == null)
            {
                iterationFilter = new int[scenarios.Iterations];
                for (int i = 0; i < scenarios.Iterations; i++)
                    iterationFilter[i] = i;
            }

            int timeCount = scenarios.TimeFrame.Count;
            var maps = new ProbabilityMap[timeCount];

            for (int t = 0; t < timeCount; t++)
                maps[t] = ProbabilityMap.Build(scenarios, t, threshold, iterationFilter);

            return new TimeAnimator(maps, scenarios, iterationFilter, threshold);
        }
    }
}
