using CSharpNumerics.Engines.GIS.Analysis;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Simulation;
using CSharpNumerics.Engines.GIS.Spread.Wildfire.Enums;
using CSharpNumerics.ML.Clustering.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread.Wildfire;

/// <summary>
/// Result of a Monte Carlo wildfire ensemble.
/// Contains per-cell burn probability and per-iteration burned area statistics.
/// </summary>
public class WildfireMonteCarloResult
{
    /// <summary>Number of Monte Carlo iterations.</summary>
    public int Iterations { get; }

    /// <summary>The spatial grid.</summary>
    public GeoGrid Grid { get; }

    /// <summary>The simulation time frame.</summary>
    public TimeFrame TimeFrame { get; }

    /// <summary>
    /// Per-cell burn probability (fraction 0–1).
    /// Length = Nx × Ny. Value = fraction of iterations where the cell burned.
    /// </summary>
    public double[] BurnProbability { get; }

    /// <summary>Burned area (hectares) per iteration.</summary>
    public double[] BurnedAreas { get; }

    /// <summary>Mean burned area across all iterations (hectares).</summary>
    public double MeanBurnedArea
    {
        get
        {
            double sum = 0;
            for (int i = 0; i < BurnedAreas.Length; i++) sum += BurnedAreas[i];
            return BurnedAreas.Length > 0 ? sum / BurnedAreas.Length : 0;
        }
    }

    /// <summary>Maximum burned area across all iterations (hectares).</summary>
    public double MaxBurnedArea
    {
        get
        {
            double max = 0;
            for (int i = 0; i < BurnedAreas.Length; i++)
                if (BurnedAreas[i] > max) max = BurnedAreas[i];
            return max;
        }
    }

    /// <summary>
    /// Per-iteration snapshot lists for optional downstream analysis (clustering).
    /// </summary>
    public IReadOnlyList<IReadOnlyList<SpreadSnapshot>> AllSnapshots { get; }

    /// <summary>
    /// Per-cell boolean mask indicating which cells are firebreaks (water,
    /// no-fuel). Length = Nx × Ny. Derived from the first snapshot of the
    /// first iteration. Frontend consumers can use this to distinguish
    /// non-burnable terrain from cells that simply did not ignite.
    /// </summary>
    public bool[] FirebreakMask
    {
        get
        {
            if (AllSnapshots == null || AllSnapshots.Count == 0 || AllSnapshots[0].Count == 0)
                return new bool[Grid.CellCount];
            var bs = AllSnapshots[0][0].Snapshot.GetLayer("burnState");
            var mask = new bool[bs.Length];
            for (int i = 0; i < bs.Length; i++)
                mask[i] = (int)bs[i] == (int)CellBurnState.Firebreak;
            return mask;
        }
    }

    /// <summary>
    /// Creates a Monte Carlo wildfire result.
    /// </summary>
    public WildfireMonteCarloResult(
        int iterations,
        GeoGrid grid,
        double[] burnProbability,
        double[] burnedAreas,
        IReadOnlyList<IReadOnlyList<SpreadSnapshot>> allSnapshots,
        TimeFrame timeFrame)
    {
        Iterations = iterations;
        Grid = grid ?? throw new ArgumentNullException(nameof(grid));
        BurnProbability = burnProbability ?? throw new ArgumentNullException(nameof(burnProbability));
        BurnedAreas = burnedAreas ?? throw new ArgumentNullException(nameof(burnedAreas));
        AllSnapshots = allSnapshots;
        TimeFrame = timeFrame ?? throw new ArgumentNullException(nameof(timeFrame));
    }

    // ═══════════════════════════════════════════════════════════════
    //  Clustering integration
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Converts this wildfire result into a <see cref="MonteCarloScenarioResult"/>
    /// that can be consumed by <see cref="ScenarioClusterAnalyzer"/>.
    /// The scenario matrix is built from the specified layer across all
    /// iterations and time steps.
    /// </summary>
    /// <param name="layerName">Layer to extract (default: "burnState").</param>
    public MonteCarloScenarioResult ToMonteCarloScenarioResult(string layerName = "burnState")
    {
        int cellCount = Grid.CellCount;
        int timeCount = AllSnapshots[0].Count;
        int cols = cellCount * timeCount;

        var matrix = new double[Iterations, cols];
        var snapshots = new List<List<GridSnapshot>>(Iterations);

        for (int i = 0; i < Iterations; i++)
        {
            var iterSnaps = AllSnapshots[i];
            var gridSnaps = new List<GridSnapshot>(iterSnaps.Count);

            for (int t = 0; t < iterSnaps.Count; t++)
            {
                var values = iterSnaps[t].Snapshot.GetLayer(layerName);
                int offset = t * cellCount;
                for (int c = 0; c < cellCount; c++)
                    matrix[i, offset + c] = values[c];
                gridSnaps.Add(iterSnaps[t].Snapshot);
            }
            snapshots.Add(gridSnaps);
        }

        return new MonteCarloScenarioResult(
            new Matrix(matrix), snapshots, Grid, TimeFrame);
    }

    /// <summary>
    /// Clusters the Monte Carlo iterations by their spatiotemporal fire patterns
    /// using the existing <see cref="ScenarioClusterAnalyzer"/> pipeline.
    /// Returns a <see cref="WildfireAnalysisResult"/> from which a probability-weighted
    /// <see cref="WildfireScenarioResult"/> can be built.
    /// </summary>
    public WildfireAnalysisResult AnalyzeWith(
        IClusteringModel algorithm,
        IClusteringEvaluator evaluator,
        int minK = 2,
        int maxK = 6)
    {
        var mcResult = ToMonteCarloScenarioResult();
        var analysis = ScenarioClusterAnalyzer
            .For(mcResult)
            .WithAlgorithm(algorithm)
            .TryClusterCounts(minK, maxK)
            .WithEvaluator(evaluator)
            .Run();
        return new WildfireAnalysisResult(this, analysis);
    }
}
