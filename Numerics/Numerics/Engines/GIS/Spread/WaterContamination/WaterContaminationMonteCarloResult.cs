using CSharpNumerics.Engines.GIS.Analysis;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Simulation;
using CSharpNumerics.ML.Clustering.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread.WaterContamination;

/// <summary>
/// Result of a Monte Carlo water contamination ensemble.
/// Contains per-cell exceedance probability and per-iteration peak concentration statistics.
/// </summary>
public class WaterContaminationMonteCarloResult
{
    /// <summary>Number of Monte Carlo iterations.</summary>
    public int Iterations { get; }

    /// <summary>The spatial grid.</summary>
    public GeoGrid Grid { get; }

    /// <summary>The simulation time frame.</summary>
    public TimeFrame TimeFrame { get; }

    /// <summary>
    /// Per-cell exceedance probability (fraction 0–1).
    /// Value = fraction of iterations where the cell's peak concentration exceeded the toxicity threshold.
    /// </summary>
    public double[] ExceedanceProbability { get; }

    /// <summary>Peak concentration (mg/L) per iteration.</summary>
    public double[] PeakConcentrations { get; }

    /// <summary>Mean peak concentration across all iterations (mg/L).</summary>
    public double MeanPeakConcentration
    {
        get
        {
            double sum = 0;
            for (int i = 0; i < PeakConcentrations.Length; i++) sum += PeakConcentrations[i];
            return PeakConcentrations.Length > 0 ? sum / PeakConcentrations.Length : 0;
        }
    }

    /// <summary>
    /// P95 peak concentration: the 95th percentile across all iterations.
    /// </summary>
    public double P95PeakConcentration
    {
        get
        {
            if (PeakConcentrations.Length == 0) return 0;
            var sorted = new double[PeakConcentrations.Length];
            Array.Copy(PeakConcentrations, sorted, PeakConcentrations.Length);
            Array.Sort(sorted);
            int idx = (int)Math.Ceiling(0.95 * sorted.Length) - 1;
            return sorted[Math.Max(0, idx)];
        }
    }

    /// <summary>Worst-case (earliest) plume arrival time across all iterations (seconds).</summary>
    public double WorstCaseArrivalTime { get; }

    /// <summary>
    /// Per-iteration snapshot lists for optional downstream analysis (clustering).
    /// </summary>
    public IReadOnlyList<IReadOnlyList<SpreadSnapshot>> AllSnapshots { get; }

    /// <summary>
    /// Creates a Monte Carlo water contamination result.
    /// </summary>
    public WaterContaminationMonteCarloResult(
        int iterations,
        GeoGrid grid,
        double[] exceedanceProbability,
        double[] peakConcentrations,
        double worstCaseArrivalTime,
        IReadOnlyList<IReadOnlyList<SpreadSnapshot>> allSnapshots,
        TimeFrame timeFrame)
    {
        Iterations = iterations;
        Grid = grid ?? throw new ArgumentNullException(nameof(grid));
        ExceedanceProbability = exceedanceProbability ?? throw new ArgumentNullException(nameof(exceedanceProbability));
        PeakConcentrations = peakConcentrations ?? throw new ArgumentNullException(nameof(peakConcentrations));
        WorstCaseArrivalTime = worstCaseArrivalTime;
        AllSnapshots = allSnapshots;
        TimeFrame = timeFrame ?? throw new ArgumentNullException(nameof(timeFrame));
    }

    // ═══════════════════════════════════════════════════════════════
    //  Clustering integration
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Converts this water contamination result into a <see cref="MonteCarloScenarioResult"/>
    /// that can be consumed by <see cref="ScenarioClusterAnalyzer"/>.
    /// </summary>
    public MonteCarloScenarioResult ToMonteCarloScenarioResult(string layerName = "concentration")
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
    /// Clusters the Monte Carlo iterations by their spatiotemporal contamination patterns
    /// and returns a <see cref="WaterContaminationAnalysisResult"/>.
    /// </summary>
    public WaterContaminationAnalysisResult AnalyzeWith(
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
        return new WaterContaminationAnalysisResult(this, analysis);
    }
}
