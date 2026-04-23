using CSharpNumerics.Engines.GIS.Analysis;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Simulation;
using CSharpNumerics.ML.Clustering;
using CSharpNumerics.ML.Clustering.Interfaces;
using CSharpNumerics.ML.Clustering.Results;
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
    /// <para>
    /// Internally filters the scenario matrix to only include river cells that
    /// carried non-zero concentration in at least one iteration, improving cluster
    /// separation in the high-dimensional feature space.
    /// </para>
    /// </summary>
    public WaterContaminationAnalysisResult AnalyzeWith(
        IClusteringModel algorithm,
        IClusteringEvaluator evaluator,
        int minK = 2,
        int maxK = 6)
    {
        var mcResult = ToMonteCarloScenarioResult();
        var clusterMatrix = BuildFilteredMatrix(mcResult);

        var experiment = ClusteringExperiment
            .For(clusterMatrix)
            .WithAlgorithm(algorithm)
            .TryClusterCounts(minK, maxK)
            .WithEvaluator(evaluator)
            .Run();

        var analysis = new ClusterAnalysisResult(mcResult, experiment);
        return new WaterContaminationAnalysisResult(this, analysis);
    }

    /// <summary>
    /// Clusters the Monte Carlo iterations using a full
    /// <see cref="ClusteringGrid"/> hyperparameter search.
    /// <para>
    /// The same active-river-cell filtering is applied before clustering
    /// to avoid diluting euclidean distances with zero-only non-river cells.
    /// </para>
    /// </summary>
    public WaterContaminationAnalysisResult AnalyzeWith(
        ClusteringGrid grid,
        IClusteringEvaluator evaluator,
        int minK = 2,
        int maxK = 6)
    {
        var mcResult = ToMonteCarloScenarioResult();
        var clusterMatrix = BuildFilteredMatrix(mcResult);

        var experiment = ClusteringExperiment
            .For(clusterMatrix)
            .WithGrid(grid)
            .TryClusterCounts(minK, maxK)
            .WithEvaluator(evaluator)
            .Run();

        var analysis = new ClusterAnalysisResult(mcResult, experiment);
        return new WaterContaminationAnalysisResult(this, analysis);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Internal — active-cell filtering
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Builds a feature matrix containing only columns for river cells
    /// that carried non-zero concentration in at least one iteration.
    /// Falls back to the full matrix if every cell is active.
    /// </summary>
    private Matrix BuildFilteredMatrix(MonteCarloScenarioResult mcResult)
    {
        int cellCount = Grid.CellCount;
        int timeCount = AllSnapshots[0].Count;

        var activeCells = new List<int>();
        for (int c = 0; c < cellCount; c++)
        {
            bool active = false;
            for (int i = 0; i < Iterations && !active; i++)
                for (int t = 0; t < timeCount && !active; t++)
                    if (mcResult.ScenarioMatrix.values[i, t * cellCount + c] > 0)
                        active = true;
            if (active) activeCells.Add(c);
        }

        if (activeCells.Count > 0 && activeCells.Count < cellCount)
        {
            int filtCols = activeCells.Count * timeCount;
            var filt = new double[Iterations, filtCols];
            for (int i = 0; i < Iterations; i++)
            {
                for (int t = 0; t < timeCount; t++)
                {
                    int srcOff = t * cellCount;
                    int dstOff = t * activeCells.Count;
                    for (int ac = 0; ac < activeCells.Count; ac++)
                        filt[i, dstOff + ac] = mcResult.ScenarioMatrix.values[i, srcOff + activeCells[ac]];
                }
            }
            return new Matrix(filt);
        }

        return mcResult.ScenarioMatrix;
    }
}
