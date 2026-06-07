using CSharpNumerics.Engines.GIS.Analysis;
using CSharpNumerics.Engines.GIS.Grid;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread.WaterContamination;

/// <summary>
/// Result of clustering a water contamination Monte Carlo ensemble.
/// Wraps the <see cref="ClusterAnalysisResult"/> and provides contamination-specific
/// accessors and a <see cref="Build"/> method that produces a probability-weighted
/// <see cref="WaterContaminationResult"/>.
/// </summary>
public class WaterContaminationAnalysisResult
{
    /// <summary>The original Monte Carlo contamination result.</summary>
    public WaterContaminationMonteCarloResult MonteCarloResult { get; }

    /// <summary>The cluster analysis (labels, dominant cluster, scores).</summary>
    public ClusterAnalysisResult ClusterAnalysis { get; }

    internal WaterContaminationAnalysisResult(
        WaterContaminationMonteCarloResult mcResult,
        ClusterAnalysisResult clusterAnalysis)
    {
        MonteCarloResult = mcResult ?? throw new ArgumentNullException(nameof(mcResult));
        ClusterAnalysis = clusterAnalysis ?? throw new ArgumentNullException(nameof(clusterAnalysis));
    }

    /// <summary>
    /// Computes per-cell exceedance probability using only the iterations
    /// assigned to the specified cluster.
    /// </summary>
    public double[] GetClusterExceedanceProbability(int cluster, double threshold)
    {
        var iterations = ClusterAnalysis.GetClusterIterations(cluster);
        int cellCount = MonteCarloResult.Grid.CellCount;

        if (iterations.Length == 0)
            return new double[cellCount];

        var exceedCount = new int[cellCount];

        foreach (int iter in iterations)
        {
            var snapshots = MonteCarloResult.AllSnapshots[iter];
            for (int t = 0; t < snapshots.Count; t++)
            {
                var conc = snapshots[t].Snapshot.GetLayer("concentration");
                for (int c = 0; c < cellCount; c++)
                {
                    if (conc[c] > threshold)
                        exceedCount[c]++;
                }
            }
        }

        // Normalize: any exceedance across time in an iteration = 1 hit.
        // Re-process to get binary per-iteration exceedance.
        var binaryExceed = new int[cellCount];
        foreach (int iter in iterations)
        {
            var snapshots = MonteCarloResult.AllSnapshots[iter];
            var peaked = new bool[cellCount];
            for (int t = 0; t < snapshots.Count; t++)
            {
                var conc = snapshots[t].Snapshot.GetLayer("concentration");
                for (int c = 0; c < cellCount; c++)
                    if (conc[c] > threshold) peaked[c] = true;
            }
            for (int c = 0; c < cellCount; c++)
                if (peaked[c]) binaryExceed[c]++;
        }

        var prob = new double[cellCount];
        for (int c = 0; c < cellCount; c++)
            prob[c] = (double)binaryExceed[c] / iterations.Length;
        return prob;
    }

    /// <summary>
    /// Builds a <see cref="WaterContaminationResult"/> from the dominant cluster.
    /// Each snapshot contains mean concentration values computed across the
    /// dominant cluster's iterations.
    /// </summary>
    public WaterContaminationResult Build(double threshold = 0.01)
    {
        var dominantIters = ClusterAnalysis.GetClusterIterations(ClusterAnalysis.DominantCluster);
        int timeCount = MonteCarloResult.AllSnapshots[0].Count;
        int cellCount = MonteCarloResult.Grid.CellCount;
        var grid = MonteCarloResult.Grid;
        double n = dominantIters.Length;

        var snapshots = new List<SpreadSnapshot>(timeCount);

        for (int t = 0; t < timeCount; t++)
        {
            var meanConc = new double[cellCount];
            var meanVel = new double[cellCount];
            var meanState = new double[cellCount];
            var meanExposure = new double[cellCount];

            foreach (int iter in dominantIters)
            {
                var snap = MonteCarloResult.AllSnapshots[iter][t];
                var conc = snap.Snapshot.GetLayer("concentration");
                var vel = snap.Snapshot.GetLayer("velocity");
                var state = snap.Snapshot.GetLayer("contaminationState");
                var exp = snap.Snapshot.GetLayer("exposureTime");

                for (int c = 0; c < cellCount; c++)
                {
                    meanConc[c] += conc[c];
                    meanVel[c] += vel[c];
                    meanState[c] += state[c];
                    meanExposure[c] += exp[c];
                }
            }

            for (int c = 0; c < cellCount; c++)
            {
                meanConc[c] /= n;
                meanVel[c] /= n;
                meanState[c] /= n;
                meanExposure[c] /= n;
            }

            double time = MonteCarloResult.AllSnapshots[0][t].Time;
            var gs = new GridSnapshot(grid, meanConc, time, t);
            gs.SetLayer("concentration", meanConc);
            gs.SetLayer("contaminationState", meanState);
            gs.SetLayer("velocity", meanVel);
            gs.SetLayer("exposureTime", meanExposure);

            snapshots.Add(new SpreadSnapshot(gs));
        }

        return new WaterContaminationResult(snapshots.AsReadOnly(), grid, threshold);
    }
}
