using CSharpNumerics.Engines.GIS.Analysis;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Spread.Wildfire.Enums;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread.Wildfire;

/// <summary>
/// Result of clustering a wildfire Monte Carlo ensemble. Wraps the
/// <see cref="ClusterAnalysisResult"/> and provides fire-specific
/// accessors and a <see cref="Build"/> method that produces a
/// probability-weighted <see cref="WildfireScenarioResult"/>.
/// </summary>
public class WildfireAnalysisResult
{
    /// <summary>The original Monte Carlo wildfire result.</summary>
    public WildfireMonteCarloResult MonteCarloResult { get; }

    /// <summary>The cluster analysis (labels, dominant cluster, scores).</summary>
    public ClusterAnalysisResult ClusterAnalysis { get; }

    internal WildfireAnalysisResult(
        WildfireMonteCarloResult mcResult,
        ClusterAnalysisResult clusterAnalysis)
    {
        MonteCarloResult = mcResult ?? throw new ArgumentNullException(nameof(mcResult));
        ClusterAnalysis = clusterAnalysis ?? throw new ArgumentNullException(nameof(clusterAnalysis));
    }

    /// <summary>
    /// Computes per-cell burn probability using only the iterations assigned
    /// to the specified cluster.
    /// </summary>
    public double[] GetClusterBurnProbability(int cluster)
    {
        var iterations = ClusterAnalysis.GetClusterIterations(cluster);
        int cellCount = MonteCarloResult.Grid.CellCount;

        if (iterations.Length == 0)
            return new double[cellCount];

        var burnCount = new int[cellCount];

        foreach (int iter in iterations)
        {
            var snapshots = MonteCarloResult.AllSnapshots[iter];
            var last = snapshots[snapshots.Count - 1];
            var bs = last.Snapshot.GetLayer("burnState");
            for (int c = 0; c < cellCount; c++)
            {
                int state = (int)bs[c];
                if (state == (int)CellBurnState.Burning || state == (int)CellBurnState.Burned)
                    burnCount[c]++;
            }
        }

        var prob = new double[cellCount];
        for (int c = 0; c < cellCount; c++)
            prob[c] = (double)burnCount[c] / iterations.Length;
        return prob;
    }

    /// <summary>
    /// Builds a <see cref="WildfireScenarioResult"/> from the dominant cluster.
    /// Each snapshot contains mean values computed across the dominant cluster's
    /// iterations. The burnState layer holds burn probability (0–1) rather than
    /// discrete states.
    /// </summary>
    public WildfireScenarioResult Build()
    {
        var dominantIters = ClusterAnalysis.GetClusterIterations(ClusterAnalysis.DominantCluster);
        int timeCount = MonteCarloResult.AllSnapshots[0].Count;
        int cellCount = MonteCarloResult.Grid.CellCount;
        var grid = MonteCarloResult.Grid;
        double n = dominantIters.Length;

        var snapshots = new List<SpreadSnapshot>(timeCount);

        for (int t = 0; t < timeCount; t++)
        {
            var meanBurnState = new double[cellCount];
            var meanFlameLength = new double[cellCount];
            var meanROS = new double[cellCount];
            var meanBurnTime = new double[cellCount];

            foreach (int iter in dominantIters)
            {
                var snap = MonteCarloResult.AllSnapshots[iter][t];
                var bs = snap.Snapshot.GetLayer("burnState");
                var fl = snap.Snapshot.GetLayer("flameLength");
                var ros = snap.Snapshot.GetLayer("rateOfSpread");
                var bt = snap.Snapshot.GetLayer("burnTime");

                for (int c = 0; c < cellCount; c++)
                {
                    int state = (int)bs[c];
                    if (state == (int)CellBurnState.Burning || state == (int)CellBurnState.Burned)
                        meanBurnState[c] += 1.0;
                    meanFlameLength[c] += fl[c];
                    meanROS[c] += ros[c];
                    meanBurnTime[c] += bt[c];
                }
            }

            for (int c = 0; c < cellCount; c++)
            {
                meanBurnState[c] /= n;
                meanFlameLength[c] /= n;
                meanROS[c] /= n;
                meanBurnTime[c] /= n;
            }

            double time = MonteCarloResult.AllSnapshots[0][t].Time;
            snapshots.Add(new SpreadSnapshot(
                grid, time, t,
                meanBurnState, meanFlameLength, meanROS, meanBurnTime));
        }

        return new WildfireScenarioResult(snapshots, grid);
    }
}
