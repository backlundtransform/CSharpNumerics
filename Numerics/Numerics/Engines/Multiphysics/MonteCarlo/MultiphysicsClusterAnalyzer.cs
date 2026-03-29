using CSharpNumerics.ML.Clustering;
using CSharpNumerics.ML.Clustering.Interfaces;
using CSharpNumerics.ML.Clustering.Results;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Engines.Multiphysics.MonteCarlo;

/// <summary>
/// Fluent entry point for clustering Monte Carlo multiphysics scenario results.
/// Wraps <see cref="ClusteringExperiment"/> to operate on
/// <see cref="MultiphysicsScenarioResult.ScenarioMatrix"/>.
/// <para>
/// Usage:
/// <code>
/// var analysis = MultiphysicsClusterAnalyzer
///     .For(mcResult)
///     .WithAlgorithm(new KMeans())
///     .TryClusterCounts(2, 5)
///     .WithEvaluator(new SilhouetteEvaluator())
///     .Run();
///
/// int dominantCluster = analysis.DominantCluster;
/// int[] representativeIters = analysis.GetClusterIterations(dominantCluster);
/// </code>
/// </para>
/// </summary>
public static class MultiphysicsClusterAnalyzer
{
    /// <summary>Start building a cluster analysis for Monte Carlo results.</summary>
    public static MultiphysicsClusterBuilder For(MultiphysicsScenarioResult scenarios)
        => new MultiphysicsClusterBuilder(scenarios);
}

/// <summary>
/// Builder that configures and runs a clustering experiment on the
/// scenario matrix from a Monte Carlo multiphysics simulation.
/// </summary>
public class MultiphysicsClusterBuilder
{
    private readonly MultiphysicsScenarioResult _scenarios;
    private readonly List<IClusteringModel> _algorithms = new();
    private readonly List<IClusteringEvaluator> _evaluators = new();
    private (int min, int max)? _clusterRange;
    private ClusteringGrid _grid;

    internal MultiphysicsClusterBuilder(MultiphysicsScenarioResult scenarios)
    {
        _scenarios = scenarios ?? throw new ArgumentNullException(nameof(scenarios));
    }

    /// <summary>Add a clustering algorithm.</summary>
    public MultiphysicsClusterBuilder WithAlgorithm(IClusteringModel model)
    {
        _algorithms.Add(model);
        return this;
    }

    /// <summary>Try cluster counts from min to max (inclusive).</summary>
    public MultiphysicsClusterBuilder TryClusterCounts(int min, int max)
    {
        _clusterRange = (min, max);
        return this;
    }

    /// <summary>Add a cluster-quality evaluator.</summary>
    public MultiphysicsClusterBuilder WithEvaluator(IClusteringEvaluator evaluator)
    {
        _evaluators.Add(evaluator);
        return this;
    }

    /// <summary>Use a full <see cref="ClusteringGrid"/> for hyper-parameter search.</summary>
    public MultiphysicsClusterBuilder WithGrid(ClusteringGrid grid)
    {
        _grid = grid;
        return this;
    }

    /// <summary>Execute the clustering experiment and return a domain-specific result.</summary>
    public MultiphysicsClusterResult Run()
    {
        var builder = ClusteringExperiment.For(_scenarios.ScenarioMatrix);

        if (_grid != null)
        {
            builder = builder.WithGrid(_grid);
        }
        else
        {
            foreach (var alg in _algorithms)
                builder = builder.WithAlgorithm(alg);
            if (_clusterRange.HasValue)
                builder = builder.TryClusterCounts(_clusterRange.Value.min, _clusterRange.Value.max);
        }

        foreach (var eval in _evaluators)
            builder = builder.WithEvaluator(eval);

        var experimentResult = builder.Run();

        return new MultiphysicsClusterResult(_scenarios, experimentResult);
    }
}

/// <summary>
/// Result of clustering Monte Carlo multiphysics scenarios.
/// Provides cluster labels, the dominant cluster, and iteration look-ups.
/// </summary>
public class MultiphysicsClusterResult
{
    /// <summary>The Monte Carlo scenario data that was clustered.</summary>
    public MultiphysicsScenarioResult Scenarios { get; }

    /// <summary>The full clustering experiment result (rankings, scores).</summary>
    public ClusteringExperimentResult ExperimentResult { get; }

    /// <summary>Cluster label assigned to each scenario iteration.</summary>
    public VectorN Labels => ExperimentResult.BestLabels;

    /// <summary>Number of clusters in the best solution.</summary>
    public int BestClusterCount => ExperimentResult.BestClusterCount;

    /// <summary>Primary evaluator score of the best solution.</summary>
    public double BestScore => ExperimentResult.BestScore;

    /// <summary>
    /// The cluster with the most assigned iterations.
    /// Represents the most probable outcome group.
    /// </summary>
    public int DominantCluster { get; }

    internal MultiphysicsClusterResult(
        MultiphysicsScenarioResult scenarios,
        ClusteringExperimentResult experimentResult)
    {
        Scenarios = scenarios;
        ExperimentResult = experimentResult;
        DominantCluster = FindDominantCluster();
    }

    /// <summary>
    /// Returns the iteration indices assigned to the given cluster.
    /// </summary>
    public int[] GetClusterIterations(int cluster)
    {
        var list = new List<int>();
        for (int i = 0; i < Labels.Length; i++)
            if ((int)Labels[i] == cluster)
                list.Add(i);
        return list.ToArray();
    }

    /// <summary>
    /// Returns the mean <see cref="SimulationResult"/> output values
    /// across all iterations in the given cluster.
    /// </summary>
    public double[] GetClusterMeanOutput(int cluster)
    {
        var iterations = GetClusterIterations(cluster);
        if (iterations.Length == 0) return Array.Empty<double>();

        int cols = Scenarios.FeatureCount;
        var mean = new double[cols];
        foreach (int iter in iterations)
        {
            for (int j = 0; j < cols; j++)
                mean[j] += Scenarios.ScenarioMatrix.values[iter, j];
        }
        for (int j = 0; j < cols; j++)
            mean[j] /= iterations.Length;

        return mean;
    }

    private int FindDominantCluster()
    {
        var counts = new Dictionary<int, int>();
        for (int i = 0; i < Labels.Length; i++)
        {
            int c = (int)Labels[i];
            if (!counts.ContainsKey(c)) counts[c] = 0;
            counts[c]++;
        }

        int dominant = 0;
        int maxCount = 0;
        foreach (var kvp in counts)
        {
            if (kvp.Value > maxCount)
            {
                maxCount = kvp.Value;
                dominant = kvp.Key;
            }
        }
        return dominant;
    }
}
