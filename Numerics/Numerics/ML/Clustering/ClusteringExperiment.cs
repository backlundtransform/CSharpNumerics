using CSharpNumerics.ML.Clustering.Evaluators;
using CSharpNumerics.ML.Clustering.Interfaces;
using CSharpNumerics.ML.Clustering.Results;
using CSharpNumerics.ML.Scalers.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace CSharpNumerics.ML.Clustering;

/// <summary>
/// Fluent entry point for clustering experiments.
/// 
/// <code>
/// var result = ClusteringExperiment
///     .For(data)
///     .WithAlgorithm(new KMeans())
///     .TryClusterCounts(2, 10)
///     .WithEvaluator(new SilhouetteEvaluator())
///     .WithScaler(new StandardScaler())
///     .Run();
/// </code>
/// </summary>
public static class ClusteringExperiment
{
    /// <summary>Start building a clustering experiment for the given data matrix.</summary>
    public static ClusteringExperimentBuilder For(Matrix data)
        => new ClusteringExperimentBuilder(data);
}

/// <summary>
/// Builder that configures and runs a clustering experiment.
/// Supports three modes:
///   1. Simple: WithAlgorithm + TryClusterCounts  (auto-expands K range)
///   2. Multi:  WithAlgorithms + TryClusterCounts  (multiple algorithms)
///   3. Grid:   WithGrid                            (full hyperparameter search)
/// </summary>
public class ClusteringExperimentBuilder
{
    private readonly Matrix _data;
    private readonly List<IClusteringModel> _algorithms = new();
    private readonly List<IClusteringEvaluator> _evaluators = new();
    private IScaler _scaler;
    private (int min, int max)? _clusterRange;
    private ClusteringGrid _grid;

    internal ClusteringExperimentBuilder(Matrix data)
    {
        _data = data;
    }

    // ── Algorithm configuration ──────────────────────────────────

    /// <summary>Add a single clustering algorithm.</summary>
    public ClusteringExperimentBuilder WithAlgorithm(IClusteringModel model)
    {
        _algorithms.Add(model);
        return this;
    }

    /// <summary>Add multiple clustering algorithms.</summary>
    public ClusteringExperimentBuilder WithAlgorithms(params IClusteringModel[] models)
    {
        _algorithms.AddRange(models);
        return this;
    }

    // ── Cluster count range ──────────────────────────────────────

    /// <summary>
    /// Try all cluster counts from min to max (inclusive).
    /// For algorithms that accept K (KMeans, Agglomerative, GMM), this
    /// generates a separate run per K value.
    /// For algorithms that discover K (DBSCAN), this is ignored.
    /// </summary>
    public ClusteringExperimentBuilder TryClusterCounts(int min, int max)
    {
        if (min < 2) min = 2;
        if (max < min) max = min;
        _clusterRange = (min, max);
        return this;
    }

    // ── Evaluators ───────────────────────────────────────────────

    /// <summary>Add a single evaluator (the first one added is the primary for ranking).</summary>
    public ClusteringExperimentBuilder WithEvaluator(IClusteringEvaluator evaluator)
    {
        _evaluators.Add(evaluator);
        return this;
    }

    /// <summary>Add multiple evaluators.</summary>
    public ClusteringExperimentBuilder WithEvaluators(params IClusteringEvaluator[] evaluators)
    {
        _evaluators.AddRange(evaluators);
        return this;
    }

    // ── Scaler ───────────────────────────────────────────────────

    /// <summary>Apply a scaler before clustering (optional).</summary>
    public ClusteringExperimentBuilder WithScaler(IScaler scaler)
    {
        _scaler = scaler;
        return this;
    }

    // ── Grid (advanced) ──────────────────────────────────────────

    /// <summary>
    /// Use a full <see cref="ClusteringGrid"/> for hyperparameter search.
    /// When set, WithAlgorithm/TryClusterCounts/WithScaler are ignored.
    /// </summary>
    public ClusteringExperimentBuilder WithGrid(ClusteringGrid grid)
    {
        _grid = grid;
        return this;
    }

    // ── Run ──────────────────────────────────────────────────────

    /// <summary>Execute the experiment and return ranked results.</summary>
    public ClusteringExperimentResult Run()
    {
        if (_evaluators.Count == 0)
            throw new InvalidOperationException(
                "At least one evaluator is required. Call WithEvaluator() before Run().");

        var totalSw = Stopwatch.StartNew();

        // Build list of pipelines to evaluate
        var pipelines = _grid != null
            ? _grid.Expand().ToList()
            : BuildPipelinesFromConfig();

        // Run each pipeline
        var results = new List<ClusteringResult>();
        InertiaEvaluator inertiaEval = _evaluators.OfType<InertiaEvaluator>().FirstOrDefault();
        var elbowCurve = new List<(int K, double Inertia)>();

        foreach (var pipeline in pipelines)
        {
            var sw = Stopwatch.StartNew();
            VectorN labels = pipeline.FitPredict(_data);
            sw.Stop();

            // Evaluate with all evaluators
            var scores = new Dictionary<string, double>();
            foreach (var evaluator in _evaluators)
                scores[evaluator.Name] = evaluator.Score(_data, labels);

            var result = new ClusteringResult
            {
                Model = pipeline.Model,
                AlgorithmName = pipeline.Model.GetType().Name,
                ClusterCount = pipeline.Model.ClusterCount,
                Labels = labels,
                Scores = scores,
                Parameters = new Dictionary<string, object>(pipeline.ModelParams),
                PipelineDescription = pipeline.ToString(),
                Duration = sw.Elapsed
            };

            results.Add(result);

            // Track elbow curve
            if (inertiaEval != null)
                elbowCurve.Add((result.ClusterCount, inertiaEval.RawInertia(_data, labels)));
        }

        totalSw.Stop();

        // Rank by primary evaluator (first one added)
        string primaryName = _evaluators[0].Name;
        var ranked = results
            .OrderByDescending(r => r.Scores.ContainsKey(primaryName) ? r.Scores[primaryName] : double.NegativeInfinity)
            .ToList();

        return new ClusteringExperimentResult
        {
            Rankings = ranked,
            PrimaryEvaluator = primaryName,
            ElbowCurve = elbowCurve
                .OrderBy(e => e.K)
                .ToList(),
            TotalDuration = totalSw.Elapsed
        };
    }

    // ═══════════════════════════════════════════════════════════════
    //  Private: build pipelines from simple/multi config
    // ═══════════════════════════════════════════════════════════════

    private List<ClusteringPipeline> BuildPipelinesFromConfig()
    {
        if (_algorithms.Count == 0)
            throw new InvalidOperationException(
                "At least one algorithm is required. Call WithAlgorithm() or WithGrid() before Run().");

        var pipelines = new List<ClusteringPipeline>();

        foreach (var algorithm in _algorithms)
        {
            if (AcceptsK(algorithm) && _clusterRange.HasValue)
            {
                // Expand into one pipeline per K value
                for (int k = _clusterRange.Value.min; k <= _clusterRange.Value.max; k++)
                {
                    var clone = algorithm.Clone();
                    var kParams = new Dictionary<string, object> { ["K"] = k };

                    pipelines.Add(new ClusteringPipeline(
                        model: clone,
                        modelParams: kParams,
                        scaler: _scaler?.Clone()
                    ));
                }
            }
            else
            {
                // Single run (DBSCAN discovers K on its own)
                pipelines.Add(new ClusteringPipeline(
                    model: algorithm.Clone(),
                    modelParams: new Dictionary<string, object>(),
                    scaler: _scaler?.Clone()
                ));
            }
        }

        return pipelines;
    }

    /// <summary>
    /// Check if the algorithm has a settable K property
    /// (KMeans, AgglomerativeClustering, GMM — yes; DBSCAN — no).
    /// </summary>
    private static bool AcceptsK(IClusteringModel model)
    {
        var prop = model.GetType().GetProperty("K");
        return prop != null && prop.CanWrite;
    }
}
