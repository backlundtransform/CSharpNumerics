using CSharpNumerics.ML.CrossValidators.Interfaces;
using CSharpNumerics.ML.CrossValidators.Result;
using CSharpNumerics.ML.Experiment.Results;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace CSharpNumerics.ML.Experiment;

/// <summary>
/// Fluent entry point for supervised ML experiments.
/// Mirrors <see cref="CSharpNumerics.ML.Clustering.ClusteringExperiment"/>
/// but for supervised models (classification + regression) using cross-validation.
///
/// <code>
/// var result = SupervisedExperiment
///     .For(X, y)
///     .WithGrid(new PipelineGrid()
///         .AddModel&lt;KNearestNeighbors&gt;(g => g
///             .Add("K", 3, 5, 7)
///             .AddScaler&lt;StandardScaler&gt;(s => { }))
///         .AddModel&lt;DecisionTree&gt;(g => g
///             .Add("MaxDepth", 3, 5, 10)))
///     .WithCrossValidators(
///         CrossValidatorConfig.KFold(folds: 5),
///         CrossValidatorConfig.StratifiedKFold(folds: 10))
///     .Run();
/// </code>
/// </summary>
public static class SupervisedExperiment
{
    /// <summary>Start building a supervised experiment for the given data.</summary>
    public static SupervisedExperimentBuilder For(Matrix X, VectorN y)
        => new SupervisedExperimentBuilder(X, y);
}

/// <summary>
/// Builder that configures and runs a supervised ML experiment.
/// Supports two modes:
///   1. Simple:  WithPipeline(s) — manually constructed pipelines
///   2. Grid:    WithGrid        — full hyperparameter search via <see cref="PipelineGrid"/>
///
/// Cross-validators are added via <see cref="WithCrossValidator"/> or
/// <see cref="WithCrossValidators"/>. The first cross-validator added is
/// the primary one used for final ranking.
/// </summary>
public class SupervisedExperimentBuilder
{
    private readonly Matrix _X;
    private readonly VectorN _y;
    private readonly List<Pipeline> _pipelines = new();
    private readonly List<CrossValidatorConfig> _cvConfigs = new();
    private PipelineGrid _grid;

    internal SupervisedExperimentBuilder(Matrix X, VectorN y)
    {
        _X = X;
        _y = y;
    }

    // ── Pipeline configuration ───────────────────────────────────

    /// <summary>Add a single manually-constructed pipeline.</summary>
    public SupervisedExperimentBuilder WithPipeline(Pipeline pipeline)
    {
        _pipelines.Add(pipeline);
        return this;
    }

    /// <summary>Add multiple manually-constructed pipelines.</summary>
    public SupervisedExperimentBuilder WithPipelines(params Pipeline[] pipelines)
    {
        _pipelines.AddRange(pipelines);
        return this;
    }

    // ── Grid (advanced) ──────────────────────────────────────────

    /// <summary>
    /// Use a full <see cref="PipelineGrid"/> for hyperparameter search.
    /// When set, manually added pipelines via WithPipeline are ignored.
    /// </summary>
    public SupervisedExperimentBuilder WithGrid(PipelineGrid grid)
    {
        _grid = grid;
        return this;
    }

    // ── Cross-validators ─────────────────────────────────────────

    /// <summary>
    /// Add a single cross-validator strategy (the first one added is the primary for ranking).
    /// </summary>
    public SupervisedExperimentBuilder WithCrossValidator(CrossValidatorConfig config)
    {
        _cvConfigs.Add(config);
        return this;
    }

    /// <summary>Add multiple cross-validator strategies.</summary>
    public SupervisedExperimentBuilder WithCrossValidators(params CrossValidatorConfig[] configs)
    {
        _cvConfigs.AddRange(configs);
        return this;
    }

    // ── Run ──────────────────────────────────────────────────────

    /// <summary>Execute the experiment and return ranked results.</summary>
    public SupervisedExperimentResult Run()
    {
        if (_cvConfigs.Count == 0)
            throw new InvalidOperationException(
                "At least one cross-validator is required. Call WithCrossValidator() before Run().");

        var totalSw = Stopwatch.StartNew();

        // Build pipeline list
        var pipelines = _grid != null
            ? _grid.Expand().ToList()
            : new List<Pipeline>(_pipelines);

        if (pipelines.Count == 0)
            throw new InvalidOperationException(
                "At least one pipeline is required. Call WithPipeline() or WithGrid() before Run().");

        // ── Run each cross-validator with all pipelines ──────────

        var cvResults = new Dictionary<string, CrossValidationResult>();
        var cvTimings = new Dictionary<string, TimeSpan>();

        foreach (var config in _cvConfigs)
        {
            var cv = config.Create(pipelines);
            var sw = Stopwatch.StartNew();
            cvResults[config.Name] = cv.Run(_X, _y);
            sw.Stop();
            cvTimings[config.Name] = sw.Elapsed;
        }

        // ── Build per-pipeline results ───────────────────────────

        var results = new List<SupervisedResult>();

        foreach (var pipeline in pipelines)
        {
            var scores = new Dictionary<string, double>();
            var detailedResults = new Dictionary<string, CrossValidationResult>();

            foreach (var config in _cvConfigs)
            {
                var cvResult = cvResults[config.Name];
                if (cvResult.Scores.TryGetValue(pipeline, out double score))
                {
                    scores[config.Name] = score;
                    detailedResults[config.Name] = cvResult;
                }
            }

            results.Add(new SupervisedResult
            {
                Pipeline = pipeline,
                ModelName = pipeline.Model.GetType().Name,
                Scores = scores,
                Parameters = new Dictionary<string, object>(pipeline.ModelParams ?? new()),
                PipelineDescription = DescribePipeline(pipeline),
                DetailedResults = detailedResults,
                Duration = TimeSpan.FromTicks(cvTimings.Values.Sum(t => t.Ticks))
            });
        }

        // ── Rank by primary cross-validator (first one added) ────

        string primaryName = _cvConfigs[0].Name;
        var ranked = results
            .OrderByDescending(r => r.Scores.ContainsKey(primaryName)
                ? r.Scores[primaryName]
                : double.NegativeInfinity)
            .ToList();

        totalSw.Stop();

        return new SupervisedExperimentResult
        {
            Rankings = ranked,
            PrimaryCrossValidator = primaryName,
            CVResults = cvResults,
            TotalDuration = totalSw.Elapsed
        };
    }

    // ═══════════════════════════════════════════════════════════════
    //  Private: describe a pipeline for display
    // ═══════════════════════════════════════════════════════════════

    private static string DescribePipeline(Pipeline pipeline)
    {
        var name = pipeline.Model.GetType().Name;

        if (pipeline.Scaler != null)
            name = $"{pipeline.Scaler.GetType().Name} → {name}";

        if (pipeline.Selector != null)
            name = $"{pipeline.Selector.GetType().Name} → {name}";

        if (pipeline.ModelParams != null && pipeline.ModelParams.Count > 0)
        {
            var parts = new List<string>();
            foreach (var kvp in pipeline.ModelParams)
                parts.Add($"{kvp.Key}={kvp.Value}");
            name += $" ({string.Join(", ", parts)})";
        }

        return name;
    }
}
