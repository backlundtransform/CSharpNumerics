using CSharpNumerics.ML.CrossValidators.Result;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.Experiment.Results;

/// <summary>
/// Result from a single pipeline evaluated across one or more cross-validators.
/// Analogous to <see cref="CSharpNumerics.ML.Clustering.Results.ClusteringResult"/>
/// but for supervised models.
/// </summary>
public class SupervisedResult
{
    /// <summary>The pipeline configuration that was evaluated.</summary>
    public Pipeline Pipeline { get; set; }

    /// <summary>Model class name (e.g. "KNearestNeighbors", "DecisionTree").</summary>
    public string ModelName { get; set; }

    /// <summary>
    /// Cross-validator name → score.
    /// Convention: higher is always better (accuracy for classification, -MSE for regression).
    /// </summary>
    public Dictionary<string, double> Scores { get; set; } = new();

    /// <summary>Hyperparameters used for this pipeline.</summary>
    public Dictionary<string, object> Parameters { get; set; } = new();

    /// <summary>Human-readable pipeline description (model + scaler + params).</summary>
    public string PipelineDescription { get; set; }

    /// <summary>Total wall-clock time across all cross-validators for this pipeline.</summary>
    public TimeSpan Duration { get; set; }

    /// <summary>
    /// Detailed <see cref="CrossValidationResult"/> per cross-validator.
    /// Key = CV name. Gives access to actual/predicted values, confusion matrix, etc.
    /// for this pipeline within each cross-validation strategy.
    /// </summary>
    public Dictionary<string, CrossValidationResult> DetailedResults { get; set; } = new();
}
