using CSharpNumerics.ML.Clustering.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.Clustering.Results;

/// <summary>
/// Result from a single clustering pipeline run: one algorithm, one configuration.
/// </summary>
public class ClusteringResult
{
    /// <summary>The fitted clustering model.</summary>
    public IClusteringModel Model { get; set; }

    /// <summary>Algorithm class name (e.g. "KMeans").</summary>
    public string AlgorithmName { get; set; }

    /// <summary>Number of clusters found or configured.</summary>
    public int ClusterCount { get; set; }

    /// <summary>Cluster assignments for each data point.</summary>
    public VectorN Labels { get; set; }

    /// <summary>Evaluator name → score. Higher is always better.</summary>
    public Dictionary<string, double> Scores { get; set; } = new();

    /// <summary>Hyperparameters used for this run.</summary>
    public Dictionary<string, object> Parameters { get; set; } = new();

    /// <summary>Pipeline description (model + scaler + params).</summary>
    public string PipelineDescription { get; set; }

    /// <summary>Wall-clock time for FitPredict.</summary>
    public TimeSpan Duration { get; set; }
}
