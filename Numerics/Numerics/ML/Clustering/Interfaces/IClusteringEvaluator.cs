
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.ML.Clustering.Interfaces;

/// <summary>
/// Scores the quality of a clustering result.
/// Convention: higher is always better.
/// Evaluators where lower is natively better (e.g. inertia, Davies-Bouldin) negate internally.
/// </summary>
public interface IClusteringEvaluator
{
    /// <summary>Display name used in result dictionaries.</summary>
    string Name { get; }

    /// <summary>
    /// Evaluate clustering quality. Higher is always better.
    /// </summary>
    /// <param name="X">The data matrix (n × d).</param>
    /// <param name="labels">Cluster assignments for each row (length n).</param>
    double Score(Matrix X, VectorN labels);
}
