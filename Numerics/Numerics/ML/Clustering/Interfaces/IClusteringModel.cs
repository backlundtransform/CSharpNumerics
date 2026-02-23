
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.ML.Clustering.Interfaces;

/// <summary>
/// Core contract for unsupervised clustering algorithms.
/// Parallel to <see cref="Models.Interfaces.IModel"/> but without a target vector.
/// </summary>
public interface IClusteringModel
{
    /// <summary>Fit the model to unlabeled data.</summary>
    void Fit(Matrix X);

    /// <summary>Assign cluster labels to rows of X.</summary>
    VectorN Predict(Matrix X);

    /// <summary>Fit and predict in one step.</summary>
    VectorN FitPredict(Matrix X);

    /// <summary>Number of clusters found or configured.</summary>
    int ClusterCount { get; }

    /// <summary>Deep copy for grid search isolation.</summary>
    IClusteringModel Clone();
}
