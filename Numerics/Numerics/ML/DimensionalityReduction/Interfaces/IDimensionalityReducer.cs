using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.ML.DimensionalityReduction.Interfaces;

/// <summary>
/// Core contract for unsupervised dimensionality reduction algorithms.
/// Transforms high-dimensional data into a lower-dimensional representation.
/// Parallel to <see cref="Scalers.Interfaces.IScaler"/> but changes the number of columns.
/// </summary>
public interface IDimensionalityReducer
{
    /// <summary>Number of output dimensions after transformation.</summary>
    int NComponents { get; set; }

    /// <summary>Fit to data and return the transformed (reduced) matrix.</summary>
    Matrix FitTransform(Matrix X);

    /// <summary>Transform new data using the already-fitted reducer.</summary>
    Matrix Transform(Matrix X);

    /// <summary>Deep copy for grid search isolation.</summary>
    IDimensionalityReducer Clone();
}
