namespace CSharpNumerics.ML.Clustering.Algorithms;

/// <summary>
/// Initialization strategy for KMeans centroid seeding.
/// </summary>
public enum KMeansInit
{
    /// <summary>Forgy method — pick K random data points as initial centroids.</summary>
    Random,

    /// <summary>KMeans++ — D²-weighted probabilistic seeding for better convergence.</summary>
    PlusPlus
}
