namespace CSharpNumerics.ML.Clustering.Algorithms;

/// <summary>
/// Linkage strategy for Agglomerative Clustering.
/// </summary>
public enum LinkageType
{
    /// <summary>Minimum distance between any two points in the two clusters.</summary>
    Single,

    /// <summary>Maximum distance between any two points in the two clusters.</summary>
    Complete,

    /// <summary>Average distance between all pairs of points across the two clusters.</summary>
    Average,

    /// <summary>Minimizes total within-cluster variance (Lance-Williams formula).</summary>
    Ward
}
