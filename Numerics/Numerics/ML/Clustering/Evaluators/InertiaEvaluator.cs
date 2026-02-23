using CSharpNumerics.ML.Clustering.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.Clustering.Evaluators;

/// <summary>
/// Inertia (within-cluster sum of squares): sum of squared distances
/// from each point to its cluster centroid.
///
/// W = Σ_k Σ_{x ∈ Ck} ‖x − μk‖²
///
/// Lower is natively better → returned as -Inertia so higher = better
/// per IClusteringEvaluator convention.
/// Used for elbow method analysis.
/// </summary>
public class InertiaEvaluator : IClusteringEvaluator
{
    public string Name => "Inertia";

    public double Score(Matrix X, VectorN labels)
    {
        int n = X.rowLength;
        int d = X.columnLength;

        // Compute centroids
        var clusterIds = new HashSet<int>();
        for (int i = 0; i < n; i++)
        {
            int c = (int)labels[i];
            if (c >= 0) clusterIds.Add(c);
        }

        var centroids = new Dictionary<int, double[]>();
        var counts = new Dictionary<int, int>();

        foreach (int c in clusterIds)
        {
            centroids[c] = new double[d];
            counts[c] = 0;
        }

        for (int i = 0; i < n; i++)
        {
            int c = (int)labels[i];
            if (c < 0) continue;
            counts[c]++;
            for (int j = 0; j < d; j++)
                centroids[c][j] += X.values[i, j];
        }

        foreach (int c in clusterIds)
        {
            if (counts[c] > 0)
                for (int j = 0; j < d; j++)
                    centroids[c][j] /= counts[c];
        }

        // Sum of squared distances
        double inertia = 0;
        for (int i = 0; i < n; i++)
        {
            int c = (int)labels[i];
            if (c < 0) continue;

            for (int j = 0; j < d; j++)
            {
                double diff = X.values[i, j] - centroids[c][j];
                inertia += diff * diff;
            }
        }

        // Negate: lower inertia is better, but convention is higher = better
        return -inertia;
    }

    /// <summary>
    /// Returns the raw (non-negated) inertia value for elbow curve plotting.
    /// </summary>
    public double RawInertia(Matrix X, VectorN labels)
    {
        return -Score(X, labels);
    }
}
