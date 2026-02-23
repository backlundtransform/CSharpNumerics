using CSharpNumerics.ML.Clustering.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.Clustering.Evaluators;

/// <summary>
/// Calinski-Harabasz Index (Variance Ratio Criterion):
/// ratio of between-cluster dispersion to within-cluster dispersion.
///
/// CH = (SS_B / (K - 1)) / (SS_W / (n - K))
///
/// Higher is better (no negation needed).
/// Fast to compute — only requires centroids.
/// </summary>
public class CalinskiHarabaszEvaluator : IClusteringEvaluator
{
    public string Name => "CalinskiHarabasz";

    public double Score(Matrix X, VectorN labels)
    {
        int n = X.rowLength;
        int d = X.columnLength;

        // Identify clusters (skip noise label -1)
        var clusterIds = new List<int>();
        var clusterIdSet = new HashSet<int>();
        int validN = 0;

        for (int i = 0; i < n; i++)
        {
            int c = (int)labels[i];
            if (c >= 0)
            {
                validN++;
                if (clusterIdSet.Add(c))
                    clusterIds.Add(c);
            }
        }

        int K = clusterIds.Count;
        if (K < 2 || validN <= K) return 0;

        // Global centroid
        var globalCentroid = new double[d];
        for (int i = 0; i < n; i++)
        {
            if ((int)labels[i] < 0) continue;
            for (int j = 0; j < d; j++)
                globalCentroid[j] += X.values[i, j];
        }
        for (int j = 0; j < d; j++)
            globalCentroid[j] /= validN;

        // Cluster centroids and counts
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

        // SS_B: between-cluster sum of squares
        // SS_B = Σ_k n_k * ‖c_k − c_global‖²
        double ssB = 0;
        foreach (int c in clusterIds)
        {
            double distSq = 0;
            for (int j = 0; j < d; j++)
            {
                double diff = centroids[c][j] - globalCentroid[j];
                distSq += diff * diff;
            }
            ssB += counts[c] * distSq;
        }

        // SS_W: within-cluster sum of squares
        // SS_W = Σ_k Σ_{x ∈ Ck} ‖x − c_k‖²
        double ssW = 0;
        for (int i = 0; i < n; i++)
        {
            int c = (int)labels[i];
            if (c < 0) continue;

            for (int j = 0; j < d; j++)
            {
                double diff = X.values[i, j] - centroids[c][j];
                ssW += diff * diff;
            }
        }

        if (ssW == 0) return double.MaxValue;

        // CH = (SS_B / (K-1)) / (SS_W / (n-K))
        return (ssB / (K - 1)) / (ssW / (validN - K));
    }
}
