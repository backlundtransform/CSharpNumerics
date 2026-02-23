using CSharpNumerics.ML.Clustering.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.Clustering.Evaluators;

/// <summary>
/// Silhouette score: measures how similar each point is to its own cluster
/// compared to the nearest other cluster.
/// 
/// s(i) = (b(i) - a(i)) / max(a(i), b(i))
/// 
/// where a(i) = mean intra-cluster distance, b(i) = mean nearest-cluster distance.
/// Range: [-1, 1], higher is better.
/// </summary>
public class SilhouetteEvaluator : IClusteringEvaluator
{
    public string Name => "Silhouette";

    public double Score(Matrix X, VectorN labels)
    {
        int n = X.rowLength;
        int d = X.columnLength;

        if (n < 2) return 0;

        // Determine cluster ids and counts
        var clusterIds = new HashSet<int>();
        for (int i = 0; i < n; i++)
            clusterIds.Add((int)labels[i]);

        // Need at least 2 clusters for silhouette to be meaningful
        if (clusterIds.Count < 2) return -1;

        // Group point indices by cluster
        var clusterPoints = new Dictionary<int, List<int>>();
        foreach (int c in clusterIds)
            clusterPoints[c] = new List<int>();
        for (int i = 0; i < n; i++)
            clusterPoints[(int)labels[i]].Add(i);

        double totalSilhouette = 0;
        int validPoints = 0;

        for (int i = 0; i < n; i++)
        {
            int ci = (int)labels[i];

            // Skip noise points (DBSCAN label -1)
            if (ci < 0) continue;

            // Skip singleton clusters
            if (clusterPoints[ci].Count <= 1)
            {
                validPoints++;
                // s(i) = 0 for singleton
                continue;
            }

            // a(i): mean distance to all other points in same cluster
            double a = 0;
            int aCount = 0;
            foreach (int j in clusterPoints[ci])
            {
                if (j == i) continue;
                a += EuclideanDistance(X.values, i, j, d);
                aCount++;
            }
            a /= aCount;

            // b(i): min over other clusters of mean distance to that cluster
            double b = double.MaxValue;
            foreach (var kvp in clusterPoints)
            {
                if (kvp.Key == ci || kvp.Key < 0) continue;

                double meanDist = 0;
                foreach (int j in kvp.Value)
                    meanDist += EuclideanDistance(X.values, i, j, d);
                meanDist /= kvp.Value.Count;

                if (meanDist < b)
                    b = meanDist;
            }

            double si = (b - a) / Math.Max(a, b);
            totalSilhouette += si;
            validPoints++;
        }

        return validPoints > 0 ? totalSilhouette / validPoints : 0;
    }

    private static double EuclideanDistance(double[,] X, int i, int j, int d)
    {
        double sum = 0;
        for (int k = 0; k < d; k++)
        {
            double diff = X[i, k] - X[j, k];
            sum += diff * diff;
        }
        return Math.Sqrt(sum);
    }
}
