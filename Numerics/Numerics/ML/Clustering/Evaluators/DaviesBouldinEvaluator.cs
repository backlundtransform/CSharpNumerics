using CSharpNumerics.ML.Clustering.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.Clustering.Evaluators;

/// <summary>
/// Davies-Bouldin Index: average similarity ratio between each cluster
/// and its most similar cluster.
///
/// DB = (1/K) Σ_i max_{j≠i} (σi + σj) / d(ci, cj)
///
/// where σi = average distance of points in cluster i to centroid,
/// d(ci, cj) = distance between centroids.
///
/// Lower is natively better → returned as -DB so higher = better
/// per IClusteringEvaluator convention.
/// </summary>
public class DaviesBouldinEvaluator : IClusteringEvaluator
{
    public string Name => "DaviesBouldin";

    public double Score(Matrix X, VectorN labels)
    {
        int n = X.rowLength;
        int d = X.columnLength;

        // Identify clusters (skip noise label -1)
        var clusterIds = new List<int>();
        var clusterIdSet = new HashSet<int>();
        for (int i = 0; i < n; i++)
        {
            int c = (int)labels[i];
            if (c >= 0 && clusterIdSet.Add(c))
                clusterIds.Add(c);
        }

        int K = clusterIds.Count;
        if (K < 2) return double.NegativeInfinity;

        // Compute centroids
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

        // σ_i = average distance from points to their centroid
        var scatter = new Dictionary<int, double>();
        foreach (int c in clusterIds)
            scatter[c] = 0;

        for (int i = 0; i < n; i++)
        {
            int c = (int)labels[i];
            if (c < 0) continue;

            double dist = 0;
            for (int j = 0; j < d; j++)
            {
                double diff = X.values[i, j] - centroids[c][j];
                dist += diff * diff;
            }
            scatter[c] += Math.Sqrt(dist);
        }

        foreach (int c in clusterIds)
        {
            if (counts[c] > 0)
                scatter[c] /= counts[c];
        }

        // DB = (1/K) Σ_i max_{j≠i} (σi + σj) / d(ci,cj)
        double db = 0;
        for (int ai = 0; ai < K; ai++)
        {
            int a = clusterIds[ai];
            double maxRatio = double.NegativeInfinity;

            for (int bi = 0; bi < K; bi++)
            {
                if (ai == bi) continue;
                int b = clusterIds[bi];

                double centroidDist = 0;
                for (int j = 0; j < d; j++)
                {
                    double diff = centroids[a][j] - centroids[b][j];
                    centroidDist += diff * diff;
                }
                centroidDist = Math.Sqrt(centroidDist);

                if (centroidDist == 0) continue;

                double ratio = (scatter[a] + scatter[b]) / centroidDist;
                if (ratio > maxRatio)
                    maxRatio = ratio;
            }

            if (maxRatio > double.NegativeInfinity)
                db += maxRatio;
        }

        db /= K;

        // Negate: lower DB is better, but convention is higher = better
        return -db;
    }
}
