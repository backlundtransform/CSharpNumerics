using CSharpNumerics.ML.Clustering.Interfaces;
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.Clustering.Algorithms;

/// <summary>
/// Agglomerative (bottom-up) hierarchical clustering.
/// Each point starts as its own cluster; the two closest clusters are merged
/// repeatedly until the desired number of clusters K is reached.
/// </summary>
public class AgglomerativeClustering : IClusteringModel, IHasHyperparameters
{
    // ── Hyperparameters ──────────────────────────────────────────
    public int K { get; set; } = 3;
    public LinkageType Linkage { get; set; } = LinkageType.Ward;

    // ── Results ──────────────────────────────────────────────────
    public int ClusterCount => K;

    /// <summary>Merge history: (clusterA, clusterB, distance) in merge order.</summary>
    public List<(int A, int B, double Distance)> Dendrogram { get; private set; }

    private int[] _labels;
    private Matrix _fitData;
    private bool _isFitted;

    // ── IHasHyperparameters ──────────────────────────────────────
    public void SetHyperParameters(Dictionary<string, object> parameters)
    {
        if (parameters == null) return;

        if (parameters.ContainsKey("K"))
            K = Convert.ToInt32(parameters["K"]);
        if (parameters.ContainsKey("Linkage"))
            Linkage = (LinkageType)parameters["Linkage"];
    }

    // ── Fit ──────────────────────────────────────────────────────
    public void Fit(Matrix X)
    {
        int n = X.rowLength;
        int d = X.columnLength;

        if (n < K)
            throw new ArgumentException($"Not enough samples ({n}) for {K} clusters.");

        // Each point starts as its own cluster
        // clusters[i] = list of point indices belonging to cluster i
        var clusters = new Dictionary<int, List<int>>();
        for (int i = 0; i < n; i++)
            clusters[i] = new List<int> { i };

        // Precompute full pairwise distance matrix (Euclidean)
        double[,] distMatrix = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
            {
                double dist = EuclideanDistance(X.values, i, j, d);
                distMatrix[i, j] = dist;
                distMatrix[j, i] = dist;
            }

        // Linkage distance cache between active clusters
        // Key: (min(id), max(id)) → distance
        var activeIds = new HashSet<int>(Enumerable.Range(0, n));
        int nextId = n; // for newly merged clusters

        Dendrogram = new List<(int, int, double)>();

        // Cluster sizes (needed for Ward)
        var sizes = new Dictionary<int, int>();
        for (int i = 0; i < n; i++)
            sizes[i] = 1;

        // Cluster centroids (needed for Ward)
        var centroids = new Dictionary<int, double[]>();
        for (int i = 0; i < n; i++)
        {
            var c = new double[d];
            for (int j = 0; j < d; j++)
                c[j] = X.values[i, j];
            centroids[i] = c;
        }

        // Merge until K clusters remain
        while (activeIds.Count > K)
        {
            // Find closest pair of active clusters
            double bestDist = double.MaxValue;
            int bestA = -1, bestB = -1;

            var ids = activeIds.ToArray();
            for (int ai = 0; ai < ids.Length; ai++)
            {
                for (int bi = ai + 1; bi < ids.Length; bi++)
                {
                    double dist = LinkageDistance(
                        ids[ai], ids[bi], clusters, distMatrix, sizes, centroids, d);

                    if (dist < bestDist)
                    {
                        bestDist = dist;
                        bestA = ids[ai];
                        bestB = ids[bi];
                    }
                }
            }

            // Merge bestA and bestB into a new cluster
            var merged = new List<int>(clusters[bestA]);
            merged.AddRange(clusters[bestB]);

            int mergedId = nextId++;
            clusters[mergedId] = merged;
            sizes[mergedId] = sizes[bestA] + sizes[bestB];

            // Update centroid for merged cluster
            int sA = sizes[bestA], sB = sizes[bestB], sM = sizes[mergedId];
            var cA = centroids[bestA];
            var cB = centroids[bestB];
            var cM = new double[d];
            for (int j = 0; j < d; j++)
                cM[j] = (cA[j] * sA + cB[j] * sB) / sM;
            centroids[mergedId] = cM;

            // Expand distance matrix if needed — for point-level lookups
            // (no need; we keep the original n×n matrix and compute linkage from point lists)

            Dendrogram.Add((bestA, bestB, bestDist));

            activeIds.Remove(bestA);
            activeIds.Remove(bestB);
            activeIds.Add(mergedId);

            clusters.Remove(bestA);
            clusters.Remove(bestB);
        }

        // Assign final labels
        _labels = new int[n];
        int label = 0;
        foreach (var clusterId in activeIds)
        {
            foreach (int idx in clusters[clusterId])
                _labels[idx] = label;
            label++;
        }

        _fitData = X;
        _isFitted = true;
    }

    // ── Predict ──────────────────────────────────────────────────
    /// <summary>
    /// Assign new points to the nearest cluster centroid found during Fit.
    /// </summary>
    public VectorN Predict(Matrix X)
    {
        if (!_isFitted)
            throw new InvalidOperationException("AgglomerativeClustering has not been fitted.");

        // Compute cluster centroids from fit data
        int d = _fitData.columnLength;
        var centroidMap = ComputeCentroids(_fitData, _labels, K, d);

        int n = X.rowLength;
        var labels = new VectorN(n);

        for (int i = 0; i < n; i++)
        {
            double bestDist = double.MaxValue;
            int bestLabel = 0;

            for (int c = 0; c < K; c++)
            {
                double dist = 0;
                for (int j = 0; j < d; j++)
                {
                    double diff = X.values[i, j] - centroidMap[c, j];
                    dist += diff * diff;
                }
                if (dist < bestDist)
                {
                    bestDist = dist;
                    bestLabel = c;
                }
            }

            labels[i] = bestLabel;
        }

        return labels;
    }

    // ── FitPredict ───────────────────────────────────────────────
    public VectorN FitPredict(Matrix X)
    {
        Fit(X);
        var result = new VectorN(X.rowLength);
        for (int i = 0; i < X.rowLength; i++)
            result[i] = _labels[i];
        return result;
    }

    // ── Clone ────────────────────────────────────────────────────
    public IClusteringModel Clone()
    {
        return new AgglomerativeClustering
        {
            K = K,
            Linkage = Linkage
        };
    }

    // ═══════════════════════════════════════════════════════════════
    //  Private helpers
    // ═══════════════════════════════════════════════════════════════

    private double LinkageDistance(
        int idA, int idB,
        Dictionary<int, List<int>> clusters,
        double[,] distMatrix,
        Dictionary<int, int> sizes,
        Dictionary<int, double[]> centroids,
        int d)
    {
        switch (Linkage)
        {
            case LinkageType.Single:
                return SingleLinkage(clusters[idA], clusters[idB], distMatrix);
            case LinkageType.Complete:
                return CompleteLinkage(clusters[idA], clusters[idB], distMatrix);
            case LinkageType.Average:
                return AverageLinkage(clusters[idA], clusters[idB], distMatrix);
            case LinkageType.Ward:
                return WardDistance(centroids[idA], centroids[idB], sizes[idA], sizes[idB], d);
            default:
                throw new ArgumentException($"Unknown linkage type: {Linkage}");
        }
    }

    private static double SingleLinkage(List<int> a, List<int> b, double[,] dist)
    {
        double min = double.MaxValue;
        foreach (int i in a)
            foreach (int j in b)
                if (dist[i, j] < min) min = dist[i, j];
        return min;
    }

    private static double CompleteLinkage(List<int> a, List<int> b, double[,] dist)
    {
        double max = double.MinValue;
        foreach (int i in a)
            foreach (int j in b)
                if (dist[i, j] > max) max = dist[i, j];
        return max;
    }

    private static double AverageLinkage(List<int> a, List<int> b, double[,] dist)
    {
        double sum = 0;
        foreach (int i in a)
            foreach (int j in b)
                sum += dist[i, j];
        return sum / (a.Count * b.Count);
    }

    /// <summary>
    /// Ward distance: increase in total within-cluster variance when merging.
    /// Δ = (nA * nB) / (nA + nB) * ‖cA - cB‖²
    /// </summary>
    private static double WardDistance(double[] cA, double[] cB, int nA, int nB, int d)
    {
        double distSq = 0;
        for (int j = 0; j < d; j++)
        {
            double diff = cA[j] - cB[j];
            distSq += diff * diff;
        }
        return (double)(nA * nB) / (nA + nB) * distSq;
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

    private static double[,] ComputeCentroids(Matrix X, int[] labels, int k, int d)
    {
        var centroids = new double[k, d];
        var counts = new int[k];

        for (int i = 0; i < X.rowLength; i++)
        {
            int c = labels[i];
            counts[c]++;
            for (int j = 0; j < d; j++)
                centroids[c, j] += X.values[i, j];
        }

        for (int c = 0; c < k; c++)
        {
            if (counts[c] > 0)
                for (int j = 0; j < d; j++)
                    centroids[c, j] /= counts[c];
        }

        return centroids;
    }
}
