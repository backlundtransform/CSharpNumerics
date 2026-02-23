using CSharpNumerics.ML.Clustering.Interfaces;
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.Clustering.Algorithms;

/// <summary>
/// Density-Based Spatial Clustering of Applications with Noise.
/// Discovers clusters of arbitrary shape without requiring K upfront.
/// Noise points are assigned label -1.
/// </summary>
public class DBSCAN : IClusteringModel, IHasHyperparameters
{
    // ── Hyperparameters ──────────────────────────────────────────
    public double Epsilon { get; set; } = 0.5;
    public int MinPoints { get; set; } = 5;

    // ── Results (available after Fit) ────────────────────────────
    public int ClusterCount { get; private set; }
    public int NoiseCount { get; private set; }

    private VectorN _labels;
    private Matrix _fitData;
    private bool _isFitted;

    // ── IHasHyperparameters ──────────────────────────────────────
    public void SetHyperParameters(Dictionary<string, object> parameters)
    {
        if (parameters == null) return;

        if (parameters.ContainsKey("Epsilon"))
            Epsilon = Convert.ToDouble(parameters["Epsilon"]);
        if (parameters.ContainsKey("MinPoints"))
            MinPoints = Convert.ToInt32(parameters["MinPoints"]);
    }

    // ── Fit ──────────────────────────────────────────────────────
    public void Fit(Matrix X)
    {
        int n = X.rowLength;
        int d = X.columnLength;
        double epsSq = Epsilon * Epsilon;

        // -2 = unvisited, -1 = noise, 0+ = cluster id
        int[] labels = new int[n];
        for (int i = 0; i < n; i++)
            labels[i] = -2;

        int clusterId = 0;

        for (int i = 0; i < n; i++)
        {
            if (labels[i] != -2)
                continue; // already processed

            var neighbors = RangeQuery(X, i, epsSq, n, d);

            if (neighbors.Count < MinPoints)
            {
                labels[i] = -1; // noise
                continue;
            }

            // Start a new cluster
            labels[i] = clusterId;
            var seeds = new Queue<int>(neighbors);

            while (seeds.Count > 0)
            {
                int q = seeds.Dequeue();

                if (labels[q] == -1)
                {
                    // Was noise, absorb into cluster
                    labels[q] = clusterId;
                    continue;
                }

                if (labels[q] != -2)
                    continue; // already assigned to a cluster

                labels[q] = clusterId;

                var qNeighbors = RangeQuery(X, q, epsSq, n, d);

                if (qNeighbors.Count >= MinPoints)
                {
                    foreach (int nb in qNeighbors)
                        seeds.Enqueue(nb);
                }
            }

            clusterId++;
        }

        // Store results
        _labels = new VectorN(n);
        int noiseCount = 0;
        for (int i = 0; i < n; i++)
        {
            _labels[i] = labels[i];
            if (labels[i] == -1) noiseCount++;
        }

        ClusterCount = clusterId;
        NoiseCount = noiseCount;
        _fitData = X;
        _isFitted = true;
    }

    // ── Predict ──────────────────────────────────────────────────
    /// <summary>
    /// Assign new points to the nearest cluster found during Fit.
    /// Points farther than Epsilon from all core points are labeled -1 (noise).
    /// </summary>
    public VectorN Predict(Matrix X)
    {
        if (!_isFitted)
            throw new InvalidOperationException("DBSCAN has not been fitted.");

        int n = X.rowLength;
        int d = X.columnLength;
        double epsSq = Epsilon * Epsilon;
        var labels = new VectorN(n);

        for (int i = 0; i < n; i++)
        {
            double bestDist = double.MaxValue;
            int bestLabel = -1;

            for (int j = 0; j < _fitData.rowLength; j++)
            {
                if (_labels[j] < 0) continue; // skip noise

                double dist = SquaredEuclidean(X.values, i, _fitData.values, j, d);
                if (dist < bestDist)
                {
                    bestDist = dist;
                    bestLabel = (int)_labels[j];
                }
            }

            labels[i] = bestDist <= epsSq ? bestLabel : -1;
        }

        return labels;
    }

    // ── FitPredict ───────────────────────────────────────────────
    public VectorN FitPredict(Matrix X)
    {
        Fit(X);
        return new VectorN(_labels.Values);
    }

    // ── Clone ────────────────────────────────────────────────────
    public IClusteringModel Clone()
    {
        return new DBSCAN
        {
            Epsilon = Epsilon,
            MinPoints = MinPoints
        };
    }

    // ═══════════════════════════════════════════════════════════════
    //  Private helpers
    // ═══════════════════════════════════════════════════════════════

    /// <summary>Find all point indices within Epsilon of point p.</summary>
    private List<int> RangeQuery(Matrix X, int p, double epsSq, int n, int d)
    {
        var result = new List<int>();
        for (int i = 0; i < n; i++)
        {
            if (i == p) continue;
            if (SquaredEuclidean(X.values, p, X.values, i, d) <= epsSq)
                result.Add(i);
        }
        return result;
    }

    private static double SquaredEuclidean(double[,] A, int i, double[,] B, int j, int d)
    {
        double sum = 0;
        for (int k = 0; k < d; k++)
        {
            double diff = A[i, k] - B[j, k];
            sum += diff * diff;
        }
        return sum;
    }
}
