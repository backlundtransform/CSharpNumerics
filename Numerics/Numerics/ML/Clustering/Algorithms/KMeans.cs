using CSharpNumerics.ML.Clustering.Interfaces;
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.Clustering.Algorithms;

/// <summary>
/// KMeans clustering with optional KMeans++ initialization.
/// Lloyd's algorithm: assign → recompute centroids → repeat until convergence.
/// </summary>
public class KMeans : IClusteringModel, IHasHyperparameters
{
    // ── Hyperparameters ──────────────────────────────────────────
    public int K { get; set; } = 3;
    public int MaxIterations { get; set; } = 300;
    public double Tolerance { get; set; } = 1e-4;
    public int? Seed { get; set; }
    public KMeansInit InitMethod { get; set; } = KMeansInit.PlusPlus;

    // ── Results (available after Fit) ────────────────────────────
    public Matrix Centroids { get; private set; }
    public double Inertia { get; private set; }
    public int Iterations { get; private set; }
    public int ClusterCount => K;

    private bool _isFitted;

    // ── IHasHyperparameters ──────────────────────────────────────
    public void SetHyperParameters(Dictionary<string, object> parameters)
    {
        if (parameters == null) return;

        if (parameters.ContainsKey("K"))
            K = Convert.ToInt32(parameters["K"]);
        if (parameters.ContainsKey("MaxIterations"))
            MaxIterations = Convert.ToInt32(parameters["MaxIterations"]);
        if (parameters.ContainsKey("Tolerance"))
            Tolerance = Convert.ToDouble(parameters["Tolerance"]);
        if (parameters.ContainsKey("Seed"))
            Seed = Convert.ToInt32(parameters["Seed"]);
        if (parameters.ContainsKey("InitMethod"))
            InitMethod = (KMeansInit)parameters["InitMethod"];
    }

    // ── Fit ──────────────────────────────────────────────────────
    public void Fit(Matrix X)
    {
        int n = X.rowLength;
        int d = X.columnLength;

        if (n < K)
            throw new ArgumentException($"Not enough samples ({n}) for {K} clusters.");

        var rng = Seed.HasValue ? new Random(Seed.Value) : new Random();

        // Initialize centroids
        double[,] centroids = InitMethod == KMeansInit.PlusPlus
            ? InitPlusPlus(X, rng)
            : InitRandom(X, rng);

        int[] assignments = new int[n];
        double prevInertia = double.MaxValue;

        for (int iter = 0; iter < MaxIterations; iter++)
        {
            // ── Assign step ──
            double inertia = 0;
            for (int i = 0; i < n; i++)
            {
                double bestDist = double.MaxValue;
                int bestCluster = 0;

                for (int c = 0; c < K; c++)
                {
                    double dist = SquaredEuclidean(X.values, i, centroids, c, d);
                    if (dist < bestDist)
                    {
                        bestDist = dist;
                        bestCluster = c;
                    }
                }

                assignments[i] = bestCluster;
                inertia += bestDist;
            }

            // ── Update step ──
            var newCentroids = new double[K, d];
            var counts = new int[K];

            for (int i = 0; i < n; i++)
            {
                int c = assignments[i];
                counts[c]++;
                for (int j = 0; j < d; j++)
                    newCentroids[c, j] += X.values[i, j];
            }

            for (int c = 0; c < K; c++)
            {
                if (counts[c] == 0)
                {
                    // Reinitialize empty cluster to a random data point
                    int idx = rng.Next(n);
                    for (int j = 0; j < d; j++)
                        newCentroids[c, j] = X.values[idx, j];
                }
                else
                {
                    for (int j = 0; j < d; j++)
                        newCentroids[c, j] /= counts[c];
                }
            }

            // ── Convergence check ──
            double shift = CentroidShift(centroids, newCentroids, K, d);
            centroids = newCentroids;
            Iterations = iter + 1;
            Inertia = inertia;

            if (shift < Tolerance)
                break;

            prevInertia = inertia;
        }

        Centroids = new Matrix(centroids);
        _isFitted = true;
    }

    // ── Predict ──────────────────────────────────────────────────
    public VectorN Predict(Matrix X)
    {
        if (!_isFitted)
            throw new InvalidOperationException("KMeans has not been fitted.");

        int n = X.rowLength;
        int d = X.columnLength;
        var labels = new VectorN(n);

        for (int i = 0; i < n; i++)
        {
            double bestDist = double.MaxValue;
            int bestCluster = 0;

            for (int c = 0; c < K; c++)
            {
                double dist = SquaredEuclidean(X.values, i, Centroids.values, c, d);
                if (dist < bestDist)
                {
                    bestDist = dist;
                    bestCluster = c;
                }
            }

            labels[i] = bestCluster;
        }

        return labels;
    }

    // ── FitPredict ───────────────────────────────────────────────
    public VectorN FitPredict(Matrix X)
    {
        Fit(X);
        return Predict(X);
    }

    // ── Clone ────────────────────────────────────────────────────
    public IClusteringModel Clone()
    {
        return new KMeans
        {
            K = K,
            MaxIterations = MaxIterations,
            Tolerance = Tolerance,
            Seed = Seed,
            InitMethod = InitMethod
        };
    }

    // ═══════════════════════════════════════════════════════════════
    //  Private helpers
    // ═══════════════════════════════════════════════════════════════

    /// <summary>Forgy initialization: pick K random distinct data points.</summary>
    private double[,] InitRandom(Matrix X, Random rng)
    {
        int n = X.rowLength;
        int d = X.columnLength;
        var centroids = new double[K, d];
        var chosen = new HashSet<int>();

        while (chosen.Count < K)
            chosen.Add(rng.Next(n));

        int c = 0;
        foreach (int idx in chosen)
        {
            for (int j = 0; j < d; j++)
                centroids[c, j] = X.values[idx, j];
            c++;
        }

        return centroids;
    }

    /// <summary>
    /// KMeans++ initialization: D²-weighted probabilistic seeding.
    /// First centroid is random; each subsequent centroid is chosen with
    /// probability proportional to squared distance from nearest existing centroid.
    /// </summary>
    private double[,] InitPlusPlus(Matrix X, Random rng)
    {
        int n = X.rowLength;
        int d = X.columnLength;
        var centroids = new double[K, d];

        // Pick first centroid at random
        int first = rng.Next(n);
        for (int j = 0; j < d; j++)
            centroids[0, j] = X.values[first, j];

        var minDistSq = new double[n];

        for (int c = 1; c < K; c++)
        {
            // Compute min squared distance from each point to nearest chosen centroid
            double totalWeight = 0;
            for (int i = 0; i < n; i++)
            {
                double dist = SquaredEuclidean(X.values, i, centroids, c - 1, d);
                if (c == 1)
                    minDistSq[i] = dist;
                else
                    minDistSq[i] = Math.Min(minDistSq[i], dist);

                totalWeight += minDistSq[i];
            }

            // Weighted random selection
            double r = rng.NextDouble() * totalWeight;
            double cumulative = 0;
            int selected = n - 1; // fallback

            for (int i = 0; i < n; i++)
            {
                cumulative += minDistSq[i];
                if (cumulative >= r)
                {
                    selected = i;
                    break;
                }
            }

            for (int j = 0; j < d; j++)
                centroids[c, j] = X.values[selected, j];
        }

        return centroids;
    }

    /// <summary>Squared Euclidean distance between row i of A and row j of B.</summary>
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

    /// <summary>Max squared shift of any centroid between iterations.</summary>
    private static double CentroidShift(double[,] old, double[,] next, int k, int d)
    {
        double maxShift = 0;
        for (int c = 0; c < k; c++)
        {
            double shift = 0;
            for (int j = 0; j < d; j++)
            {
                double diff = old[c, j] - next[c, j];
                shift += diff * diff;
            }
            if (shift > maxShift)
                maxShift = shift;
        }
        return maxShift;
    }
}
