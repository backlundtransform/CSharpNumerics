using CSharpNumerics.ML.DimensionalityReduction.Interfaces;
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.DimensionalityReduction.Algorithms;

/// <summary>
/// Principal Component Analysis (PCA).
/// Projects data onto the top <see cref="NComponents"/> eigenvectors of the covariance matrix.
/// Uses the power iteration method with deflation for eigendecomposition.
/// </summary>
public class PCA : IDimensionalityReducer, IHasHyperparameters
{
    // ── Hyperparameters ──────────────────────────────────────────
    public int NComponents { get; set; } = 2;
    public int MaxIterations { get; set; } = 1000;
    public double Tolerance { get; set; } = 1e-8;
    public int? Seed { get; set; }

    // ── Results (available after Fit) ────────────────────────────

    /// <summary>Principal components (NComponents × d), one eigenvector per row.</summary>
    public double[,] Components { get; private set; }

    /// <summary>Variance explained by each component.</summary>
    public double[] ExplainedVariance { get; private set; }

    /// <summary>Fraction of total variance explained by each component (sums to ≤ 1).</summary>
    public double[] ExplainedVarianceRatio { get; private set; }

    /// <summary>Column means used for centering (length d).</summary>
    public double[] Mean { get; private set; }

    private bool _isFitted;

    // ── IHasHyperparameters ──────────────────────────────────────
    public void SetHyperParameters(Dictionary<string, object> parameters)
    {
        if (parameters == null) return;

        if (parameters.ContainsKey("NComponents"))
            NComponents = Convert.ToInt32(parameters["NComponents"]);
        if (parameters.ContainsKey("MaxIterations"))
            MaxIterations = Convert.ToInt32(parameters["MaxIterations"]);
        if (parameters.ContainsKey("Tolerance"))
            Tolerance = Convert.ToDouble(parameters["Tolerance"]);
        if (parameters.ContainsKey("Seed"))
            Seed = Convert.ToInt32(parameters["Seed"]);
    }

    // ── Fit ──────────────────────────────────────────────────────
    public void Fit(Matrix X)
    {
        int n = X.rowLength;
        int d = X.columnLength;

        if (NComponents > d)
            throw new ArgumentException(
                $"NComponents ({NComponents}) cannot exceed number of features ({d}).");
        if (NComponents < 1)
            throw new ArgumentException("NComponents must be at least 1.");

        // Center data
        Mean = ComputeMean(X.values, n, d);
        double[,] centered = CenterData(X.values, Mean, n, d);

        // Covariance matrix (d × d)
        double[,] cov = ComputeCovariance(centered, n, d);

        // Eigen decomposition via power iteration with deflation
        Components = new double[NComponents, d];
        ExplainedVariance = new double[NComponents];

        var rng = Seed.HasValue ? new Random(Seed.Value) : new Random();

        for (int k = 0; k < NComponents; k++)
        {
            var (eigenvalue, eigenvector) = PowerIteration(cov, d, rng);
            ExplainedVariance[k] = eigenvalue;

            for (int j = 0; j < d; j++)
                Components[k, j] = eigenvector[j];

            // Deflate: remove the found component from the covariance matrix
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                    cov[i, j] -= eigenvalue * eigenvector[i] * eigenvector[j];
        }

        // Compute explained variance ratio
        double totalVariance = 0;
        double[,] fullCov = ComputeCovariance(centered, n, d);
        for (int i = 0; i < d; i++)
            totalVariance += fullCov[i, i];

        ExplainedVarianceRatio = new double[NComponents];
        for (int k = 0; k < NComponents; k++)
            ExplainedVarianceRatio[k] = totalVariance > 0
                ? ExplainedVariance[k] / totalVariance
                : 0;

        _isFitted = true;
    }

    // ── FitTransform ─────────────────────────────────────────────
    public Matrix FitTransform(Matrix X)
    {
        Fit(X);
        return Transform(X);
    }

    // ── Transform ────────────────────────────────────────────────
    public Matrix Transform(Matrix X)
    {
        if (!_isFitted)
            throw new InvalidOperationException("PCA has not been fitted.");

        int n = X.rowLength;
        int d = X.columnLength;
        var result = new double[n, NComponents];

        for (int i = 0; i < n; i++)
            for (int k = 0; k < NComponents; k++)
            {
                double dot = 0;
                for (int j = 0; j < d; j++)
                    dot += (X.values[i, j] - Mean[j]) * Components[k, j];
                result[i, k] = dot;
            }

        return new Matrix(result);
    }

    // ── Clone ────────────────────────────────────────────────────
    public IDimensionalityReducer Clone()
    {
        var clone = new PCA
        {
            NComponents = NComponents,
            MaxIterations = MaxIterations,
            Tolerance = Tolerance,
            Seed = Seed
        };

        if (_isFitted)
        {
            clone.Mean = (double[])Mean.Clone();
            clone.ExplainedVariance = (double[])ExplainedVariance.Clone();
            clone.ExplainedVarianceRatio = (double[])ExplainedVarianceRatio.Clone();
            clone.Components = (double[,])Components.Clone();
            clone._isFitted = true;
        }

        return clone;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Private helpers
    // ═══════════════════════════════════════════════════════════════

    private static double[] ComputeMean(double[,] data, int n, int d)
    {
        var mean = new double[d];
        for (int j = 0; j < d; j++)
        {
            double sum = 0;
            for (int i = 0; i < n; i++)
                sum += data[i, j];
            mean[j] = sum / n;
        }
        return mean;
    }

    private static double[,] CenterData(double[,] data, double[] mean, int n, int d)
    {
        var centered = new double[n, d];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < d; j++)
                centered[i, j] = data[i, j] - mean[j];
        return centered;
    }

    private static double[,] ComputeCovariance(double[,] centered, int n, int d)
    {
        var cov = new double[d, d];
        for (int i = 0; i < d; i++)
            for (int j = i; j < d; j++)
            {
                double sum = 0;
                for (int k = 0; k < n; k++)
                    sum += centered[k, i] * centered[k, j];
                cov[i, j] = sum / (n - 1);
                cov[j, i] = cov[i, j];
            }
        return cov;
    }

    private (double eigenvalue, double[] eigenvector) PowerIteration(
        double[,] matrix, int d, Random rng)
    {
        // Random initial vector
        var v = new double[d];
        for (int i = 0; i < d; i++)
            v[i] = rng.NextDouble() - 0.5;
        Normalize(v);

        double eigenvalue = 0;

        for (int iter = 0; iter < MaxIterations; iter++)
        {
            // w = A * v
            var w = new double[d];
            for (int i = 0; i < d; i++)
            {
                double sum = 0;
                for (int j = 0; j < d; j++)
                    sum += matrix[i, j] * v[j];
                w[i] = sum;
            }

            // Rayleigh quotient: eigenvalue = v^T * w
            double newEigenvalue = 0;
            for (int i = 0; i < d; i++)
                newEigenvalue += v[i] * w[i];

            // Normalize
            Normalize(w);

            if (Math.Abs(newEigenvalue - eigenvalue) < Tolerance)
            {
                eigenvalue = newEigenvalue;
                v = w;
                break;
            }

            eigenvalue = newEigenvalue;
            v = w;
        }

        return (eigenvalue, v);
    }

    private static void Normalize(double[] v)
    {
        double norm = 0;
        for (int i = 0; i < v.Length; i++)
            norm += v[i] * v[i];
        norm = Math.Sqrt(norm);

        if (norm > 1e-15)
            for (int i = 0; i < v.Length; i++)
                v[i] /= norm;
    }
}
