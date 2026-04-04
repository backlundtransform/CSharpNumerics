using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Statistics.Fitting;

/// <summary>
/// Internal helper for linear algebra operations used across fitting algorithms.
/// </summary>
internal static class FittingSolver
{
    /// <summary>Solves the linear system Ax = b using Gaussian elimination with partial pivoting.</summary>
    internal static double[] Solve(double[,] A, double[] b)
    {
        int n = b.Length;
        double[,] aug = new double[n, n + 1];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                aug[i, j] = A[i, j];
            aug[i, n] = b[i];
        }

        for (int col = 0; col < n; col++)
        {
            int maxRow = col;
            double maxVal = Math.Abs(aug[col, col]);
            for (int row = col + 1; row < n; row++)
            {
                double absVal = Math.Abs(aug[row, col]);
                if (absVal > maxVal)
                {
                    maxVal = absVal;
                    maxRow = row;
                }
            }

            if (maxRow != col)
            {
                for (int j = 0; j <= n; j++)
                {
                    double tmp = aug[col, j];
                    aug[col, j] = aug[maxRow, j];
                    aug[maxRow, j] = tmp;
                }
            }

            double pivot = aug[col, col];
            if (Math.Abs(pivot) < 1e-14)
                throw new InvalidOperationException("Singular matrix encountered during fitting.");

            for (int row = col + 1; row < n; row++)
            {
                double factor = aug[row, col] / pivot;
                for (int j = col; j <= n; j++)
                    aug[row, j] -= factor * aug[col, j];
            }
        }

        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--)
        {
            x[i] = aug[i, n];
            for (int j = i + 1; j < n; j++)
                x[i] -= aug[i, j] * x[j];
            x[i] /= aug[i, i];
        }

        return x;
    }

    /// <summary>Inverts a matrix using Gauss-Jordan elimination with partial pivoting.</summary>
    internal static double[,] Invert(double[,] A, int n)
    {
        double[,] aug = new double[n, 2 * n];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                aug[i, j] = A[i, j];
            aug[i, n + i] = 1.0;
        }

        for (int col = 0; col < n; col++)
        {
            int maxRow = col;
            double maxVal = Math.Abs(aug[col, col]);
            for (int row = col + 1; row < n; row++)
            {
                double absVal = Math.Abs(aug[row, col]);
                if (absVal > maxVal)
                {
                    maxVal = absVal;
                    maxRow = row;
                }
            }

            if (maxRow != col)
            {
                for (int j = 0; j < 2 * n; j++)
                {
                    double tmp = aug[col, j];
                    aug[col, j] = aug[maxRow, j];
                    aug[maxRow, j] = tmp;
                }
            }

            double pivot = aug[col, col];
            if (Math.Abs(pivot) < 1e-14)
                throw new InvalidOperationException("Singular matrix cannot be inverted.");

            for (int j = 0; j < 2 * n; j++)
                aug[col, j] /= pivot;

            for (int row = 0; row < n; row++)
            {
                if (row == col) continue;
                double factor = aug[row, col];
                for (int j = 0; j < 2 * n; j++)
                    aug[row, j] -= factor * aug[col, j];
            }
        }

        double[,] inv = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                inv[i, j] = aug[i, n + j];

        return inv;
    }

    /// <summary>Builds a Vandermonde design matrix [1, x, x², …, x^degree].</summary>
    internal static double[,] BuildPolynomialDesignMatrix(double[] x, int degree)
    {
        int n = x.Length;
        int cols = degree + 1;
        double[,] X = new double[n, cols];
        for (int i = 0; i < n; i++)
        {
            X[i, 0] = 1.0;
            for (int j = 1; j < cols; j++)
                X[i, j] = X[i, j - 1] * x[i];
        }
        return X;
    }

    /// <summary>Builds a Vandermonde design matrix from a VectorN.</summary>
    internal static double[,] BuildPolynomialDesignMatrix(VectorN x, int degree)
    {
        return BuildPolynomialDesignMatrix(x.Values, degree);
    }

    /// <summary>Computes A^T · A.</summary>
    internal static double[,] MultiplyATA(double[,] A, int rows, int cols)
    {
        double[,] result = new double[cols, cols];
        for (int i = 0; i < cols; i++)
            for (int j = i; j < cols; j++)
            {
                double sum = 0;
                for (int k = 0; k < rows; k++)
                    sum += A[k, i] * A[k, j];
                result[i, j] = sum;
                result[j, i] = sum;
            }
        return result;
    }

    /// <summary>Computes A^T · b.</summary>
    internal static double[] MultiplyATb(double[,] A, double[] b, int rows, int cols)
    {
        double[] result = new double[cols];
        for (int j = 0; j < cols; j++)
        {
            double sum = 0;
            for (int i = 0; i < rows; i++)
                sum += A[i, j] * b[i];
            result[j] = sum;
        }
        return result;
    }

    /// <summary>Computes A^T · diag(w) · A.</summary>
    internal static double[,] MultiplyATWA(double[,] A, double[] w, int rows, int cols)
    {
        double[,] result = new double[cols, cols];
        for (int i = 0; i < cols; i++)
            for (int j = i; j < cols; j++)
            {
                double sum = 0;
                for (int k = 0; k < rows; k++)
                    sum += A[k, i] * w[k] * A[k, j];
                result[i, j] = sum;
                result[j, i] = sum;
            }
        return result;
    }

    /// <summary>Computes A^T · diag(w) · b.</summary>
    internal static double[] MultiplyATWb(double[,] A, double[] w, double[] b, int rows, int cols)
    {
        double[] result = new double[cols];
        for (int j = 0; j < cols; j++)
        {
            double sum = 0;
            for (int i = 0; i < rows; i++)
                sum += A[i, j] * w[i] * b[i];
            result[j] = sum;
        }
        return result;
    }

    /// <summary>Computes fitted values ŷ = X · β.</summary>
    internal static double[] ComputeFitted(double[,] X, double[] beta, int rows, int cols)
    {
        double[] fitted = new double[rows];
        for (int i = 0; i < rows; i++)
        {
            double sum = 0;
            for (int j = 0; j < cols; j++)
                sum += X[i, j] * beta[j];
            fitted[i] = sum;
        }
        return fitted;
    }

    /// <summary>Computes SE[i] = sqrt(residualVariance · diag((X^TX)^−1)[i]).</summary>
    internal static double[] ComputeStandardErrors(double[,] XtXInverse, double residualVariance, int p)
    {
        double[] se = new double[p];
        for (int i = 0; i < p; i++)
        {
            double val = residualVariance * XtXInverse[i, i];
            se[i] = val > 0 ? Math.Sqrt(val) : 0.0;
        }
        return se;
    }

    /// <summary>Computes the median of an array (non-destructive).</summary>
    internal static double Median(double[] values)
    {
        double[] sorted = (double[])values.Clone();
        Array.Sort(sorted);
        int n = sorted.Length;
        if (n % 2 == 0)
            return (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0;
        return sorted[n / 2];
    }

    /// <summary>Computes the median absolute deviation (MAD).</summary>
    internal static double MAD(double[] values)
    {
        double med = Median(values);
        double[] absDevs = new double[values.Length];
        for (int i = 0; i < values.Length; i++)
            absDevs[i] = Math.Abs(values[i] - med);
        return Median(absDevs);
    }

    /// <summary>Inverse normal CDF (Acklam's rational approximation, relative error &lt; 1.15e-9).</summary>
    internal static double InverseNormalCDF(double p)
    {
        const double a1 = -3.969683028665376e+01;
        const double a2 = 2.209460984245205e+02;
        const double a3 = -2.759285104469687e+02;
        const double a4 = 1.383577518672690e+02;
        const double a5 = -3.066479806614716e+01;
        const double a6 = 2.506628277459239e+00;

        const double b1 = -5.447609879822406e+01;
        const double b2 = 1.615858368580409e+02;
        const double b3 = -1.556989798598866e+02;
        const double b4 = 6.680131188771972e+01;
        const double b5 = -1.328068155288572e+01;

        const double c1 = -7.784894002430293e-03;
        const double c2 = -3.223964580411365e-01;
        const double c3 = -2.400758277161838e+00;
        const double c4 = -2.549732539343734e+00;
        const double c5 = 4.374664141464968e+00;
        const double c6 = 2.938163982698783e+00;

        const double d1 = 7.784695709041462e-03;
        const double d2 = 3.224671290700398e-01;
        const double d3 = 2.445134137142996e+00;
        const double d4 = 3.754408661907416e+00;

        const double pLow = 0.02425;
        const double pHigh = 1.0 - pLow;

        if (p < pLow)
        {
            double q = Math.Sqrt(-2.0 * Math.Log(p));
            return (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
                   ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
        }

        if (p <= pHigh)
        {
            double q = p - 0.5;
            double r = q * q;
            return (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
                   (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0);
        }

        {
            double q = Math.Sqrt(-2.0 * Math.Log(1.0 - p));
            return -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
                    ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
        }
    }

    /// <summary>
    /// Approximation of the Student's t-distribution quantile using Cornish-Fisher expansion.
    /// </summary>
    internal static double TQuantile(double p, int df)
    {
        if (p < 0.5)
            return -TQuantile(1.0 - p, df);

        if (df == 1)
            return Math.Tan(Math.PI * (p - 0.5));

        double z = InverseNormalCDF(p);

        if (df >= 1000)
            return z;

        double v = df;
        double z2 = z * z;
        double z3 = z * z2;
        double z5 = z3 * z2;
        double z7 = z5 * z2;
        double z9 = z7 * z2;

        double g1 = (z3 + z) / 4.0;
        double g2 = (5 * z5 + 16 * z3 + 3 * z) / 96.0;
        double g3 = (3 * z7 + 19 * z5 + 17 * z3 - 15 * z) / 384.0;
        double g4 = (79 * z9 + 776 * z7 + 1482 * z5 - 1920 * z3 - 945 * z) / 92160.0;

        return z + g1 / v + g2 / (v * v) + g3 / (v * v * v) + g4 / (v * v * v * v);
    }
}
