using System;
using System.Collections.Generic;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Numerics.LinearAlgebra.Decompositions;

/// <summary>
/// QR decomposition via Householder reflections: A = Q·R where Q has orthonormal columns
/// and R is upper triangular. For an m×n matrix with m ≥ n, Solve computes the least squares
/// solution min ‖A x − b‖, which makes QR the numerically stable choice for overdetermined
/// systems (curve fitting, linear regression) — no normal equations needed.
/// The factorization is computed once and can be reused for multiple solves.
/// </summary>
public sealed class QrDecomposition
{
    private readonly double[,] qr;
    private readonly double[] rDiag;
    private readonly int m;
    private readonly int n;

    public QrDecomposition(Matrix matrix)
    {
        if (matrix.rowLength < matrix.columnLength)
        {
            throw new Exception("The matrix must have at least as many rows as columns");
        }

        m = matrix.rowLength;
        n = matrix.columnLength;
        qr = (double[,])matrix.values.Clone();
        rDiag = new double[n];

        for (var k = 0; k < n; k++)
        {
            var norm = 0.0;
            for (var i = k; i < m; i++)
            {
                norm = Hypot(norm, qr[i, k]);
            }

            if (norm != 0.0)
            {
                if (qr[k, k] < 0)
                {
                    norm = -norm;
                }

                for (var i = k; i < m; i++)
                {
                    qr[i, k] /= norm;
                }
                qr[k, k] += 1.0;

                for (var j = k + 1; j < n; j++)
                {
                    var s = 0.0;
                    for (var i = k; i < m; i++)
                    {
                        s += qr[i, k] * qr[i, j];
                    }
                    s = -s / qr[k, k];
                    for (var i = k; i < m; i++)
                    {
                        qr[i, j] += s * qr[i, k];
                    }
                }
            }

            rDiag[k] = -norm;
        }
    }

    /// <summary>
    /// True if the columns of A are linearly independent. Diagonal entries of R are compared
    /// against a tolerance relative to the largest one, since Householder rounding leaves
    /// rank-deficient columns at ~1e-15 rather than exactly zero.
    /// </summary>
    public bool IsFullRank
    {
        get
        {
            var maxDiag = 0.0;
            for (var j = 0; j < n; j++)
            {
                maxDiag = Math.Max(maxDiag, Math.Abs(rDiag[j]));
            }

            if (maxDiag == 0.0)
            {
                return false;
            }

            var threshold = maxDiag * Math.Max(m, n) * Math.Pow(2.0, -52.0);
            for (var j = 0; j < n; j++)
            {
                if (Math.Abs(rDiag[j]) <= threshold)
                {
                    return false;
                }
            }
            return true;
        }
    }

    /// <summary>
    /// The economy-size orthogonal factor Q (m×n) with orthonormal columns.
    /// </summary>
    public Matrix Q
    {
        get
        {
            var q = new double[m, n];
            for (var k = n - 1; k >= 0; k--)
            {
                for (var i = 0; i < m; i++)
                {
                    q[i, k] = 0.0;
                }
                q[k, k] = 1.0;

                for (var j = k; j < n; j++)
                {
                    if (qr[k, k] != 0.0)
                    {
                        var s = 0.0;
                        for (var i = k; i < m; i++)
                        {
                            s += qr[i, k] * q[i, j];
                        }
                        s = -s / qr[k, k];
                        for (var i = k; i < m; i++)
                        {
                            q[i, j] += s * qr[i, k];
                        }
                    }
                }
            }
            return new Matrix(q);
        }
    }

    /// <summary>
    /// The upper triangular factor R (n×n).
    /// </summary>
    public Matrix R
    {
        get
        {
            var r = new double[n, n];
            for (var i = 0; i < n; i++)
            {
                for (var j = i; j < n; j++)
                {
                    r[i, j] = i == j ? rDiag[i] : qr[i, j];
                }
            }
            return new Matrix(r);
        }
    }

    /// <summary>
    /// Solves A x = b in the least squares sense: returns x minimizing ‖A x − b‖.
    /// For square full-rank systems this is the exact solution.
    /// </summary>
    public VectorN Solve(VectorN b)
    {
        if (b.Length != m)
        {
            throw new Exception("The vector length must match the matrix row length");
        }

        var x = new double[m];
        for (var i = 0; i < m; i++)
        {
            x[i] = b[i];
        }

        SolveInPlace(x);

        var result = new double[n];
        Array.Copy(x, result, n);
        return new VectorN(result);
    }

    /// <summary>
    /// Solves A x = b in the least squares sense.
    /// </summary>
    public List<double> Solve(List<double> b)
    {
        if (b.Count != m)
        {
            throw new Exception("The vector length must match the matrix row length");
        }

        var x = b.ToArray();
        SolveInPlace(x);

        var result = new List<double>(n);
        for (var i = 0; i < n; i++)
        {
            result.Add(x[i]);
        }
        return result;
    }

    /// <summary>
    /// Solves A X = B in the least squares sense for multiple right-hand sides (one per column of B).
    /// </summary>
    public Matrix Solve(Matrix b)
    {
        if (b.rowLength != m)
        {
            throw new Exception("The row length of B must match the matrix row length");
        }

        var columns = b.columnLength;
        var x = new double[n, columns];

        for (var c = 0; c < columns; c++)
        {
            var column = new double[m];
            for (var i = 0; i < m; i++)
            {
                column[i] = b.values[i, c];
            }

            SolveInPlace(column);

            for (var i = 0; i < n; i++)
            {
                x[i, c] = column[i];
            }
        }

        return new Matrix(x);
    }

    private void SolveInPlace(double[] x)
    {
        if (!IsFullRank)
        {
            throw new Exception("The matrix is rank deficient");
        }

        for (var k = 0; k < n; k++)
        {
            var s = 0.0;
            for (var i = k; i < m; i++)
            {
                s += qr[i, k] * x[i];
            }
            s = -s / qr[k, k];
            for (var i = k; i < m; i++)
            {
                x[i] += s * qr[i, k];
            }
        }

        for (var k = n - 1; k >= 0; k--)
        {
            x[k] /= rDiag[k];
            for (var i = 0; i < k; i++)
            {
                x[i] -= x[k] * qr[i, k];
            }
        }
    }

    private static double Hypot(double a, double b)
    {
        if (Math.Abs(a) > Math.Abs(b))
        {
            var r = b / a;
            return Math.Abs(a) * Math.Sqrt(1 + r * r);
        }
        if (b != 0)
        {
            var r = a / b;
            return Math.Abs(b) * Math.Sqrt(1 + r * r);
        }
        return 0.0;
    }
}
