using System;
using System.Collections.Generic;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Numerics.LinearAlgebra.Decompositions;

/// <summary>
/// Cholesky decomposition A = L·Lᵀ for symmetric positive definite matrices,
/// where L is lower triangular. Roughly twice as fast as LU for such systems.
/// The factorization is computed once and can be reused for multiple solves.
/// </summary>
public sealed class CholeskyDecomposition
{
    private const double SymmetryTolerance = 1e-10;

    private readonly double[,] l;
    private readonly int n;

    public CholeskyDecomposition(Matrix matrix)
    {
        if (matrix.rowLength != matrix.columnLength)
        {
            throw new Exception("Is not a NxN matrix");
        }

        n = matrix.rowLength;
        l = new double[n, n];
        IsPositiveDefinite = IsSymmetric(matrix);

        for (var j = 0; j < n; j++)
        {
            var d = matrix.values[j, j];
            for (var k = 0; k < j; k++)
            {
                d -= l[j, k] * l[j, k];
            }

            if (d <= 0.0)
            {
                IsPositiveDefinite = false;
                return;
            }

            l[j, j] = Math.Sqrt(d);

            for (var i = j + 1; i < n; i++)
            {
                var sum = matrix.values[i, j];
                for (var k = 0; k < j; k++)
                {
                    sum -= l[i, k] * l[j, k];
                }
                l[i, j] = sum / l[j, j];
            }
        }
    }

    /// <summary>
    /// True if the matrix is symmetric and all pivots are positive, i.e. the factorization A = L·Lᵀ exists.
    /// </summary>
    public bool IsPositiveDefinite { get; }

    /// <summary>
    /// The lower triangular factor L.
    /// </summary>
    public Matrix Lower
    {
        get
        {
            EnsurePositiveDefinite();
            return new Matrix((double[,])l.Clone());
        }
    }

    /// <summary>
    /// Computes the determinant from the factorization as ∏ L[i,i]².
    /// </summary>
    public double Determinant()
    {
        EnsurePositiveDefinite();
        var determinant = 1.0;
        for (var j = 0; j < n; j++)
        {
            determinant *= l[j, j] * l[j, j];
        }
        return determinant;
    }

    /// <summary>
    /// Solves A x = b using the cached factorization (forward and back substitution).
    /// </summary>
    public VectorN Solve(VectorN b)
    {
        if (b.Length != n)
        {
            throw new Exception("The vector length must match the matrix dimension");
        }

        var x = new double[n];
        for (var i = 0; i < n; i++)
        {
            x[i] = b[i];
        }

        SolveInPlace(x);
        return new VectorN(x);
    }

    /// <summary>
    /// Solves A x = b using the cached factorization.
    /// </summary>
    public List<double> Solve(List<double> b)
    {
        if (b.Count != n)
        {
            throw new Exception("The vector length must match the matrix dimension");
        }

        var x = b.ToArray();
        SolveInPlace(x);
        return new List<double>(x);
    }

    /// <summary>
    /// Solves A X = B for multiple right-hand sides (one per column of B).
    /// </summary>
    public Matrix Solve(Matrix b)
    {
        if (b.rowLength != n)
        {
            throw new Exception("The row length of B must match the matrix dimension");
        }

        var columns = b.columnLength;
        var x = new double[n, columns];

        for (var c = 0; c < columns; c++)
        {
            var column = new double[n];
            for (var i = 0; i < n; i++)
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

    /// <summary>
    /// Computes the inverse by solving A X = I. Reuses the cached factorization.
    /// </summary>
    public Matrix Inverse()
    {
        var identity = new double[n, n];
        for (var i = 0; i < n; i++)
        {
            identity[i, i] = 1.0;
        }
        return Solve(new Matrix(identity));
    }

    private void SolveInPlace(double[] x)
    {
        EnsurePositiveDefinite();

        for (var i = 0; i < n; i++)
        {
            for (var k = 0; k < i; k++)
            {
                x[i] -= l[i, k] * x[k];
            }
            x[i] /= l[i, i];
        }

        for (var i = n - 1; i >= 0; i--)
        {
            for (var k = i + 1; k < n; k++)
            {
                x[i] -= l[k, i] * x[k];
            }
            x[i] /= l[i, i];
        }
    }

    private void EnsurePositiveDefinite()
    {
        if (!IsPositiveDefinite)
        {
            throw new Exception("The matrix is not symmetric positive definite");
        }
    }

    private static bool IsSymmetric(Matrix matrix)
    {
        for (var i = 0; i < matrix.rowLength; i++)
        {
            for (var j = i + 1; j < matrix.columnLength; j++)
            {
                var a = matrix.values[i, j];
                var b = matrix.values[j, i];
                var scale = Math.Max(1.0, Math.Max(Math.Abs(a), Math.Abs(b)));

                if (Math.Abs(a - b) > SymmetryTolerance * scale)
                {
                    return false;
                }
            }
        }
        return true;
    }
}
