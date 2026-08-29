using System;
using System.Collections.Generic;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Numerics.LinearAlgebra.Decompositions;

/// <summary>
/// LU decomposition with partial (row) pivoting: P·A = L·U where L is unit lower triangular,
/// U is upper triangular and P is a permutation matrix.
/// The factorization is computed once and can be reused for multiple solves.
/// </summary>
public sealed class LuDecomposition
{
    private readonly double[,] lu;
    private readonly int[] pivot;
    private readonly int pivotSign;
    private readonly int n;

    public LuDecomposition(Matrix matrix)
    {
        if (matrix.rowLength != matrix.columnLength)
        {
            throw new Exception("Is not a NxN matrix");
        }

        n = matrix.rowLength;
        lu = (double[,])matrix.values.Clone();
        pivot = new int[n];
        pivotSign = 1;

        for (var i = 0; i < n; i++)
        {
            pivot[i] = i;
        }

        for (var j = 0; j < n; j++)
        {
            for (var i = 0; i < n; i++)
            {
                var kMax = Math.Min(i, j);
                var sum = 0.0;
                for (var k = 0; k < kMax; k++)
                {
                    sum += lu[i, k] * lu[k, j];
                }
                lu[i, j] -= sum;
            }

            var p = j;
            for (var i = j + 1; i < n; i++)
            {
                if (Math.Abs(lu[i, j]) > Math.Abs(lu[p, j]))
                {
                    p = i;
                }
            }

            if (p != j)
            {
                for (var k = 0; k < n; k++)
                {
                    (lu[p, k], lu[j, k]) = (lu[j, k], lu[p, k]);
                }

                (pivot[p], pivot[j]) = (pivot[j], pivot[p]);
                pivotSign = -pivotSign;
            }

            if (lu[j, j] != 0.0)
            {
                for (var i = j + 1; i < n; i++)
                {
                    lu[i, j] /= lu[j, j];
                }
            }
        }
    }

    /// <summary>
    /// True if the matrix is singular (a zero pivot was encountered) and cannot be solved.
    /// </summary>
    public bool IsSingular
    {
        get
        {
            for (var j = 0; j < n; j++)
            {
                if (lu[j, j] == 0.0)
                {
                    return true;
                }
            }
            return false;
        }
    }

    /// <summary>
    /// The unit lower triangular factor L.
    /// </summary>
    public Matrix Lower
    {
        get
        {
            var l = new double[n, n];
            for (var i = 0; i < n; i++)
            {
                for (var j = 0; j < n; j++)
                {
                    l[i, j] = i > j ? lu[i, j] : (i == j ? 1.0 : 0.0);
                }
            }
            return new Matrix(l);
        }
    }

    /// <summary>
    /// The upper triangular factor U.
    /// </summary>
    public Matrix Upper
    {
        get
        {
            var u = new double[n, n];
            for (var i = 0; i < n; i++)
            {
                for (var j = i; j < n; j++)
                {
                    u[i, j] = lu[i, j];
                }
            }
            return new Matrix(u);
        }
    }

    /// <summary>
    /// The row permutation applied by pivoting: row i of P·A is row Permutation[i] of A.
    /// </summary>
    public int[] Permutation => (int[])pivot.Clone();

    /// <summary>
    /// The permutation matrix P such that P·A = L·U.
    /// </summary>
    public Matrix PermutationMatrix
    {
        get
        {
            var p = new double[n, n];
            for (var i = 0; i < n; i++)
            {
                p[i, pivot[i]] = 1.0;
            }
            return new Matrix(p);
        }
    }

    /// <summary>
    /// Computes the determinant from the factorization as sign(P) · ∏ U[i,i].
    /// </summary>
    public double Determinant()
    {
        var determinant = (double)pivotSign;
        for (var j = 0; j < n; j++)
        {
            determinant *= lu[j, j];
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
            x[i] = b[pivot[i]];
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

        var x = new double[n];
        for (var i = 0; i < n; i++)
        {
            x[i] = b[pivot[i]];
        }

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
                column[i] = b.values[pivot[i], c];
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
        if (IsSingular)
        {
            throw new Exception("This matrix is not invertible");
        }

        for (var i = 1; i < n; i++)
        {
            for (var k = 0; k < i; k++)
            {
                x[i] -= lu[i, k] * x[k];
            }
        }

        for (var i = n - 1; i >= 0; i--)
        {
            for (var k = i + 1; k < n; k++)
            {
                x[i] -= lu[i, k] * x[k];
            }
            x[i] /= lu[i, i];
        }
    }
}
