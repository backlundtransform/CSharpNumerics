using System;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Cyclic Jacobi eigenvalue solver for real symmetric matrices. Unlike power iteration with
/// deflation, it returns <em>all</em> eigenvalues and orthonormal eigenvectors accurately,
/// including the smallest ones — which is exactly what the Schrödinger Hamiltonian needs (the
/// low-lying energy levels are the smallest eigenvalues). Eigenvalues are returned in ascending
/// order with matching eigenvectors.
/// </summary>
internal static class SymmetricEigenSolver
{
    /// <summary>
    /// Diagonalises a symmetric matrix. Returns ascending eigenvalues and the corresponding
    /// eigenvectors as <c>vectors[k]</c> (length n).
    /// </summary>
    public static (double[] values, double[][] vectors) Solve(double[,] matrix, int maxSweeps = 100, double tolerance = 1e-12)
    {
        int n = matrix.GetLength(0);
        if (matrix.GetLength(1) != n)
            throw new ArgumentException("Matrix must be square.", nameof(matrix));

        var a = (double[,])matrix.Clone();
        var v = new double[n, n];
        for (int i = 0; i < n; i++) v[i, i] = 1.0;

        for (int sweep = 0; sweep < maxSweeps; sweep++)
        {
            double off = OffDiagonalNorm(a, n);
            if (off < tolerance)
                break;

            for (int p = 0; p < n - 1; p++)
            {
                for (int q = p + 1; q < n; q++)
                {
                    if (Math.Abs(a[p, q]) < 1e-300)
                        continue;

                    // Rotation angle that zeroes the (p,q) entry of JᵀAJ:
                    // tan(2φ) = 2·a_pq / (a_qq − a_pp).
                    double phi = 0.5 * Math.Atan2(2.0 * a[p, q], a[q, q] - a[p, p]);
                    double c = Math.Cos(phi);
                    double s = Math.Sin(phi);

                    // A := A·J  (rotate columns p, q)
                    for (int i = 0; i < n; i++)
                    {
                        double aip = a[i, p], aiq = a[i, q];
                        a[i, p] = (c * aip) - (s * aiq);
                        a[i, q] = (s * aip) + (c * aiq);
                    }
                    // A := Jᵀ·A  (rotate rows p, q)
                    for (int j = 0; j < n; j++)
                    {
                        double apj = a[p, j], aqj = a[q, j];
                        a[p, j] = (c * apj) - (s * aqj);
                        a[q, j] = (s * apj) + (c * aqj);
                    }
                    // Accumulate eigenvectors V := V·J
                    for (int i = 0; i < n; i++)
                    {
                        double vip = v[i, p], viq = v[i, q];
                        v[i, p] = (c * vip) - (s * viq);
                        v[i, q] = (s * vip) + (c * viq);
                    }
                }
            }
        }

        // Extract eigenvalues (diagonal) and sort ascending.
        var values = new double[n];
        for (int i = 0; i < n; i++) values[i] = a[i, i];

        var order = new int[n];
        for (int i = 0; i < n; i++) order[i] = i;
        Array.Sort(order, (x, y) => values[x].CompareTo(values[y]));

        var sortedValues = new double[n];
        var vectors = new double[n][];
        for (int k = 0; k < n; k++)
        {
            int col = order[k];
            sortedValues[k] = values[col];
            var vec = new double[n];
            for (int i = 0; i < n; i++) vec[i] = v[i, col];
            vectors[k] = vec;
        }

        return (sortedValues, vectors);
    }

    private static double OffDiagonalNorm(double[,] a, int n)
    {
        double sum = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                sum += a[i, j] * a[i, j];
        return Math.Sqrt(2.0 * sum);
    }
}
