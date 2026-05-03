namespace CSharpNumerics.Numerics.Objects;

using System;
using System.Collections.Generic;

/// <summary>
/// Sparse matrix stored in Compressed Sparse Row (CSR) format.
/// Efficient for matrix–vector products and iterative solvers on large, sparse systems.
/// </summary>
public sealed class SparseMatrix
{
    /// <summary>Non-zero values, row by row.</summary>
    public double[] Values { get; }

    /// <summary>Column index for each entry in <see cref="Values"/>.</summary>
    public int[] ColIndices { get; }

    /// <summary>Row pointers: RowPointers[i] is the index into Values where row i starts.
    /// RowPointers[Rows] == NonZeroCount.</summary>
    public int[] RowPointers { get; }

    /// <summary>Number of rows.</summary>
    public int Rows { get; }

    /// <summary>Number of columns.</summary>
    public int Cols { get; }

    /// <summary>Number of stored non-zero entries.</summary>
    public int NonZeroCount => Values.Length;

    private SparseMatrix(int rows, int cols, double[] values, int[] colIndices, int[] rowPointers)
    {
        Rows = rows;
        Cols = cols;
        Values = values;
        ColIndices = colIndices;
        RowPointers = rowPointers;
    }

    /// <summary>
    /// Builds a CSR sparse matrix from a list of (row, col, value) triplets.
    /// Duplicate (row, col) entries are summed — standard FEM assembly behaviour.
    /// </summary>
    /// <param name="rows">Number of rows.</param>
    /// <param name="cols">Number of columns.</param>
    /// <param name="triplets">List of (row, col, value) entries.</param>
    public static SparseMatrix FromTriplets(int rows, int cols, List<(int row, int col, double val)> triplets)
    {
        if (rows <= 0 || cols <= 0)
            throw new ArgumentException("Matrix dimensions must be positive.");

        // Sort by row, then by column
        triplets.Sort((a, b) =>
        {
            int cmp = a.row.CompareTo(b.row);
            return cmp != 0 ? cmp : a.col.CompareTo(b.col);
        });

        // Merge duplicates and build CSR arrays
        var vals = new List<double>();
        var colIdx = new List<int>();
        var rowPtr = new int[rows + 1];

        int currentRow = 0;
        rowPtr[0] = 0;

        for (int i = 0; i < triplets.Count; i++)
        {
            var (r, c, v) = triplets[i];
            if (r < 0 || r >= rows || c < 0 || c >= cols)
                throw new ArgumentOutOfRangeException($"Triplet ({r},{c}) is out of bounds for {rows}×{cols} matrix.");

            // Fill empty rows up to this triplet's row
            while (currentRow < r)
            {
                currentRow++;
                rowPtr[currentRow] = vals.Count;
            }

            // Accumulate duplicates at the same (row, col)
            if (vals.Count > 0 && colIdx[colIdx.Count - 1] == c && currentRow == r
                && vals.Count > rowPtr[r])
            {
                vals[vals.Count - 1] += v;
            }
            else
            {
                vals.Add(v);
                colIdx.Add(c);
            }
        }

        // Fill remaining row pointers
        currentRow++;
        while (currentRow <= rows)
        {
            rowPtr[currentRow] = vals.Count;
            currentRow++;
        }

        return new SparseMatrix(rows, cols, vals.ToArray(), colIdx.ToArray(), rowPtr);
    }

    /// <summary>
    /// Sparse matrix–vector product: y = A · x.
    /// </summary>
    public VectorN Multiply(VectorN x)
    {
        if (x.Length != Cols)
            throw new ArgumentException($"Vector length {x.Length} does not match matrix columns {Cols}.");

        var result = new double[Rows];
        for (int i = 0; i < Rows; i++)
        {
            double sum = 0.0;
            int start = RowPointers[i];
            int end = RowPointers[i + 1];
            for (int k = start; k < end; k++)
                sum += Values[k] * x[ColIndices[k]];
            result[i] = sum;
        }

        return new VectorN(result);
    }

    /// <summary>
    /// Extracts the main diagonal as a vector. Zero for rows with no diagonal entry.
    /// Used for Jacobi preconditioning.
    /// </summary>
    public VectorN Diagonal()
    {
        int n = Math.Min(Rows, Cols);
        var diag = new double[Rows];

        for (int i = 0; i < Rows; i++)
        {
            int start = RowPointers[i];
            int end = RowPointers[i + 1];
            for (int k = start; k < end; k++)
            {
                if (ColIndices[k] == i)
                {
                    diag[i] = Values[k];
                    break;
                }
            }
        }

        return new VectorN(diag);
    }

    /// <summary>
    /// Applies Dirichlet boundary conditions by row/column elimination.
    /// For each fixed DOF, the corresponding row and column are zeroed and the diagonal set to 1.
    /// The right-hand side vector is modified accordingly.
    /// Returns a new <see cref="SparseMatrix"/> — the original is not mutated.
    /// </summary>
    /// <param name="fixedDofs">Dictionary mapping DOF index → prescribed value.</param>
    /// <param name="rhs">Right-hand side vector (modified in place).</param>
    public SparseMatrix ApplyDirichlet(Dictionary<int, double> fixedDofs, VectorN rhs)
    {
        if (Rows != Cols)
            throw new InvalidOperationException("Dirichlet BCs require a square matrix.");

        // Build a set of fixed DOFs for fast lookup
        var fixedSet = new HashSet<int>(fixedDofs.Keys);

        // First pass: subtract column contributions from RHS
        for (int i = 0; i < Rows; i++)
        {
            if (fixedSet.Contains(i)) continue;

            int start = RowPointers[i];
            int end = RowPointers[i + 1];
            for (int k = start; k < end; k++)
            {
                int j = ColIndices[k];
                if (fixedSet.Contains(j))
                    rhs[i] -= Values[k] * fixedDofs[j];
            }
        }

        // Second pass: rebuild CSR with zeroed rows/columns for fixed DOFs
        var newTriplets = new List<(int row, int col, double val)>();

        for (int i = 0; i < Rows; i++)
        {
            if (fixedSet.Contains(i))
            {
                // Fixed row: only diagonal = 1
                newTriplets.Add((i, i, 1.0));
                rhs[i] = fixedDofs[i];
                continue;
            }

            int start = RowPointers[i];
            int end = RowPointers[i + 1];
            for (int k = start; k < end; k++)
            {
                int j = ColIndices[k];
                if (!fixedSet.Contains(j))
                    newTriplets.Add((i, j, Values[k]));
            }
        }

        return FromTriplets(Rows, Cols, newTriplets);
    }

    /// <summary>
    /// Solves A·x = b using the Preconditioned Conjugate Gradient method
    /// with diagonal (Jacobi) preconditioning.
    /// Requires A to be symmetric positive-definite.
    /// </summary>
    /// <param name="b">Right-hand side vector.</param>
    /// <param name="tolerance">Convergence tolerance on the residual norm.</param>
    /// <param name="maxIterations">Maximum number of CG iterations.</param>
    /// <returns>Solution vector x.</returns>
    public VectorN SolvePCG(VectorN b, double tolerance = 1e-10, int maxIterations = 10000)
    {
        if (Rows != Cols)
            throw new InvalidOperationException("PCG requires a square matrix.");
        if (b.Length != Rows)
            throw new ArgumentException("RHS length must match matrix size.");

        int n = Rows;

        // Diagonal preconditioner M⁻¹
        var diag = Diagonal();
        var invDiag = new double[n];
        for (int i = 0; i < n; i++)
            invDiag[i] = Math.Abs(diag[i]) > 1e-30 ? 1.0 / diag[i] : 1.0;

        // Initial guess x = 0
        var x = new VectorN(n);
        var r = new VectorN((double[])b.Values.Clone());
        var z = PreconditionMultiply(invDiag, r);
        var p = new VectorN((double[])z.Values.Clone());
        double rz = r.Dot(z);

        double bNorm = b.Norm();
        if (bNorm < 1e-30)
            return x; // zero RHS → zero solution

        for (int iter = 0; iter < maxIterations; iter++)
        {
            var ap = Multiply(p);
            double pAp = p.Dot(ap);
            if (Math.Abs(pAp) < 1e-30)
                break;

            double alpha = rz / pAp;

            // x = x + alpha * p
            for (int i = 0; i < n; i++)
                x[i] += alpha * p[i];

            // r = r - alpha * Ap
            for (int i = 0; i < n; i++)
                r[i] -= alpha * ap[i];

            double rNorm = r.Norm();
            if (rNorm / bNorm < tolerance)
                break;

            z = PreconditionMultiply(invDiag, r);
            double rzNew = r.Dot(z);
            double beta = rzNew / rz;
            rz = rzNew;

            // p = z + beta * p
            for (int i = 0; i < n; i++)
                p[i] = z[i] + beta * p[i];
        }

        return x;
    }

    private static VectorN PreconditionMultiply(double[] invDiag, VectorN v)
    {
        var result = new double[v.Length];
        for (int i = 0; i < v.Length; i++)
            result[i] = invDiag[i] * v[i];
        return new VectorN(result);
    }

    /// <summary>
    /// Returns the value at (row, col). O(nnz_per_row) lookup — use sparingly.
    /// </summary>
    public double Get(int row, int col)
    {
        int start = RowPointers[row];
        int end = RowPointers[row + 1];
        for (int k = start; k < end; k++)
        {
            if (ColIndices[k] == col)
                return Values[k];
        }
        return 0.0;
    }
}
