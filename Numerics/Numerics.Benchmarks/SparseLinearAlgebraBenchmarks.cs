using BenchmarkDotNet.Attributes;
using CSharpNumerics.Numerics.Objects;

namespace Numerics.Benchmarks;

/// <summary>
/// Baseline for the sparse path: SpMV and the preconditioned conjugate gradient solver
/// that <c>Assembler2D</c> runs on.
/// </summary>
/// <remarks>
/// The system is the 5-point Laplacian on a <see cref="Grid"/> × <see cref="Grid"/> mesh —
/// the same shape a 2-D finite element or finite difference problem produces, so the
/// numbers transfer to real FEM workloads.
/// </remarks>
[MemoryDiagnoser]
public class SparseLinearAlgebraBenchmarks
{
    [Params(32, 128)]
    public int Grid;

    private SparseMatrix laplacian = default!;
    private VectorN rightHandSide = default!;

    /// <summary>Degrees of freedom in the assembled system — the interior nodes.</summary>
    public int Unknowns => (Grid - 2) * (Grid - 2);

    [GlobalSetup]
    public void Setup()
    {
        laplacian = BuildLaplacian(Grid);

        var values = new double[Unknowns];
        for (var i = 0; i < values.Length; i++)
        {
            values[i] = 1.0;
        }

        rightHandSide = new VectorN(values);
    }

    [Benchmark(Description = "SpMV (sparse matrix × vector)")]
    public VectorN SparseMatrixVectorMultiply() => laplacian.Multiply(rightHandSide);

    /// <summary>
    /// Conjugate gradient with a Jacobi preconditioner. Iteration count grows with the
    /// grid, so this scales worse than SpMV alone — that growth is what stronger
    /// preconditioners would attack.
    /// </summary>
    [Benchmark(Description = "PCG solve")]
    public VectorN SolvePcg() => laplacian.SolvePCG(rightHandSide, tolerance: 1e-8);

    /// <summary>
    /// Assembles the 5-point Laplacian over the interior nodes only, with the Dirichlet
    /// boundary folded into the right-hand side rather than kept as rows.
    /// </summary>
    /// <remarks>
    /// Eliminating the boundary rather than pinning it with identity rows is what keeps the
    /// matrix symmetric: an identity row zeroes a row but leaves the matching column entry
    /// in its interior neighbours, and conjugate gradient has no convergence guarantee on a
    /// non-symmetric system.
    /// </remarks>
    private static SparseMatrix BuildLaplacian(int grid)
    {
        var interior = grid - 2;
        var triplets = new List<(int row, int col, double val)>();
        int Index(int x, int y) => y * interior + x;

        for (var y = 0; y < interior; y++)
        {
            for (var x = 0; x < interior; x++)
            {
                var row = Index(x, y);
                triplets.Add((row, row, 4.0));

                if (x > 0)
                {
                    triplets.Add((row, Index(x - 1, y), -1.0));
                }

                if (x < interior - 1)
                {
                    triplets.Add((row, Index(x + 1, y), -1.0));
                }

                if (y > 0)
                {
                    triplets.Add((row, Index(x, y - 1), -1.0));
                }

                if (y < interior - 1)
                {
                    triplets.Add((row, Index(x, y + 1), -1.0));
                }
            }
        }

        return SparseMatrix.FromTriplets(interior * interior, interior * interior, triplets);
    }
}
