using BenchmarkDotNet.Attributes;
using CSharpNumerics.Numerics.LinearAlgebra;
using CSharpNumerics.Numerics.LinearAlgebra.Decompositions;
using CSharpNumerics.Numerics.Objects;

namespace Numerics.Benchmarks;

/// <summary>
/// Baseline for the dense linear algebra kernels: the triple-loop multiply and the
/// decompositions built in v4.1.
/// </summary>
/// <remarks>
/// These are the hot paths that v4.3 routes more callers through and that v4.4 plans to
/// vectorise, so this is the measurement both are judged against.
/// </remarks>
[MemoryDiagnoser]
public class DenseLinearAlgebraBenchmarks
{
    [Params(64, 256)]
    public int N;

    private Matrix general = default!;
    private Matrix other = default!;
    private Matrix symmetricPositiveDefinite = default!;
    private VectorN rightHandSide = default!;
    private LuDecomposition luFactorization = default!;

    [GlobalSetup]
    public void Setup()
    {
        // Fixed seed: the same matrices every run, so results are comparable across commits.
        var random = new Random(20260930);

        general = RandomMatrix(random, N, N);
        other = RandomMatrix(random, N, N);
        symmetricPositiveDefinite = MakeSymmetricPositiveDefinite(random, N);

        var values = new double[N];
        for (var i = 0; i < N; i++)
        {
            values[i] = random.NextDouble();
        }

        rightHandSide = new VectorN(values);
        luFactorization = general.Lu();
    }

    [Benchmark(Description = "Matrix × Matrix")]
    public Matrix MatrixMultiply() => general * other;

    [Benchmark(Description = "Matrix × Vector")]
    public VectorN MatrixVectorMultiply() => general * rightHandSide;

    [Benchmark(Description = "LU factorize")]
    public LuDecomposition LuFactorize() => general.Lu();

    [Benchmark(Description = "QR factorize")]
    public QrDecomposition QrFactorize() => general.Qr();

    [Benchmark(Description = "Cholesky factorize")]
    public CholeskyDecomposition CholeskyFactorize() => symmetricPositiveDefinite.Cholesky();

    /// <summary>Factor and solve — what a caller pays who does not keep the factorization.</summary>
    [Benchmark(Description = "LU factorize + solve")]
    public VectorN LuFactorizeAndSolve() => general.Lu().Solve(rightHandSide);

    /// <summary>
    /// Solve only, against a factorization computed once. The gap between this and
    /// <see cref="LuFactorizeAndSolve"/> is what the v4.3 migration is meant to recover
    /// wherever several right-hand sides share one matrix.
    /// </summary>
    [Benchmark(Description = "LU solve (reusing factorization)")]
    public VectorN LuSolveReusingFactorization() => luFactorization.Solve(rightHandSide);

    [Benchmark(Description = "Matrix.Inverse")]
    public Matrix Inverse() => general.Inverse();

    private static Matrix RandomMatrix(Random random, int rows, int columns)
    {
        var values = new double[rows, columns];
        for (var i = 0; i < rows; i++)
        {
            for (var j = 0; j < columns; j++)
            {
                values[i, j] = random.NextDouble() * 2.0 - 1.0;
            }
        }

        return new Matrix(values);
    }

    /// <summary>
    /// Builds BᵀB + nI, which is symmetric and diagonally dominant, so Cholesky is defined.
    /// </summary>
    private static Matrix MakeSymmetricPositiveDefinite(Random random, int n)
    {
        var b = RandomMatrix(random, n, n);
        var values = new double[n, n];

        for (var i = 0; i < n; i++)
        {
            for (var j = 0; j < n; j++)
            {
                var sum = 0.0;
                for (var k = 0; k < n; k++)
                {
                    sum += b.values[k, i] * b.values[k, j];
                }

                values[i, j] = sum;
            }

            values[i, i] += n;
        }

        return new Matrix(values);
    }
}
