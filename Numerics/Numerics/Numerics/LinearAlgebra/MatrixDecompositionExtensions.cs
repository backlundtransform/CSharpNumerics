using CSharpNumerics.Numerics.LinearAlgebra.Decompositions;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Numerics.LinearAlgebra;

public static class MatrixDecompositionExtensions
{
    /// <summary>
    /// Computes the LU decomposition with partial pivoting: P·A = L·U.
    /// The returned object caches the factorization for reuse across multiple solves.
    /// </summary>
    public static LuDecomposition Lu(this Matrix matrix) => new LuDecomposition(matrix);

    /// <summary>
    /// Computes the Cholesky decomposition A = L·Lᵀ for a symmetric positive definite matrix.
    /// The returned object caches the factorization for reuse across multiple solves.
    /// </summary>
    public static CholeskyDecomposition Cholesky(this Matrix matrix) => new CholeskyDecomposition(matrix);
}
