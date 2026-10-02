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

    /// <summary>
    /// Computes the QR decomposition via Householder reflections: A = Q·R.
    /// For overdetermined systems (rows &gt; columns) Solve gives the least squares solution.
    /// </summary>
    public static QrDecomposition Qr(this Matrix matrix) => new QrDecomposition(matrix);

    /// <summary>
    /// Computes the eigenvalue decomposition A·V = V·D.
    /// Symmetric matrices give real, ascending eigenvalues with orthonormal eigenvectors.
    /// </summary>
    public static EigenDecomposition Eigen(this Matrix matrix) => new EigenDecomposition(matrix);

    /// <summary>
    /// Solves A·X = B for a matrix expected to be symmetric positive definite, one solve per
    /// column of <paramref name="rightHandSides"/>.
    /// </summary>
    /// <remarks>
    /// Uses Cholesky, which is about twice as fast as LU and preserves symmetry, and falls back
    /// to LU if the matrix turns out not to be positive definite. The fallback matters for
    /// quantities that are positive definite in theory but can drift in floating point —
    /// a recursively updated covariance, for instance.
    /// </remarks>
    public static Matrix SolveSymmetricPositiveDefinite(this Matrix matrix, Matrix rightHandSides)
    {
        var cholesky = new CholeskyDecomposition(matrix);

        return cholesky.IsPositiveDefinite
            ? cholesky.Solve(rightHandSides)
            : new LuDecomposition(matrix).Solve(rightHandSides);
    }
}
