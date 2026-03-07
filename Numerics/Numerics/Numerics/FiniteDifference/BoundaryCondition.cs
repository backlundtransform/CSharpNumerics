namespace CSharpNumerics.Numerics.FiniteDifference
{
    /// <summary>
    /// Specifies how boundary cells are treated when computing finite-difference stencils.
    /// </summary>
    public enum BoundaryCondition
    {
        /// <summary>
        /// u = 0 outside the domain (fixed value).
        /// </summary>
        Dirichlet,

        /// <summary>
        /// ∂u/∂n = 0 at boundaries (zero-flux). Implemented by mirroring the nearest interior cell.
        /// </summary>
        Neumann,

        /// <summary>
        /// The domain wraps around: u[-1] = u[N-1], u[N] = u[0].
        /// </summary>
        Periodic
    }
}
