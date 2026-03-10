namespace CSharpNumerics.Physics.Waves
{
    /// <summary>
    /// Boundary condition types for wave simulations.
    /// Maps to <see cref="CSharpNumerics.Numerics.FiniteDifference.BoundaryCondition"/>
    /// for finite-difference operators, with the addition of absorbing boundaries.
    /// </summary>
    public enum BoundaryType
    {
        /// <summary>u = 0 at boundary (fixed end, e.g. clamped string).</summary>
        Fixed,

        /// <summary>∂u/∂n = 0 at boundary (free end, zero slope).</summary>
        Free,

        /// <summary>Domain wraps around: u(0) = u(L).</summary>
        Periodic,

        /// <summary>Waves are absorbed at the boundary (no reflection).</summary>
        Absorbing
    }
}
