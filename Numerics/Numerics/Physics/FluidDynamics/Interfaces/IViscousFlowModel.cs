namespace CSharpNumerics.Physics.FluidDynamics.Interfaces;

/// <summary>
/// Abstraction for viscous pipe flow physics used by pipe flow solvers.
/// Provides kinematic viscosity, driving force, and cylindrical
/// Laplacian diffusion computations.
/// </summary>
public interface IViscousFlowModel
{
    /// <summary>
    /// Computes the body-force driving term per unit mass from a pressure gradient:
    /// f = −(dP/dx) / ρ.
    /// </summary>
    double DrivingForce(double pressureGradient, double density);

    /// <summary>
    /// Computes the cylindrical Laplacian diffusion term for the radial
    /// momentum equation: ν · [d²v/dr² + (1/r)(dv/dr)].
    /// </summary>
    double CylindricalDiffusion(double nu, double d2vdr2, double dvdr, double r);

    /// <summary>
    /// Computes the cylindrical Laplacian at the symmetry axis (r → 0)
    /// using L'Hôpital's rule: ν · 2 · d²v/dr².
    /// </summary>
    double SymmetryAxisDiffusion(double nu, double d2vdr2);
}
