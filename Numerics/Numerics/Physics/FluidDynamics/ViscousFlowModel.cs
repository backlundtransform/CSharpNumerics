using CSharpNumerics.Physics.FluidDynamics.Interfaces;

namespace CSharpNumerics.Physics.FluidDynamics;

/// <summary>
/// Default viscous pipe flow model using Hagen–Poiseuille physics.
/// </summary>
public class ViscousFlowModel : IViscousFlowModel
{
    public double DrivingForce(double pressureGradient, double density)
    {
        return -pressureGradient / density;
    }

    public double CylindricalDiffusion(double nu, double d2vdr2, double dvdr, double r)
    {
        return nu * (d2vdr2 + dvdr / r);
    }

    public double SymmetryAxisDiffusion(double nu, double d2vdr2)
    {
        return nu * 2.0 * d2vdr2;
    }
}
