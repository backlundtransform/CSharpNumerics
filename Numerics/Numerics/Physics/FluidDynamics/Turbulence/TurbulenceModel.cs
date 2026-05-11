using System;

namespace CSharpNumerics.Physics.FluidDynamics.Turbulence;

/// <summary>
/// Smagorinsky subgrid-scale (SGS) turbulence model for Large Eddy Simulation (LES).
/// Computes an eddy (turbulent) viscosity from the local strain rate:
///   ν_t = (Cs · Δ)² · |S|
/// where Cs is the Smagorinsky coefficient, Δ is the grid spacing (filter width),
/// and |S| = √(2 · Sij · Sij) is the magnitude of the strain-rate tensor.
/// 
/// This is the simplest and most widely used SGS model, ideal for game-quality
/// visual turbulence where computational cost matters more than quantitative accuracy.
/// </summary>
public class TurbulenceModel
{
    /// <summary>Smagorinsky coefficient. Typical range: 0.1–0.2. Default 0.17.</summary>
    public double Cs { get; set; } = 0.17;

    /// <summary>Grid spacing / filter width in metres.</summary>
    public double Delta { get; set; }

    /// <summary>
    /// Creates a Smagorinsky turbulence model.
    /// </summary>
    /// <param name="gridSpacing">Grid spacing Δ in metres.</param>
    /// <param name="cs">Smagorinsky coefficient (default 0.17).</param>
    public TurbulenceModel(double gridSpacing, double cs = 0.17)
    {
        Delta = gridSpacing;
        Cs = cs;
    }

    /// <summary>
    /// Computes the turbulent (eddy) viscosity ν_t = (Cs·Δ)²·|S|.
    /// </summary>
    /// <param name="strainRateMagnitude">|S| = √(2·Sij·Sij), the magnitude of the strain-rate tensor.</param>
    /// <returns>Turbulent kinematic viscosity in m²/s.</returns>
    public double EddyViscosity(double strainRateMagnitude)
    {
        double csDelta = Cs * Delta;
        return csDelta * csDelta * Math.Abs(strainRateMagnitude);
    }

    /// <summary>
    /// Computes the effective (total) viscosity: ν_eff = ν + ν_t.
    /// </summary>
    /// <param name="molecularViscosity">Molecular kinematic viscosity ν in m²/s.</param>
    /// <param name="strainRateMagnitude">|S|, magnitude of strain-rate tensor.</param>
    /// <returns>Total effective viscosity in m²/s.</returns>
    public double EffectiveViscosity(double molecularViscosity, double strainRateMagnitude)
    {
        return molecularViscosity + EddyViscosity(strainRateMagnitude);
    }

    /// <summary>
    /// Computes the strain-rate magnitude |S| from the 2D velocity gradients on a staggered grid.
    /// |S| = √(2·(S11² + S22² + 2·S12²))
    /// where S11 = du/dx, S22 = dv/dy, S12 = ½(du/dy + dv/dx).
    /// </summary>
    /// <param name="dudx">∂u/∂x</param>
    /// <param name="dudy">∂u/∂y</param>
    /// <param name="dvdx">∂v/∂x</param>
    /// <param name="dvdy">∂v/∂y</param>
    public static double StrainRate2D(double dudx, double dudy, double dvdx, double dvdy)
    {
        double s11 = dudx;
        double s22 = dvdy;
        double s12 = 0.5 * (dudy + dvdx);
        return Math.Sqrt(2.0 * (s11 * s11 + s22 * s22 + 2.0 * s12 * s12));
    }

    /// <summary>
    /// Computes the strain-rate magnitude |S| from 3D velocity gradients.
    /// |S| = √(2·(S11² + S22² + S33² + 2·S12² + 2·S13² + 2·S23²))
    /// </summary>
    public static double StrainRate3D(
        double dudx, double dudy, double dudz,
        double dvdx, double dvdy, double dvdz,
        double dwdx, double dwdy, double dwdz)
    {
        double s11 = dudx;
        double s22 = dvdy;
        double s33 = dwdz;
        double s12 = 0.5 * (dudy + dvdx);
        double s13 = 0.5 * (dudz + dwdx);
        double s23 = 0.5 * (dvdz + dwdy);
        return Math.Sqrt(2.0 * (s11 * s11 + s22 * s22 + s33 * s33
            + 2.0 * s12 * s12 + 2.0 * s13 * s13 + 2.0 * s23 * s23));
    }
}
