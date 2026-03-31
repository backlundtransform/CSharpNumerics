using System;

namespace CSharpNumerics.Physics.SolidMechanics;

/// <summary>
/// Extension methods for stress and strain calculations:
/// Hooke's law, normal/shear stress, and shear modulus.
/// </summary>
public static class StressStrainExtensions
{
    // ═══════════════════════════════════════════════════════════════
    //  Stress & Strain — Hooke's Law
    // ═══════════════════════════════════════════════════════════════

    #region Stress & Strain

    /// <summary>
    /// Normal (axial) stress: σ = F / A.
    /// </summary>
    /// <param name="force">Applied force F in newtons.</param>
    /// <param name="area">Cross-sectional area A in m².</param>
    /// <returns>Normal stress σ in Pa.</returns>
    public static double NormalStress(this double force, double area)
    {
        if (area <= 0) throw new ArgumentException("Area must be greater than zero.");
        return force / area;
    }

    /// <summary>
    /// Normal strain from Hooke's law: ε = σ / E.
    /// </summary>
    /// <param name="stress">Normal stress σ in Pa.</param>
    /// <param name="youngsModulus">Young's modulus E in Pa.</param>
    /// <returns>Dimensionless strain ε.</returns>
    public static double NormalStrain(this double stress, double youngsModulus)
    {
        if (youngsModulus <= 0) throw new ArgumentException("Young's modulus must be greater than zero.");
        return stress / youngsModulus;
    }

    /// <summary>
    /// Hooke's law: σ = E · ε.
    /// </summary>
    /// <param name="youngsModulus">Young's modulus E in Pa.</param>
    /// <param name="strain">Dimensionless strain ε.</param>
    /// <returns>Normal stress σ in Pa.</returns>
    public static double HookesLaw(this double youngsModulus, double strain)
    {
        return youngsModulus * strain;
    }

    /// <summary>
    /// Average shear stress: τ = V / A.
    /// </summary>
    /// <param name="shearForce">Transverse shear force V in newtons.</param>
    /// <param name="area">Cross-sectional area A in m².</param>
    /// <returns>Shear stress τ in Pa.</returns>
    public static double ShearStress(this double shearForce, double area)
    {
        if (area <= 0) throw new ArgumentException("Area must be greater than zero.");
        return shearForce / area;
    }

    /// <summary>
    /// Shear modulus from Young's modulus and Poisson's ratio: G = E / (2(1 + ν)).
    /// </summary>
    /// <param name="youngsModulus">Young's modulus E in Pa.</param>
    /// <param name="poissonsRatio">Poisson's ratio ν (dimensionless, typically 0–0.5).</param>
    /// <returns>Shear modulus G in Pa.</returns>
    public static double ShearModulus(this double youngsModulus, double poissonsRatio)
    {
        return youngsModulus / (2.0 * (1.0 + poissonsRatio));
    }

    #endregion
}
