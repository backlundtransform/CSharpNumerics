using System;

namespace CSharpNumerics.Physics.FluidDynamics;

/// <summary>
/// Extension methods for Kármán vortex street analysis:
/// vortex shedding frequency, Strouhal–Reynolds correlations,
/// Roshko number, vortex spacing, convection velocity,
/// and regime classification.
/// </summary>
public static class KarmanVortexStreetExtensions
{
    // ═══════════════════════════════════════════════════════════════
    //  Von Kármán stability ratio
    //
    //  h/a = (1/π) · arccosh(√2) ≈ 0.2806
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Von Kármán's theoretical spacing ratio h/a = (1/π)·arccosh(√2)
    /// for a stable staggered vortex street.
    /// </summary>
    public static readonly double StableSpacingRatio =
        Math.Log(Math.Sqrt(2) + 1) / Math.PI;   // arccosh(√2) = ln(√2 + 1)

    // ═══════════════════════════════════════════════════════════════
    //  Vortex shedding frequency
    //
    //  f = St · U / D
    // ═══════════════════════════════════════════════════════════════

    #region Shedding Frequency

    /// <summary>
    /// Computes the vortex shedding frequency from the Strouhal number:
    /// f = St · U / D.
    /// </summary>
    /// <param name="strouhalNumber">Strouhal number St (dimensionless).</param>
    /// <param name="freeStreamSpeed">Free-stream velocity U in m/s.</param>
    /// <param name="diameter">Cylinder (or body) diameter D in m.</param>
    public static double VortexSheddingFrequency(
        this double strouhalNumber, double freeStreamSpeed, double diameter)
    {
        if (diameter <= 0) throw new ArgumentException("Diameter must be positive.");
        if (freeStreamSpeed <= 0) throw new ArgumentException("Free-stream speed must be positive.");
        return strouhalNumber * freeStreamSpeed / diameter;
    }

    /// <summary>
    /// Computes the vortex shedding frequency using the Roshko empirical
    /// Strouhal–Reynolds correlation for a circular cylinder:
    /// f = St(Re) · U / D, where St ≈ 0.198 · (1 − 19.7 / Re).
    /// Valid for 250 &lt; Re &lt; 2 × 10⁵.
    /// </summary>
    /// <param name="reynoldsNumber">Reynolds number Re (dimensionless).</param>
    /// <param name="freeStreamSpeed">Free-stream velocity U in m/s.</param>
    /// <param name="diameter">Cylinder diameter D in m.</param>
    public static double VortexSheddingFrequencyFromReynolds(
        this double reynoldsNumber, double freeStreamSpeed, double diameter)
    {
        double st = reynoldsNumber.RoshkoStrouhalNumber();
        return st.VortexSheddingFrequency(freeStreamSpeed, diameter);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Strouhal–Reynolds correlation (Roshko 1954)
    //
    //  St ≈ 0.198 (1 − 19.7 / Re)   for 250 < Re < 2×10⁵
    // ═══════════════════════════════════════════════════════════════

    #region Strouhal–Reynolds Correlation

    /// <summary>
    /// Computes the Strouhal number for a circular cylinder using
    /// Roshko's empirical correlation (1954):
    /// St ≈ 0.198 · (1 − 19.7 / Re).
    /// Valid for 250 &lt; Re &lt; 2 × 10⁵.
    /// </summary>
    /// <param name="reynoldsNumber">Reynolds number Re (dimensionless).</param>
    public static double RoshkoStrouhalNumber(this double reynoldsNumber)
    {
        if (reynoldsNumber <= 0)
            throw new ArgumentException("Reynolds number must be positive.");
        return 0.198 * (1.0 - 19.7 / reynoldsNumber);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Roshko number
    //
    //  Ro = f · D² / ν = St · Re
    // ═══════════════════════════════════════════════════════════════

    #region Roshko Number

    /// <summary>
    /// Computes the Roshko number: Ro = f · D² / ν = St · Re.
    /// Combines the shedding frequency with viscous scaling.
    /// </summary>
    /// <param name="sheddingFrequency">Vortex shedding frequency f in Hz.</param>
    /// <param name="diameter">Cylinder diameter D in m.</param>
    /// <param name="kinematicViscosity">Kinematic viscosity ν in m²/s.</param>
    public static double RoshkoNumber(
        this double sheddingFrequency, double diameter, double kinematicViscosity)
    {
        if (kinematicViscosity <= 0)
            throw new ArgumentException("Kinematic viscosity must be positive.");
        return sheddingFrequency * diameter * diameter / kinematicViscosity;
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Vortex street geometry
    //
    //  a = U / f         (streamwise spacing)
    //  h = (h/a) · a     (lateral spacing, h/a ≈ 0.281)
    // ═══════════════════════════════════════════════════════════════

    #region Vortex Street Geometry

    /// <summary>
    /// Computes the streamwise (longitudinal) spacing between consecutive
    /// vortices of the same row: a = U / f.
    /// </summary>
    /// <param name="freeStreamSpeed">Free-stream velocity U in m/s.</param>
    /// <param name="sheddingFrequency">Vortex shedding frequency f in Hz.</param>
    public static double VortexStreamwiseSpacing(
        this double freeStreamSpeed, double sheddingFrequency)
    {
        if (sheddingFrequency <= 0)
            throw new ArgumentException("Shedding frequency must be positive.");
        return freeStreamSpeed / sheddingFrequency;
    }

    /// <summary>
    /// Computes the lateral (cross-stream) spacing between the two rows
    /// of vortices using von Kármán's stability ratio:
    /// h = (h/a) · a ≈ 0.281 · a.
    /// </summary>
    /// <param name="streamwiseSpacing">Streamwise spacing a in m.</param>
    public static double VortexLateralSpacing(this double streamwiseSpacing)
    {
        return StableSpacingRatio * streamwiseSpacing;
    }

    /// <summary>
    /// Computes the wavelength of the vortex street, which equals the
    /// streamwise spacing: λ = U / f.
    /// </summary>
    /// <param name="freeStreamSpeed">Free-stream velocity U in m/s.</param>
    /// <param name="sheddingFrequency">Vortex shedding frequency f in Hz.</param>
    public static double VortexStreetWavelength(
        this double freeStreamSpeed, double sheddingFrequency)
    {
        return freeStreamSpeed.VortexStreamwiseSpacing(sheddingFrequency);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Convection velocity
    //
    //  U_v ≈ (1 − 1/Re^0.5 · k) · U   (empirical ≈ 0.85–0.90 U)
    // ═══════════════════════════════════════════════════════════════

    #region Convection Velocity

    /// <summary>
    /// Estimates the downstream convection velocity of shed vortices.
    /// Empirically, vortices travel at approximately 85–90% of the
    /// free-stream velocity: U_v ≈ convectionRatio · U.
    /// </summary>
    /// <param name="freeStreamSpeed">Free-stream velocity U in m/s.</param>
    /// <param name="convectionRatio">Fraction of U (default 0.875, typical for moderate Re).</param>
    public static double VortexConvectionVelocity(
        this double freeStreamSpeed, double convectionRatio = 0.875)
    {
        return convectionRatio * freeStreamSpeed;
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Shedding regime classification
    // ═══════════════════════════════════════════════════════════════

    #region Regime Classification

    /// <summary>
    /// Returns true if the Reynolds number falls within the periodic
    /// vortex shedding regime for a circular cylinder (≈ 47 &lt; Re &lt; 2 × 10⁵).
    /// </summary>
    /// <param name="reynoldsNumber">Reynolds number Re (dimensionless).</param>
    public static bool IsPeriodicSheddingRegime(this double reynoldsNumber)
    {
        return reynoldsNumber > 47.0 && reynoldsNumber < 2e5;
    }

    /// <summary>
    /// Classifies the cylinder wake regime based on the Reynolds number.
    /// Returns a descriptive string: "Creeping", "Steady separation",
    /// "Laminar shedding", "Turbulent transition", or "Turbulent wake".
    /// </summary>
    /// <param name="reynoldsNumber">Reynolds number Re (dimensionless).</param>
    public static string CylinderWakeRegime(this double reynoldsNumber)
    {
        if (reynoldsNumber < 5) return "Creeping";
        if (reynoldsNumber < 47) return "Steady separation";
        if (reynoldsNumber < 200) return "Laminar shedding";
        if (reynoldsNumber < 2e5) return "Turbulent transition";
        return "Turbulent wake";
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Drag coefficient from vortex street (von Kármán)
    //
    //  Cd = (2h/a) · ((U_v/U)² − 1) + (2h/a)·(h/(π·a))
    //  Simplified empirical: Cd ≈ f(Re) for a circular cylinder
    // ═══════════════════════════════════════════════════════════════

    #region Drag from Vortex Street

    /// <summary>
    /// Computes the drag coefficient from von Kármán's vortex street
    /// momentum analysis:
    /// C_d = (2h / (a·D)) · [U∞² − U_v²] / U∞²   (per unit span).
    /// This is the idealised inviscid wake-drag formula; real drag is higher
    /// due to viscous contributions.
    /// </summary>
    /// <param name="lateralSpacing">Lateral spacing h between vortex rows in m.</param>
    /// <param name="streamwiseSpacing">Streamwise spacing a between consecutive vortices in m.</param>
    /// <param name="freeStreamSpeed">Free-stream velocity U in m/s.</param>
    /// <param name="convectionVelocity">Vortex convection velocity U_v in m/s.</param>
    /// <param name="diameter">Cylinder diameter D in m.</param>
    public static double VortexStreetDragCoefficient(
        this double lateralSpacing, double streamwiseSpacing,
        double freeStreamSpeed, double convectionVelocity, double diameter)
    {
        if (streamwiseSpacing <= 0 || diameter <= 0 || freeStreamSpeed <= 0)
            throw new ArgumentException("Spacing, diameter, and speed must be positive.");

        double ratio = lateralSpacing / streamwiseSpacing;
        double speedRatio = convectionVelocity / freeStreamSpeed;
        return 2.0 * ratio * (1.0 - speedRatio * speedRatio)
             + 2.0 * ratio * ratio * lateralSpacing / (Math.PI * streamwiseSpacing);
    }

    #endregion
}
