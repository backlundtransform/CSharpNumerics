using System;

namespace CSharpNumerics.Physics.Environmental.Water;

/// <summary>
/// Longitudinal dispersion in rivers and open channels — Fischer et al. (1979)
/// dispersion coefficient, shear velocity, first-order decay, and retardation.
/// <para>
/// All methods are pure functions. No grid or simulation state.
/// </para>
/// </summary>
public static class LongitudinalDispersion
{
    private const double Gravity = 9.80665; // m/s²

    // ═══════════════════════════════════════════════════════════════
    //  Fischer dispersion coefficient
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Fischer (1979) longitudinal dispersion coefficient EL (m²/s):
    /// <para>EL = 0.011 · u² · W² / (H · u*)</para>
    /// where u is mean velocity, W is channel width, H is mean depth,
    /// and u* is shear velocity.
    /// Returns 0 if any input is non-positive.
    /// </summary>
    /// <param name="velocity">Cross-sectional mean velocity u (m/s).</param>
    /// <param name="width">Channel width W (m).</param>
    /// <param name="depth">Mean flow depth H (m).</param>
    /// <param name="shearVelocity">Shear velocity u* (m/s).</param>
    public static double FischerCoefficient(double velocity, double width, double depth, double shearVelocity)
    {
        if (velocity <= 0 || width <= 0 || depth <= 0 || shearVelocity <= 0)
            return 0;

        return 0.011 * velocity * velocity * width * width / (depth * shearVelocity);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Shear velocity
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Shear velocity u* (m/s) = √(g · Rh · S₀).
    /// </summary>
    /// <param name="hydraulicRadius">Hydraulic radius Rh (m).</param>
    /// <param name="bedSlope">Bed slope S₀ (m/m).</param>
    public static double ShearVelocity(double hydraulicRadius, double bedSlope)
    {
        if (hydraulicRadius <= 0 || bedSlope <= 0)
            return 0;

        return Math.Sqrt(Gravity * hydraulicRadius * bedSlope);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Decay & retardation
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// First-order decay constant λ (1/s) from half-life:
    /// <para>λ = ln(2) / t½</para>
    /// Returns 0 if half-life ≤ 0 or infinite (conservative tracer).
    /// </summary>
    /// <param name="halfLifeSeconds">Half-life t½ (s).</param>
    public static double DecayConstant(double halfLifeSeconds)
    {
        if (halfLifeSeconds <= 0 || double.IsPositiveInfinity(halfLifeSeconds))
            return 0;

        return Math.Log(2) / halfLifeSeconds;
    }

    /// <summary>
    /// Retardation factor Rf (dimensionless):
    /// <para>Rf = 1 + ρs · Kd / n</para>
    /// where ρs is bulk sediment density, Kd is partition coefficient,
    /// and n is porosity. Rf = 1 means no retardation (conservative tracer).
    /// Returns 1 if Kd ≤ 0.
    /// </summary>
    /// <param name="bulkDensity">Bulk sediment density ρs (kg/m³). Typical: 1400–1800.</param>
    /// <param name="partitionCoefficient">Solid–water partition coefficient Kd (L/kg).</param>
    /// <param name="porosity">Bed sediment porosity n (fraction, 0–1). Typical: 0.3–0.5.</param>
    public static double RetardationFactor(double bulkDensity, double partitionCoefficient, double porosity)
    {
        if (partitionCoefficient <= 0 || porosity <= 0 || bulkDensity <= 0)
            return 1.0;

        // Kd is in L/kg = 10⁻³ m³/kg, but the standard formula uses consistent units:
        // Rf = 1 + (ρs/n) · Kd  where Kd in L/kg and ρs in kg/L → convert ρs from kg/m³ to kg/L
        double rhoS_kgPerL = bulkDensity / 1000.0;
        return 1.0 + rhoS_kgPerL * partitionCoefficient / porosity;
    }
}
