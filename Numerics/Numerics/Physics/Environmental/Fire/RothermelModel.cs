using System;
using CSharpNumerics.Physics.Materials.Fire;

namespace CSharpNumerics.Physics.Environmental.Fire;

/// <summary>
/// Rothermel (1972) surface fire spread model — the standard used by
/// FARSITE, FlamMap, and BehavePlus.
/// <para>
/// All methods are pure functions operating on <see cref="FuelModel"/> parameters.
/// No grid or simulation state — this is the physics layer.
/// </para>
/// </summary>
public static class RothermelModel
{
    // ═══════════════════════════════════════════════════════════════
    //  Primary equation: Rate of Spread
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Computes the Rothermel rate of spread (m/min):
    /// <para>R = IR · ξ · (1 + φw + φs) / (ρb · ε · Qig)</para>
    /// Returns 0 when moisture ≥ moisture of extinction.
    /// </summary>
    /// <param name="fuel">Fuel model parameters.</param>
    /// <param name="moistureContent">Dead fuel moisture content (fraction, 0–1).</param>
    /// <param name="midflameWindSpeed">Mid-flame wind speed (m/s).</param>
    /// <param name="slopeRadians">Terrain slope angle (radians, 0 = flat).</param>
    public static double RateOfSpread(FuelModel fuel, double moistureContent, double midflameWindSpeed, double slopeRadians)
    {
        if (fuel.OvendryFuelLoad <= 0 || fuel.FuelBedDepth <= 0)
            return 0;

        // Moisture ratio — if M >= Mx, fire cannot spread
        if (moistureContent >= fuel.MoistureOfExtinction)
            return 0;

        double IR = ReactionIntensity(fuel, moistureContent);
        double xi = PropagatingFluxRatio(fuel);
        double phiW = WindFactor(fuel, midflameWindSpeed);
        double phiS = SlopeFactor(PackingRatio(fuel), slopeRadians);

        double rhoB = BulkDensity(fuel);
        double eps = EffectiveHeatingNumber(fuel);
        double Qig = HeatOfPreignition(moistureContent);

        double denominator = rhoB * eps * Qig;
        if (denominator <= 0)
            return 0;

        return IR * xi * (1.0 + phiW + phiS) / denominator;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Sub-equations
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Reaction intensity IR (kJ/m²·min) — the rate of energy release per unit area
    /// of the flaming fire front.
    /// <para>
    /// IR = Γ' · wn · h · ηM · ηs
    /// </para>
    /// </summary>
    public static double ReactionIntensity(FuelModel fuel, double moistureContent)
    {
        if (fuel.OvendryFuelLoad <= 0 || fuel.FuelBedDepth <= 0)
            return 0;

        double beta = PackingRatio(fuel);
        double betaOp = OptimalPackingRatio(fuel);
        double sigma = fuel.SurfaceAreaToVolumeRatio;

        // Optimum reaction velocity Γ'max (1/min)
        // Γ'max = σ^1.5 / (495 + 0.0594·σ^1.5)
        double sigma15 = Math.Pow(sigma, 1.5);
        double gammaMax = sigma15 / (495.0 + 0.0594 * sigma15);

        // Γ' = Γ'max · (β/β_op)^A · exp(A·(1 - β/β_op))
        // A = 1 / (4.774·σ^0.1 - 7.27)
        double A = 1.0 / (4.774 * Math.Pow(sigma, 0.1) - 7.27);
        double ratio = beta / betaOp;
        double gamma = gammaMax * Math.Pow(ratio, A) * Math.Exp(A * (1.0 - ratio));

        // Net fuel load wn = w₀ · (1 - ST) where ST = total mineral content ≈ 0.0555
        double ST = 0.0555;
        double wn = fuel.OvendryFuelLoad * (1.0 - ST);

        // Moisture damping ηM = 1 - 2.59·rM + 5.11·rM² - 3.52·rM³
        // rM = M / Mx  (moisture ratio)
        double rM = moistureContent / fuel.MoistureOfExtinction;
        rM = Math.Min(rM, 1.0);
        double etaM = 1.0 - 2.59 * rM + 5.11 * rM * rM - 3.52 * rM * rM * rM;
        etaM = Math.Max(etaM, 0);

        // Mineral damping ηs = 0.174·Se^(-0.19)  where Se = effective mineral content ≈ 0.01
        double Se = 0.01;
        double etaS = 0.174 * Math.Pow(Se, -0.19);
        etaS = Math.Min(etaS, 1.0);

        return gamma * wn * fuel.LowHeatContent * etaM * etaS;
    }

    /// <summary>
    /// Wind correction factor φw (dimensionless).
    /// <para>
    /// φw = C · (β/β_op)^(-E) · U^B
    /// </para>
    /// where U is mid-flame wind speed in m/min.
    /// </summary>
    /// <param name="fuel">Fuel model parameters.</param>
    /// <param name="midflameWindSpeed">Mid-flame wind speed (m/s).</param>
    public static double WindFactor(FuelModel fuel, double midflameWindSpeed)
    {
        if (midflameWindSpeed <= 0)
            return 0;

        double sigma = fuel.SurfaceAreaToVolumeRatio;
        double beta = PackingRatio(fuel);
        double betaOp = OptimalPackingRatio(fuel);

        // C = 7.47 · exp(-0.133 · σ^0.55)
        double C = 7.47 * Math.Exp(-0.133 * Math.Pow(sigma, 0.55));

        // B = 0.02526 · σ^0.54
        double B = 0.02526 * Math.Pow(sigma, 0.54);

        // E = 0.715 · exp(-3.59e-4 · σ)
        double E = 0.715 * Math.Exp(-3.59e-4 * sigma);

        // Convert wind speed m/s → m/min (Rothermel uses m/min internally)
        double U = midflameWindSpeed * 60.0;

        return C * Math.Pow(beta / betaOp, -E) * Math.Pow(U, B);
    }

    /// <summary>
    /// Slope correction factor φs (dimensionless).
    /// <para>φs = 5.275 · β^(-0.3) · tan²(θ)</para>
    /// Returns 0 for flat or downhill slopes.
    /// </summary>
    /// <param name="packingRatio">Packing ratio β = ρb / ρp.</param>
    /// <param name="slopeRadians">Terrain slope angle in radians.</param>
    public static double SlopeFactor(double packingRatio, double slopeRadians)
    {
        if (slopeRadians <= 0 || packingRatio <= 0)
            return 0;

        double tanSlope = Math.Tan(slopeRadians);
        return 5.275 * Math.Pow(packingRatio, -0.3) * tanSlope * tanSlope;
    }

    /// <summary>
    /// Propagating flux ratio ξ (dimensionless) — fraction of reaction intensity
    /// that heats adjacent fuel particles.
    /// <para>ξ = exp((0.792 + 0.681·σ^0.5)·(β + 0.1)) / (192 + 0.2595·σ)</para>
    /// </summary>
    public static double PropagatingFluxRatio(FuelModel fuel)
    {
        double sigma = fuel.SurfaceAreaToVolumeRatio;
        double beta = PackingRatio(fuel);

        double numerator = Math.Exp((0.792 + 0.681 * Math.Sqrt(sigma)) * (beta + 0.1));
        double denominator = 192.0 + 0.2595 * sigma;

        return numerator / denominator;
    }

    /// <summary>
    /// Heat of pre-ignition Qig (kJ/kg) — energy required to bring fuel to ignition.
    /// <para>Qig = 581 + 2594 · Mf</para>
    /// where Mf is the fuel moisture content (fraction).
    /// </summary>
    public static double HeatOfPreignition(double moistureContent)
    {
        return 581.0 + 2594.0 * moistureContent;
    }

    /// <summary>
    /// Effective heating number ε (dimensionless) — fraction of fuel that
    /// participates in ignition.
    /// <para>ε = exp(-138 / σ)</para>
    /// </summary>
    public static double EffectiveHeatingNumber(FuelModel fuel)
    {
        if (fuel.SurfaceAreaToVolumeRatio <= 0) return 0;
        return Math.Exp(-138.0 / fuel.SurfaceAreaToVolumeRatio);
    }

    /// <summary>
    /// Packing ratio β = ρb / ρp (dimensionless).
    /// ρb = w₀ / δ (ovendry bulk density).
    /// </summary>
    public static double PackingRatio(FuelModel fuel)
    {
        if (fuel.FuelBedDepth <= 0 || fuel.ParticleDensity <= 0)
            return 0;

        double rhoB = fuel.OvendryFuelLoad / fuel.FuelBedDepth;
        return rhoB / fuel.ParticleDensity;
    }

    /// <summary>
    /// Optimal packing ratio β_op = 3.348 · σ^(-0.8189).
    /// </summary>
    public static double OptimalPackingRatio(FuelModel fuel)
    {
        if (fuel.SurfaceAreaToVolumeRatio <= 0) return 0;
        return 3.348 * Math.Pow(fuel.SurfaceAreaToVolumeRatio, -0.8189);
    }

    /// <summary>
    /// Ovendry bulk density ρb = w₀ / δ (kg/m³).
    /// </summary>
    public static double BulkDensity(FuelModel fuel)
    {
        if (fuel.FuelBedDepth <= 0) return 0;
        return fuel.OvendryFuelLoad / fuel.FuelBedDepth;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Derived outputs
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Byram's fireline intensity (kW/m) and the correlated flame length (m).
    /// <para>
    /// I_B = IR · R / 60  (kW/m, with R in m/min and IR in kJ/m²·min)
    /// </para>
    /// <para>
    /// L = 0.0775 · I_B^0.46  (Byram 1959 flame length correlation, metres)
    /// </para>
    /// </summary>
    /// <param name="reactionIntensity">IR in kJ/m²·min.</param>
    /// <param name="rateOfSpread">R in m/min.</param>
    /// <returns>Flame length in metres.</returns>
    public static double FlameLength(double reactionIntensity, double rateOfSpread)
    {
        if (reactionIntensity <= 0 || rateOfSpread <= 0) return 0;

        // Byram's fireline intensity: IB = IR * R / 60 (kW/m)
        // Note: IR is kJ/(m²·min), R is m/min → IR*R is kJ/(m·min) → /60 gives kW/m
        double IB = reactionIntensity * rateOfSpread / 60.0;

        // Flame length (m) — Byram (1959) correlation
        return 0.0775 * Math.Pow(IB, 0.46);
    }
}
