using CSharpNumerics.Physics.Constants;
using System;

namespace CSharpNumerics.Physics.Environmental.Fire;

/// <summary>
/// Heskestad's axisymmetric fire plume: flame height, centreline temperature and
/// velocity, and the rate at which the plume entrains air as it rises.
/// <para>
/// This is the correlation that says how fast smoke leaves a fire and how much
/// air it drags with it. It is what turns "smoke rises" into a rate — and
/// therefore what a stand-in constant rise velocity has to be replaced by before
/// any timing can be believed.
/// </para>
/// <para>
/// Valid above the flame tip for buoyancy-driven fires in an unconfined space.
/// It says nothing about a plume striking a deckhead, being deflected into a
/// ceiling jet, or filling a compartment — those are separate correlations.
/// </para>
/// <para>
/// Source: Heskestad, "Fire Plumes, Flame Height, and Air Entrainment",
/// SFPE Handbook of Fire Protection Engineering.
/// </para>
/// </summary>
public static class FirePlume
{
    /// <summary>Ambient air temperature used by default, in kelvin (20 °C).</summary>
    public const double DefaultAmbientTemperature = 293.15;

    /// <summary>Specific heat of air at constant pressure, in kJ/(kg·K).</summary>
    public const double AirSpecificHeat = 1.00;

    /// <summary>Ambient air density used by default, in kg/m³.</summary>
    public const double DefaultAmbientDensity = 1.2;

    // ═══════════════════════════════════════════════════════════════
    //  Flame geometry
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Mean flame height in metres above the fuel surface:
    /// <para>L = 0.235 · Q^(2/5) − 1.02 · D</para>
    /// Uses the total heat release rate, not the convective part. A negative
    /// result means the fire is too spread out to form a coherent flame and is
    /// returned as zero.
    /// </summary>
    /// <param name="totalKw">Total heat release rate in kW.</param>
    /// <param name="fireDiameter">Effective diameter of the fire base in metres.</param>
    public static double FlameHeight(double totalKw, double fireDiameter)
    {
        if (totalKw <= 0) return 0;
        if (fireDiameter <= 0) throw new ArgumentException("Fire diameter must be positive.", nameof(fireDiameter));

        double height = 0.235 * Math.Pow(totalKw, 0.4) - 1.02 * fireDiameter;
        return height > 0 ? height : 0;
    }

    /// <summary>
    /// Height of the virtual origin in metres, measured from the fuel surface:
    /// <para>z₀ = 0.083 · Q^(2/5) − 1.02 · D</para>
    /// The plume behaves as though it issued from a point source at this height,
    /// which may be below the fuel surface (negative) for a wide, weak fire.
    /// </summary>
    /// <param name="totalKw">Total heat release rate in kW.</param>
    /// <param name="fireDiameter">Effective diameter of the fire base in metres.</param>
    public static double VirtualOrigin(double totalKw, double fireDiameter)
    {
        if (totalKw <= 0) return 0;
        if (fireDiameter <= 0) throw new ArgumentException("Fire diameter must be positive.", nameof(fireDiameter));

        return 0.083 * Math.Pow(totalKw, 0.4) - 1.02 * fireDiameter;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Centreline conditions
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Centreline temperature rise above ambient in kelvin, at a height above
    /// the flame tip:
    /// <para>ΔT₀ = 9.1 · (T∞ / (g · c_p² · ρ∞²))^(1/3) · Qc^(2/3) · (z − z₀)^(−5/3)</para>
    /// For standard air this reduces to roughly 25 · Qc^(2/3) · (z − z₀)^(−5/3).
    /// Returns 0 at or below the virtual origin, where the correlation does not
    /// apply.
    /// </summary>
    /// <param name="convectiveKw">Convective heat release rate in kW.</param>
    /// <param name="height">Height above the fuel surface in metres.</param>
    /// <param name="virtualOrigin">Virtual origin height in metres, from <see cref="VirtualOrigin"/>.</param>
    /// <param name="ambientTemperature">Ambient temperature in kelvin.</param>
    /// <param name="ambientDensity">Ambient air density in kg/m³.</param>
    public static double CentrelineTemperatureRise(
        double convectiveKw,
        double height,
        double virtualOrigin,
        double ambientTemperature = DefaultAmbientTemperature,
        double ambientDensity = DefaultAmbientDensity)
    {
        if (convectiveKw <= 0) return 0;

        double z = height - virtualOrigin;
        if (z <= 0) return 0;

        double g = PhysicsConstants.GravitationalAcceleration;
        double factor = Math.Pow(
            ambientTemperature / (g * AirSpecificHeat * AirSpecificHeat * ambientDensity * ambientDensity),
            1.0 / 3.0);

        return 9.1 * factor * Math.Pow(convectiveKw, 2.0 / 3.0) * Math.Pow(z, -5.0 / 3.0);
    }

    /// <summary>
    /// Centreline vertical velocity in m/s, at a height above the flame tip:
    /// <para>u₀ = 3.4 · (g / (c_p · ρ∞ · T∞))^(1/3) · Qc^(1/3) · (z − z₀)^(−1/3)</para>
    /// For standard air this reduces to roughly 1.03 · (Qc / (z − z₀))^(1/3).
    /// <para>
    /// This is the quantity a buoyancy model wants: how fast the plume is
    /// actually moving upward at a given height above the fire.
    /// </para>
    /// </summary>
    /// <param name="convectiveKw">Convective heat release rate in kW.</param>
    /// <param name="height">Height above the fuel surface in metres.</param>
    /// <param name="virtualOrigin">Virtual origin height in metres, from <see cref="VirtualOrigin"/>.</param>
    /// <param name="ambientTemperature">Ambient temperature in kelvin.</param>
    /// <param name="ambientDensity">Ambient air density in kg/m³.</param>
    public static double CentrelineVelocity(
        double convectiveKw,
        double height,
        double virtualOrigin,
        double ambientTemperature = DefaultAmbientTemperature,
        double ambientDensity = DefaultAmbientDensity)
    {
        if (convectiveKw <= 0) return 0;

        double z = height - virtualOrigin;
        if (z <= 0) return 0;

        double g = PhysicsConstants.GravitationalAcceleration;
        double factor = Math.Pow(
            g / (AirSpecificHeat * ambientDensity * ambientTemperature),
            1.0 / 3.0);

        return 3.4 * factor * Math.Pow(convectiveKw, 1.0 / 3.0) * Math.Pow(z, -1.0 / 3.0);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Entrainment
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Mass flow rate in the plume in kg/s at a given height, above the flame tip:
    /// <para>ṁ = 0.071 · Qc^(1/3) · (z − z₀)^(5/3) + 0.0018 · Qc</para>
    /// <para>
    /// Almost all of this is entrained air rather than combustion products, which
    /// is why a compartment fills with smoke far faster than the fire produces
    /// it. Returns 0 at or below the virtual origin.
    /// </para>
    /// </summary>
    /// <param name="convectiveKw">Convective heat release rate in kW.</param>
    /// <param name="height">Height above the fuel surface in metres.</param>
    /// <param name="virtualOrigin">Virtual origin height in metres, from <see cref="VirtualOrigin"/>.</param>
    public static double MassFlowRate(double convectiveKw, double height, double virtualOrigin)
    {
        if (convectiveKw <= 0) return 0;

        double z = height - virtualOrigin;
        if (z <= 0) return 0;

        return 0.071 * Math.Pow(convectiveKw, 1.0 / 3.0) * Math.Pow(z, 5.0 / 3.0)
             + 0.0018 * convectiveKw;
    }
}
