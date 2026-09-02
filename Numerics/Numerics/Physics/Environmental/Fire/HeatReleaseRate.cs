using CSharpNumerics.Physics.Environmental.Enums;
using System;

namespace CSharpNumerics.Physics.Environmental.Fire;

/// <summary>
/// Heat release rate of a compartment fire over time, using the t-squared growth
/// model that underpins NFPA 72, NFPA 92 and EN 1991-1-2.
/// <para>
/// Heat release rate is the single most important descriptor of a fire: plume
/// behaviour, hot gas layer temperature and the onset of flashover all follow
/// from it. Everything here is a pure function of time and a handful of fire
/// parameters — no grid, no simulation state.
/// </para>
/// </summary>
public static class HeatReleaseRate
{
    /// <summary>
    /// Reference heat release rate used to define the standard growth classes,
    /// in kW. Equal to 1000 BTU/s.
    /// </summary>
    public const double ReferenceHeatReleaseRate = 1055.0;

    /// <summary>
    /// Fraction of the total heat release carried by the plume as convected heat.
    /// The remainder leaves as radiation. Typically 0.6 to 0.8.
    /// </summary>
    public const double DefaultConvectiveFraction = 0.7;

    // ═══════════════════════════════════════════════════════════════
    //  Growth coefficients
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Time in seconds for a standard growth class to reach
    /// <see cref="ReferenceHeatReleaseRate"/>.
    /// </summary>
    public static double TimeToReference(FireGrowthRate rate)
    {
        switch (rate)
        {
            case FireGrowthRate.Slow: return 600.0;
            case FireGrowthRate.Medium: return 300.0;
            case FireGrowthRate.Fast: return 150.0;
            case FireGrowthRate.UltraFast: return 75.0;
            default: throw new ArgumentOutOfRangeException(nameof(rate));
        }
    }

    /// <summary>
    /// Growth coefficient α in kW/s² for a standard growth class:
    /// <para>α = 1055 / t_ref²</para>
    /// Gives 0.00293, 0.01172, 0.0469 and 0.1876 kW/s² for slow through
    /// ultra-fast.
    /// </summary>
    public static double GrowthCoefficient(FireGrowthRate rate)
    {
        double tRef = TimeToReference(rate);
        return ReferenceHeatReleaseRate / (tRef * tRef);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Growth phase
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Heat release rate in kW during unrestricted t-squared growth:
    /// <para>Q(t) = α · t²</para>
    /// Returns 0 before ignition.
    /// </summary>
    /// <param name="alpha">Growth coefficient in kW/s².</param>
    /// <param name="timeSeconds">Time since ignition in seconds.</param>
    public static double Growth(double alpha, double timeSeconds)
    {
        if (alpha < 0) throw new ArgumentException("Growth coefficient must be non-negative.", nameof(alpha));
        if (timeSeconds <= 0) return 0;
        return alpha * timeSeconds * timeSeconds;
    }

    /// <summary>
    /// Time in seconds for a t-squared fire to reach a given heat release rate:
    /// <para>t = √(Q / α)</para>
    /// </summary>
    /// <param name="alpha">Growth coefficient in kW/s².</param>
    /// <param name="targetKw">Target heat release rate in kW.</param>
    public static double TimeToReach(double alpha, double targetKw)
    {
        if (alpha <= 0) throw new ArgumentException("Growth coefficient must be positive.", nameof(alpha));
        if (targetKw < 0) throw new ArgumentException("Target must be non-negative.", nameof(targetKw));
        return Math.Sqrt(targetKw / alpha);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Full design curve
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Heat release rate in kW over the whole life of a design fire: t-squared
    /// growth up to <paramref name="peakKw"/>, a steady burning period, then
    /// linear decay to zero.
    /// <para>
    /// This is the shape used for design fires in NFPA 92 and the fire
    /// engineering literature. The steady period represents burning limited by
    /// fuel surface area or ventilation rather than by growth.
    /// </para>
    /// </summary>
    /// <param name="alpha">Growth coefficient in kW/s².</param>
    /// <param name="peakKw">Peak heat release rate in kW.</param>
    /// <param name="steadyDurationSeconds">How long the fire burns at its peak, in seconds.</param>
    /// <param name="decayDurationSeconds">How long decay from peak to zero takes, in seconds.</param>
    /// <param name="timeSeconds">Time since ignition in seconds.</param>
    public static double Curve(
        double alpha,
        double peakKw,
        double steadyDurationSeconds,
        double decayDurationSeconds,
        double timeSeconds)
    {
        if (alpha <= 0) throw new ArgumentException("Growth coefficient must be positive.", nameof(alpha));
        if (peakKw <= 0) throw new ArgumentException("Peak must be positive.", nameof(peakKw));
        if (steadyDurationSeconds < 0) throw new ArgumentException("Steady duration must be non-negative.", nameof(steadyDurationSeconds));
        if (decayDurationSeconds < 0) throw new ArgumentException("Decay duration must be non-negative.", nameof(decayDurationSeconds));

        if (timeSeconds <= 0) return 0;

        double growthEnd = TimeToReach(alpha, peakKw);
        if (timeSeconds < growthEnd)
            return alpha * timeSeconds * timeSeconds;

        double steadyEnd = growthEnd + steadyDurationSeconds;
        if (timeSeconds <= steadyEnd)
            return peakKw;

        if (decayDurationSeconds <= 0) return 0;

        double intoDecay = timeSeconds - steadyEnd;
        if (intoDecay >= decayDurationSeconds) return 0;

        return peakKw * (1.0 - intoDecay / decayDurationSeconds);
    }

    /// <summary>
    /// Convective part of the heat release rate in kW — the fraction that goes
    /// into the plume and drives buoyancy. Plume correlations take this rather
    /// than the total.
    /// </summary>
    /// <param name="totalKw">Total heat release rate in kW.</param>
    /// <param name="convectiveFraction">Convected fraction, typically 0.6 to 0.8.</param>
    public static double Convective(double totalKw, double convectiveFraction = DefaultConvectiveFraction)
    {
        if (convectiveFraction < 0 || convectiveFraction > 1)
            throw new ArgumentOutOfRangeException(nameof(convectiveFraction), "Convective fraction must be between 0 and 1.");
        return totalKw <= 0 ? 0 : totalKw * convectiveFraction;
    }
}
