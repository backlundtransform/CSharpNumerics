using System;

namespace CSharpNumerics.Physics.Environmental.Fire;

/// <summary>
/// Purser's fractional effective dose model for the toxic and obscuring effects
/// of fire smoke, plus visibility through it.
/// <para>
/// This is what turns a concentration field into an answer someone can act on.
/// Fire rarely kills by heat: it incapacitates through carbon monoxide, hydrogen
/// cyanide and oxygen depletion, and it prevents escape by destroying visibility
/// long before that. A dose reaching 1.0 marks incapacitation of a susceptible
/// person, at which point self-rescue has stopped.
/// </para>
/// <para>
/// Doses accumulate over time, so the methods here give a rate per minute of
/// exposure at a given concentration. Integrate them along a time history to get
/// the dose. Concentrations are in parts per million by volume unless stated.
/// </para>
/// <para>
/// Source: Purser, "Assessment of Hazards to Occupants from Smoke, Toxic Gases
/// and Heat", SFPE Handbook of Fire Protection Engineering; ISO 13571.
/// </para>
/// </summary>
public static class FractionalEffectiveDose
{
    /// <summary>
    /// Respiratory minute volume in litres per minute for an adult at light
    /// activity, as used in the standard correlation.
    /// </summary>
    public const double DefaultRespiratoryMinuteVolume = 25.0;

    /// <summary>
    /// Carboxyhaemoglobin concentration taken as incapacitating, in percent.
    /// </summary>
    public const double DefaultIncapacitatingCohb = 30.0;

    /// <summary>
    /// Mass extinction coefficient of smoke from flaming combustion, in m²/kg.
    /// </summary>
    public const double FlamingMassExtinctionCoefficient = 7600.0;

    /// <summary>
    /// Contrast factor for reflecting signs. Light-emitting signs use 8.
    /// </summary>
    public const double ReflectingSignContrast = 3.0;

    // ═══════════════════════════════════════════════════════════════
    //  Asphyxiant gases
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Fractional effective dose accumulated per minute of exposure to carbon
    /// monoxide:
    /// <para>rate = 3.317×10⁻⁵ · [CO]^1.036 · RMV / D</para>
    /// At the default respiratory rate and incapacitating COHb this is
    /// approximately [CO]^1.036 / 36000 per minute.
    /// </summary>
    /// <param name="coPpm">Carbon monoxide concentration in ppm.</param>
    /// <param name="respiratoryMinuteVolume">Breathing rate in L/min.</param>
    /// <param name="incapacitatingCohb">Incapacitating carboxyhaemoglobin level in percent.</param>
    public static double CarbonMonoxideRate(
        double coPpm,
        double respiratoryMinuteVolume = DefaultRespiratoryMinuteVolume,
        double incapacitatingCohb = DefaultIncapacitatingCohb)
    {
        if (coPpm <= 0) return 0;
        if (respiratoryMinuteVolume <= 0) throw new ArgumentException("Breathing rate must be positive.", nameof(respiratoryMinuteVolume));
        if (incapacitatingCohb <= 0) throw new ArgumentException("Incapacitating COHb must be positive.", nameof(incapacitatingCohb));

        return 3.317e-5 * Math.Pow(coPpm, 1.036) * respiratoryMinuteVolume / incapacitatingCohb;
    }

    /// <summary>
    /// Fractional effective dose accumulated per minute of exposure to hydrogen
    /// cyanide:
    /// <para>rate = exp(0.023 · [HCN] − 5.396)</para>
    /// <para>
    /// Cyanide acts far faster than carbon monoxide and is released in quantity
    /// by burning polyurethane foam and other nitrogen-bearing plastics, which
    /// is why it dominates in furnished spaces.
    /// </para>
    /// </summary>
    /// <param name="hcnPpm">Hydrogen cyanide concentration in ppm.</param>
    public static double HydrogenCyanideRate(double hcnPpm)
    {
        if (hcnPpm <= 0) return 0;
        return Math.Exp(0.023 * hcnPpm - 5.396);
    }

    /// <summary>
    /// Fractional effective dose accumulated per minute from oxygen depletion:
    /// <para>rate = 1 / exp(8.13 − 0.54 · (20.9 − [O₂]))</para>
    /// Returns 0 at or above the normal ambient 20.9 percent.
    /// </summary>
    /// <param name="oxygenPercent">Oxygen concentration in percent by volume.</param>
    public static double OxygenDepletionRate(double oxygenPercent)
    {
        if (oxygenPercent >= 20.9) return 0;
        return 1.0 / Math.Exp(8.13 - 0.54 * (20.9 - oxygenPercent));
    }

    /// <summary>
    /// Multiplier on asphyxiant uptake caused by the hyperventilation that carbon
    /// dioxide induces:
    /// <para>VCO₂ = exp(0.1903 · [CO₂] + 2.0004) / 7.1</para>
    /// <para>
    /// Carbon dioxide is not itself very toxic at fire concentrations, but it
    /// makes people breathe harder and so take in everything else faster.
    /// Returns 1 at ambient.
    /// </para>
    /// </summary>
    /// <param name="co2Percent">Carbon dioxide concentration in percent by volume.</param>
    public static double HyperventilationFactor(double co2Percent)
    {
        if (co2Percent <= 0) return 1.0;
        double factor = Math.Exp(0.1903 * co2Percent + 2.0004) / 7.1;
        return factor < 1.0 ? 1.0 : factor;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Combined dose
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Total fractional effective dose accumulated per minute, combining the
    /// asphyxiant gases and the effect of carbon dioxide on breathing rate:
    /// <para>rate = (rate_CO + rate_HCN) · VCO₂ + rate_O₂</para>
    /// A cumulative dose of 1.0 marks incapacitation.
    /// </summary>
    /// <param name="coPpm">Carbon monoxide concentration in ppm.</param>
    /// <param name="hcnPpm">Hydrogen cyanide concentration in ppm.</param>
    /// <param name="co2Percent">Carbon dioxide concentration in percent.</param>
    /// <param name="oxygenPercent">Oxygen concentration in percent.</param>
    public static double Rate(
        double coPpm,
        double hcnPpm = 0,
        double co2Percent = 0,
        double oxygenPercent = 20.9)
    {
        double asphyxiants = CarbonMonoxideRate(coPpm) + HydrogenCyanideRate(hcnPpm);
        return asphyxiants * HyperventilationFactor(co2Percent) + OxygenDepletionRate(oxygenPercent);
    }

    /// <summary>
    /// Time in minutes to reach a fractional effective dose of 1.0 while held at
    /// a constant exposure. Returns <see cref="double.PositiveInfinity"/> when the
    /// exposure never incapacitates.
    /// </summary>
    /// <param name="coPpm">Carbon monoxide concentration in ppm.</param>
    /// <param name="hcnPpm">Hydrogen cyanide concentration in ppm.</param>
    /// <param name="co2Percent">Carbon dioxide concentration in percent.</param>
    /// <param name="oxygenPercent">Oxygen concentration in percent.</param>
    public static double TimeToIncapacitation(
        double coPpm,
        double hcnPpm = 0,
        double co2Percent = 0,
        double oxygenPercent = 20.9)
    {
        double rate = Rate(coPpm, hcnPpm, co2Percent, oxygenPercent);
        return rate <= 0 ? double.PositiveInfinity : 1.0 / rate;
    }

    /// <summary>
    /// Accumulates a fractional effective dose over a time history of exposures.
    /// </summary>
    /// <param name="ratesPerMinute">Dose rate per minute at each step, from <see cref="Rate"/>.</param>
    /// <param name="stepSeconds">Length of each step in seconds.</param>
    public static double Accumulate(double[] ratesPerMinute, double stepSeconds)
    {
        if (ratesPerMinute == null) throw new ArgumentNullException(nameof(ratesPerMinute));
        if (stepSeconds <= 0) throw new ArgumentException("Step must be positive.", nameof(stepSeconds));

        double stepMinutes = stepSeconds / 60.0;
        double dose = 0;
        for (int i = 0; i < ratesPerMinute.Length; i++)
            dose += ratesPerMinute[i] * stepMinutes;
        return dose;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Visibility
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Visibility distance in metres through smoke of a given soot concentration:
    /// <para>S = C / (K_m · ρ_soot)</para>
    /// <para>
    /// Visibility is normally what stops escape first: people slow down, turn
    /// back or miss exits at distances where the toxic dose is still far below
    /// incapacitating. Ten metres is the usual tenability limit for a large
    /// space, five for a small one.
    /// </para>
    /// </summary>
    /// <param name="sootConcentrationKgPerM3">Soot mass concentration in kg/m³.</param>
    /// <param name="contrastFactor">3 for reflecting signs, 8 for light-emitting ones.</param>
    /// <param name="massExtinctionCoefficient">Mass extinction coefficient in m²/kg.</param>
    public static double Visibility(
        double sootConcentrationKgPerM3,
        double contrastFactor = ReflectingSignContrast,
        double massExtinctionCoefficient = FlamingMassExtinctionCoefficient)
    {
        if (sootConcentrationKgPerM3 <= 0) return double.PositiveInfinity;
        if (massExtinctionCoefficient <= 0) throw new ArgumentException("Extinction coefficient must be positive.", nameof(massExtinctionCoefficient));

        return contrastFactor / (massExtinctionCoefficient * sootConcentrationKgPerM3);
    }
}
