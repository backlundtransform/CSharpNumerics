using System;

namespace CSharpNumerics.Physics.Mechanics.Propulsion;

/// <summary>
/// Tsiolkovsky rocket equation utilities.
/// Provides ΔV, burn time, and mass ratio calculations for ideal rocket propulsion.
/// </summary>
public static class TsiolkovskyExtensions
{
    private const double G0 = 9.80665; // m/s²

    /// <summary>
    /// Computes ideal ΔV from the Tsiolkovsky rocket equation:
    /// ΔV = Isp · g₀ · ln(m₀ / mf)
    /// </summary>
    /// <param name="isp">Specific impulse in seconds.</param>
    /// <param name="initialMass">Initial total mass (wet) in kg.</param>
    /// <param name="finalMass">Final mass (dry + remaining propellant) in kg.</param>
    /// <returns>Ideal delta-V in m/s.</returns>
    public static double DeltaV(double isp, double initialMass, double finalMass)
    {
        if (isp <= 0) throw new ArgumentOutOfRangeException(nameof(isp));
        if (initialMass <= 0) throw new ArgumentOutOfRangeException(nameof(initialMass));
        if (finalMass <= 0) throw new ArgumentOutOfRangeException(nameof(finalMass));
        if (finalMass > initialMass) throw new ArgumentException("Final mass cannot exceed initial mass.");

        return isp * G0 * Math.Log(initialMass / finalMass);
    }

    /// <summary>
    /// Computes the mass ratio required for a given ΔV:
    /// m₀/mf = exp(ΔV / (Isp · g₀))
    /// </summary>
    /// <param name="deltaV">Required delta-V in m/s.</param>
    /// <param name="isp">Specific impulse in seconds.</param>
    /// <returns>Mass ratio m₀/mf (always ≥ 1).</returns>
    public static double MassRatio(double deltaV, double isp)
    {
        if (isp <= 0) throw new ArgumentOutOfRangeException(nameof(isp));
        if (deltaV < 0) throw new ArgumentOutOfRangeException(nameof(deltaV));

        return Math.Exp(deltaV / (isp * G0));
    }

    /// <summary>
    /// Computes the propellant mass needed for a given ΔV.
    /// m_prop = m_dry · (MassRatio - 1)
    /// </summary>
    /// <param name="deltaV">Required delta-V in m/s.</param>
    /// <param name="isp">Specific impulse in seconds.</param>
    /// <param name="dryMass">Dry mass (payload + structure) in kg.</param>
    /// <returns>Required propellant mass in kg.</returns>
    public static double PropellantMass(double deltaV, double isp, double dryMass)
    {
        if (dryMass <= 0) throw new ArgumentOutOfRangeException(nameof(dryMass));
        double ratio = MassRatio(deltaV, isp);
        return dryMass * (ratio - 1.0);
    }

    /// <summary>
    /// Computes burn time for constant-thrust propulsion:
    /// t_burn = m_prop / ṁ = m_prop · Isp · g₀ / T
    /// </summary>
    /// <param name="propellantMass">Propellant mass to burn in kg.</param>
    /// <param name="thrust">Constant thrust in Newtons.</param>
    /// <param name="isp">Specific impulse in seconds.</param>
    /// <returns>Burn duration in seconds.</returns>
    public static double BurnTime(double propellantMass, double thrust, double isp)
    {
        if (propellantMass <= 0) throw new ArgumentOutOfRangeException(nameof(propellantMass));
        if (thrust <= 0) throw new ArgumentOutOfRangeException(nameof(thrust));
        if (isp <= 0) throw new ArgumentOutOfRangeException(nameof(isp));

        double massFlowRate = thrust / (isp * G0);
        return propellantMass / massFlowRate;
    }

    /// <summary>
    /// Computes exhaust velocity from specific impulse:
    /// Ve = Isp · g₀
    /// </summary>
    /// <param name="isp">Specific impulse in seconds.</param>
    /// <returns>Effective exhaust velocity in m/s.</returns>
    public static double ExhaustVelocity(double isp)
    {
        if (isp <= 0) throw new ArgumentOutOfRangeException(nameof(isp));
        return isp * G0;
    }

    /// <summary>
    /// Computes payload fraction for a single stage:
    /// λ = m_payload / m₀ = (1/MR) - ε / (1 - ε)
    /// where ε = structural mass fraction = m_struct / (m_struct + m_prop)
    /// </summary>
    /// <param name="deltaV">Required delta-V in m/s.</param>
    /// <param name="isp">Specific impulse in seconds.</param>
    /// <param name="structuralFraction">Structural fraction ε (typically 0.05–0.15).</param>
    /// <returns>Payload fraction (may be negative if ΔV is unachievable).</returns>
    public static double PayloadFraction(double deltaV, double isp, double structuralFraction)
    {
        if (structuralFraction <= 0 || structuralFraction >= 1)
            throw new ArgumentOutOfRangeException(nameof(structuralFraction));

        double mr = MassRatio(deltaV, isp);
        return (1.0 / mr - structuralFraction) / (1.0 - structuralFraction);
    }
}
