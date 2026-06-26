using System;
using CSharpNumerics.Physics.Constants;

namespace CSharpNumerics.Physics.Relativity;

/// <summary>
/// Special-relativistic kinematics and dynamics in flat spacetime.
/// All speeds are in m/s and must be strictly below the speed of light.
/// These quantities are the foundation for the general-relativistic models
/// (Schwarzschild geometry, GPS-style clock corrections, etc.).
/// </summary>
public static class SpecialRelativity
{
    private const double C = PhysicsConstants.SpeedOfLight;
    private const double C2 = C * C;

    /// <summary>
    /// Lorentz factor γ = 1 / √(1 − v²/c²).
    /// </summary>
    /// <param name="velocity">Speed in m/s (sign is irrelevant).</param>
    /// <returns>γ ≥ 1.</returns>
    /// <exception cref="ArgumentOutOfRangeException">If |v| ≥ c.</exception>
    public static double LorentzFactor(double velocity)
    {
        double beta2 = velocity * velocity / C2;
        if (beta2 >= 1.0)
            throw new ArgumentOutOfRangeException(nameof(velocity), "Speed must be below the speed of light.");

        return 1.0 / Math.Sqrt(1.0 - beta2);
    }

    /// <summary>
    /// Time dilation: the coordinate time Δt elapsed in the rest frame while a
    /// clock moving at <paramref name="velocity"/> measures proper time
    /// <paramref name="properTime"/>. Δt = γ·Δτ (the moving clock runs slow).
    /// </summary>
    public static double TimeDilation(double properTime, double velocity)
    {
        return LorentzFactor(velocity) * properTime;
    }

    /// <summary>
    /// Length contraction: the length L = L₀/γ measured for an object whose
    /// proper (rest-frame) length is <paramref name="properLength"/> when it
    /// moves at <paramref name="velocity"/> along the measured direction.
    /// </summary>
    public static double LengthContraction(double properLength, double velocity)
    {
        return properLength / LorentzFactor(velocity);
    }

    /// <summary>
    /// Relativistic momentum p = γ·m·v.
    /// </summary>
    public static double Momentum(double mass, double velocity)
    {
        return LorentzFactor(velocity) * mass * velocity;
    }

    /// <summary>
    /// Rest energy E₀ = m·c².
    /// </summary>
    public static double RestEnergy(double mass)
    {
        return mass * C2;
    }

    /// <summary>
    /// Total relativistic energy E = γ·m·c².
    /// </summary>
    public static double TotalEnergy(double mass, double velocity)
    {
        return LorentzFactor(velocity) * mass * C2;
    }

    /// <summary>
    /// Relativistic kinetic energy E_k = (γ − 1)·m·c².
    /// Evaluated as m·c²·β²/(s·(1 + s)) with s = √(1 − β²), which is
    /// algebraically identical to (γ − 1)·m·c² but avoids the catastrophic
    /// cancellation of computing γ − 1 directly at low speeds (where it
    /// correctly reduces to the classical ½·m·v²).
    /// </summary>
    public static double KineticEnergy(double mass, double velocity)
    {
        double beta2 = velocity * velocity / C2;
        if (beta2 >= 1.0)
            throw new ArgumentOutOfRangeException(nameof(velocity), "Speed must be below the speed of light.");

        double s = Math.Sqrt(1.0 - beta2);
        return mass * C2 * beta2 / (s * (1.0 + s));
    }

    /// <summary>
    /// Total energy from the energy–momentum relation E = √((p·c)² + (m·c²)²).
    /// </summary>
    public static double EnergyFromMomentum(double mass, double momentum)
    {
        double pc = momentum * C;
        double mc2 = mass * C2;
        return Math.Sqrt(pc * pc + mc2 * mc2);
    }

    /// <summary>
    /// Relativistic velocity addition: combines two collinear velocities
    /// s = (u + v) / (1 + u·v/c²). The result is always below c when both
    /// inputs are.
    /// </summary>
    public static double AddVelocities(double u, double v)
    {
        return (u + v) / (1.0 + u * v / C2);
    }

    /// <summary>
    /// Relativistic longitudinal Doppler factor f_observed / f_source for purely
    /// radial motion: √((1 − β)/(1 + β)), with β = v/c.
    /// </summary>
    /// <param name="velocity">Radial speed in m/s. Positive = receding (redshift),
    /// negative = approaching (blueshift).</param>
    public static double DopplerFactor(double velocity)
    {
        double beta = velocity / C;
        if (Math.Abs(beta) >= 1.0)
            throw new ArgumentOutOfRangeException(nameof(velocity), "Speed must be below the speed of light.");

        return Math.Sqrt((1.0 - beta) / (1.0 + beta));
    }

    /// <summary>
    /// Observed frequency under the relativistic longitudinal Doppler effect.
    /// </summary>
    /// <param name="sourceFrequency">Emitted frequency (Hz).</param>
    /// <param name="velocity">Radial speed in m/s. Positive = receding.</param>
    public static double ObservedFrequency(double sourceFrequency, double velocity)
    {
        return sourceFrequency * DopplerFactor(velocity);
    }

    /// <summary>
    /// Rapidity φ = artanh(v/c), the additive measure of velocity
    /// (rapidities add linearly under collinear boosts).
    /// </summary>
    public static double Rapidity(double velocity)
    {
        double beta = velocity / C;
        if (Math.Abs(beta) >= 1.0)
            throw new ArgumentOutOfRangeException(nameof(velocity), "Speed must be below the speed of light.");

        return 0.5 * Math.Log((1.0 + beta) / (1.0 - beta));
    }
}
