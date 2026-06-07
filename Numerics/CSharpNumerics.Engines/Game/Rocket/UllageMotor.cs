using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Engines.Game.Rocket;

/// <summary>
/// Small solid motor fired after stage separation to provide
/// a brief acceleration that settles propellant to the tank bottom
/// before main engine ignition (ullage settling).
/// 
/// Prevents vapor ingestion and ensures proper feed to turbopumps.
/// </summary>
public class UllageMotor
{
    /// <summary>Thrust of the ullage motor in Newtons.</summary>
    public double Thrust { get; }

    /// <summary>Burn duration of the ullage motor in seconds.</summary>
    public double BurnDuration { get; }

    /// <summary>Mass of the ullage motor (dry + propellant) in kg.</summary>
    public double Mass { get; }

    /// <summary>Number of ullage motors fired simultaneously.</summary>
    public int Count { get; }

    /// <summary>Time at which ullage motor was activated (set when Fire is called).</summary>
    public double IgnitionTime { get; private set; } = -1;

    /// <summary>True if the ullage motor has been fired.</summary>
    public bool HasFired => IgnitionTime >= 0;

    /// <summary>
    /// Creates an ullage motor specification.
    /// </summary>
    /// <param name="thrust">Thrust per motor in Newtons (typical 1–10 kN).</param>
    /// <param name="burnDuration">Burn time in seconds (typical 1–4s).</param>
    /// <param name="mass">Mass per motor in kg.</param>
    /// <param name="count">Number of motors (typical 2–4).</param>
    public UllageMotor(double thrust, double burnDuration, double mass, int count = 2)
    {
        if (thrust <= 0) throw new ArgumentOutOfRangeException(nameof(thrust));
        if (burnDuration <= 0) throw new ArgumentOutOfRangeException(nameof(burnDuration));

        Thrust = thrust;
        BurnDuration = burnDuration;
        Mass = mass;
        Count = count;
    }

    /// <summary>Total thrust from all ullage motors (N).</summary>
    public double TotalThrust => Thrust * Count;

    /// <summary>Total mass of all ullage motors (kg).</summary>
    public double TotalMass => Mass * Count;

    /// <summary>Total impulse delivered by all motors (N·s).</summary>
    public double TotalImpulse => TotalThrust * BurnDuration;

    /// <summary>
    /// Fires the ullage motor at the given simulation time.
    /// </summary>
    /// <param name="time">Current simulation time in seconds.</param>
    public void Fire(double time)
    {
        IgnitionTime = time;
    }

    /// <summary>
    /// Returns current thrust at the given time. Zero if not fired or burn complete.
    /// </summary>
    /// <param name="currentTime">Current simulation time in seconds.</param>
    public double CurrentThrust(double currentTime)
    {
        if (!HasFired) return 0;
        double elapsed = currentTime - IgnitionTime;
        if (elapsed < 0 || elapsed > BurnDuration) return 0;
        return TotalThrust;
    }

    /// <summary>
    /// True if the ullage burn is currently active.
    /// </summary>
    /// <param name="currentTime">Current simulation time in seconds.</param>
    public bool IsBurning(double currentTime)
    {
        if (!HasFired) return 0 > 0;
        double elapsed = currentTime - IgnitionTime;
        return elapsed >= 0 && elapsed <= BurnDuration;
    }

    /// <summary>
    /// Resets the motor to unfired state.
    /// </summary>
    public void Reset()
    {
        IgnitionTime = -1;
    }
}
