using System;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Mechanics;
using CSharpNumerics.Physics.OrbitalMechanics;

namespace CSharpNumerics.Engines.Game.Rocket;

/// <summary>
/// Detects when a rocket's trajectory becomes a closed (bound) orbit.
/// An orbit is considered "inserted" when eccentricity &lt; 1 and periapsis is above the atmosphere.
/// </summary>
public class OrbitalInsertionDetector
{
    /// <summary>Minimum periapsis altitude (meters) to consider orbit valid. Default 100 km.</summary>
    public double MinPeriapsisAltitude { get; set; } = 100000.0;

    /// <summary>Maximum eccentricity to consider orbit "inserted". Default 1.0 (any bound orbit).</summary>
    public double MaxEccentricity { get; set; } = 1.0;

    /// <summary>Whether insertion has been detected.</summary>
    public bool IsInserted { get; private set; }

    /// <summary>Time at which insertion was detected (seconds).</summary>
    public double InsertionTime { get; private set; }

    /// <summary>Orbital elements at moment of insertion.</summary>
    public OrbitalElements InsertionElements { get; private set; }

    /// <summary>
    /// Checks the current state for orbital insertion conditions.
    /// </summary>
    /// <param name="position">Position in ECI (meters).</param>
    /// <param name="velocity">Velocity in ECI (m/s).</param>
    /// <param name="time">Current simulation time (seconds).</param>
    /// <returns>True if orbit insertion conditions are met.</returns>
    public bool Check(Vector position, Vector velocity, double time)
    {
        if (IsInserted) return true;

        var elements = StateToElements.FromStateVector(position, velocity, EarthModel.GM);

        // Check if orbit is bound (elliptical)
        if (elements.Eccentricity >= MaxEccentricity) return false;

        // Check periapsis is above atmosphere
        double periapsisAlt = elements.Periapsis - EarthModel.SemiMajorAxis;
        if (periapsisAlt < MinPeriapsisAltitude) return false;

        IsInserted = true;
        InsertionTime = time;
        InsertionElements = elements;
        return true;
    }

    /// <summary>Resets the detector state.</summary>
    public void Reset()
    {
        IsInserted = false;
        InsertionTime = 0;
        InsertionElements = default;
    }
}
