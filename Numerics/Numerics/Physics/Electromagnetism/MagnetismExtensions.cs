using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics;
using System;

namespace CSharpNumerics.Physics.Electromagnetism;

/// <summary>
/// Extension methods for magnetism: Lorentz force and Biot–Savart magnetic field sources
/// (infinite wire, circular loop, solenoid).
/// </summary>
public static class MagnetismExtensions
{
    // ═══════════════════════════════════════════════════════════════
    //  Lorentz force
    // ═══════════════════════════════════════════════════════════════

    #region Lorentz Force

    /// <summary>
    /// Computes the Lorentz force on a charged particle: F = q(E + v × B).
    /// </summary>
    /// <param name="charge">Charge in coulombs.</param>
    /// <param name="velocity">Velocity of the particle in m/s.</param>
    /// <param name="electricField">Electric field at the particle position.</param>
    /// <param name="magneticField">Magnetic field at the particle position.</param>
    public static Vector LorentzForce(this double charge, Vector velocity, Vector electricField, Vector magneticField)
    {
        return charge * (electricField + velocity.Cross(magneticField));
    }

    /// <summary>
    /// Computes only the electric part of the Lorentz force: F = qE.
    /// </summary>
    /// <param name="charge">Charge in coulombs.</param>
    /// <param name="electricField">Electric field at the particle position.</param>
    public static Vector ElectricForce(this double charge, Vector electricField)
    {
        return charge * electricField;
    }

    /// <summary>
    /// Computes only the magnetic part of the Lorentz force: F = q(v × B).
    /// </summary>
    /// <param name="charge">Charge in coulombs.</param>
    /// <param name="velocity">Velocity of the particle in m/s.</param>
    /// <param name="magneticField">Magnetic field at the particle position.</param>
    public static Vector MagneticForce(this double charge, Vector velocity, Vector magneticField)
    {
        return charge * velocity.Cross(magneticField);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Magnetic field sources (Biot–Savart)
    // ═══════════════════════════════════════════════════════════════

    #region Biot–Savart

    /// <summary>
    /// Computes the magnetic field magnitude at perpendicular distance d
    /// from an infinite straight wire: B = μ₀I / (2πd).
    /// </summary>
    /// <param name="current">Current in amperes.</param>
    /// <param name="distance">Perpendicular distance from the wire in meters.</param>
    public static double MagneticFieldFromWire(this double current, double distance)
    {
        if (distance <= 0) throw new ArgumentException("Distance must be greater than zero.");
        return PhysicsConstants.VacuumPermeability * Math.Abs(current) / (2.0 * Math.PI * distance);
    }

    /// <summary>
    /// Computes the magnetic field vector from an infinite straight wire
    /// at a given field point: B = (μ₀I / 2πd) · (ŵ × ρ̂),
    /// where ŵ is the wire direction and ρ̂ is the perpendicular direction
    /// from the wire to the field point.
    /// </summary>
    /// <param name="current">Current in amperes (sign determines direction along wire).</param>
    /// <param name="wireDirection">Direction of current flow (need not be unit length).</param>
    /// <param name="wirePoint">Any point on the wire.</param>
    /// <param name="fieldPoint">Position where the field is evaluated.</param>
    public static Vector MagneticFieldFromWire(this double current, Vector wireDirection, Vector wirePoint, Vector fieldPoint)
    {
        var wireUnit = wireDirection.GetUnitVector();
        var r = fieldPoint - wirePoint;
        var rParallel = r.Dot(wireUnit) * wireUnit;
        var rPerp = r - rParallel;
        var d = rPerp.GetMagnitude();
        if (d == 0) throw new ArgumentException("Field point lies on the wire.");

        var bDirection = wireUnit.Cross(rPerp.GetUnitVector());
        var bMagnitude = PhysicsConstants.VacuumPermeability * current / (2.0 * Math.PI * d);

        return bMagnitude * bDirection;
    }

    /// <summary>
    /// Computes the magnetic field at the center of a circular current loop:
    /// B = μ₀I / (2R), directed along the loop normal.
    /// </summary>
    /// <param name="current">Current in amperes.</param>
    /// <param name="radius">Radius of the loop in meters.</param>
    /// <param name="normal">Normal direction of the loop (right-hand rule).</param>
    public static Vector MagneticFieldFromLoop(this double current, double radius, Vector normal)
    {
        if (radius <= 0) throw new ArgumentException("Radius must be greater than zero.");
        var bMagnitude = PhysicsConstants.VacuumPermeability * current / (2.0 * radius);
        return bMagnitude * normal.GetUnitVector();
    }

    /// <summary>
    /// Computes the magnetic field inside an ideal solenoid: B = μ₀nI,
    /// where n = turnsPerUnitLength.
    /// </summary>
    /// <param name="current">Current in amperes.</param>
    /// <param name="turnsPerUnitLength">Number of turns per meter.</param>
    /// <param name="direction">Direction along the solenoid axis.</param>
    public static Vector MagneticFieldFromSolenoid(this double current, double turnsPerUnitLength, Vector direction)
    {
        var bMagnitude = PhysicsConstants.VacuumPermeability * turnsPerUnitLength * current;
        return bMagnitude * direction.GetUnitVector();
    }

    #endregion
}
