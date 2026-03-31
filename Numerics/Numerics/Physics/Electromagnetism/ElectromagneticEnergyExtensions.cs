using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Electromagnetism;

/// <summary>
/// Extension methods for electromagnetic energy and momentum:
/// energy density, Poynting vector, and radiation pressure.
/// </summary>
public static class ElectromagneticEnergyExtensions
{
    // ═══════════════════════════════════════════════════════════════
    //  Energy & momentum
    // ═══════════════════════════════════════════════════════════════

    #region Energy & Momentum

    /// <summary>
    /// Computes the electromagnetic energy density:
    /// u = ½(ε₀|E|² + |B|²/μ₀).
    /// </summary>
    /// <param name="electricField">Electric field vector at the point.</param>
    /// <param name="magneticField">Magnetic field vector at the point.</param>
    public static double EnergyDensity(this Vector electricField, Vector magneticField)
    {
        var eMag2 = electricField.Dot(electricField);
        var bMag2 = magneticField.Dot(magneticField);
        return 0.5 * (PhysicsConstants.VacuumPermittivity * eMag2
                    + bMag2 / PhysicsConstants.VacuumPermeability);
    }

    /// <summary>
    /// Computes the Poynting vector: S = (1/μ₀)(E × B).
    /// Represents the directional energy flux (power per unit area).
    /// </summary>
    /// <param name="electricField">Electric field vector at the point.</param>
    /// <param name="magneticField">Magnetic field vector at the point.</param>
    public static Vector PoyntingVector(this Vector electricField, Vector magneticField)
    {
        return (1.0 / PhysicsConstants.VacuumPermeability) * electricField.Cross(magneticField);
    }

    /// <summary>
    /// Computes the radiation pressure from the Poynting vector magnitude:
    /// P = |S| / c (for full absorption) or P = 2|S| / c (for full reflection).
    /// </summary>
    /// <param name="electricField">Electric field vector.</param>
    /// <param name="magneticField">Magnetic field vector.</param>
    /// <param name="reflected">If true, assumes total reflection (pressure doubled).</param>
    public static double RadiationPressure(this Vector electricField, Vector magneticField, bool reflected = false)
    {
        var sMag = electricField.PoyntingVector(magneticField).GetMagnitude();
        return (reflected ? 2.0 : 1.0) * sMag / PhysicsConstants.SpeedOfLight;
    }

    #endregion
}
