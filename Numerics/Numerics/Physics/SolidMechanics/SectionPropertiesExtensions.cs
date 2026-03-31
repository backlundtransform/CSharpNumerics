using System;

namespace CSharpNumerics.Physics.SolidMechanics;

/// <summary>
/// Extension methods for cross-section properties:
/// second moment of area (area moment of inertia) for rectangular,
/// circular, and tubular sections.
/// </summary>
public static class SectionPropertiesExtensions
{
    // ═══════════════════════════════════════════════════════════════
    //  Second Moment of Area (Area Moment of Inertia)
    // ═══════════════════════════════════════════════════════════════

    #region Second Moment of Area

    /// <summary>
    /// Second moment of area for a solid rectangle about the centroidal axis: I = bh³ / 12.
    /// </summary>
    /// <param name="width">Width b in metres.</param>
    /// <param name="height">Height h in metres.</param>
    /// <returns>Second moment of area I in m⁴.</returns>
    public static double RectangularSecondMoment(this double width, double height)
    {
        if (width <= 0) throw new ArgumentException("Width must be greater than zero.");
        if (height <= 0) throw new ArgumentException("Height must be greater than zero.");
        return width * height * height * height / 12.0;
    }

    /// <summary>
    /// Second moment of area for a solid circle about a diametral axis: I = πr⁴ / 4.
    /// </summary>
    /// <param name="radius">Radius r in metres.</param>
    /// <returns>Second moment of area I in m⁴.</returns>
    public static double CircularSecondMoment(this double radius)
    {
        if (radius <= 0) throw new ArgumentException("Radius must be greater than zero.");
        return Math.PI * radius * radius * radius * radius / 4.0;
    }

    /// <summary>
    /// Second moment of area for a hollow tube: I = π(R⁴ − r⁴) / 4.
    /// </summary>
    /// <param name="outerRadius">Outer radius R in metres.</param>
    /// <param name="innerRadius">Inner radius r in metres.</param>
    /// <returns>Second moment of area I in m⁴.</returns>
    public static double TubularSecondMoment(this double outerRadius, double innerRadius)
    {
        if (outerRadius <= 0) throw new ArgumentException("Outer radius must be greater than zero.");
        if (innerRadius < 0 || innerRadius >= outerRadius)
            throw new ArgumentException("Inner radius must be non-negative and less than outer radius.");
        double R4 = outerRadius * outerRadius * outerRadius * outerRadius;
        double r4 = innerRadius * innerRadius * innerRadius * innerRadius;
        return Math.PI * (R4 - r4) / 4.0;
    }

    #endregion
}
