using System;

namespace CSharpNumerics.Physics.Environmental.Water;

/// <summary>
/// Manning's equation for open-channel flow — computes cross-sectional
/// mean velocity and discharge from channel geometry and bed slope.
/// <para>
/// All methods are pure functions. No grid or simulation state.
/// </para>
/// </summary>
public static class ManningEquation
{
    // ═══════════════════════════════════════════════════════════════
    //  Primary equations
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Cross-sectional mean velocity (m/s) via Manning's equation:
    /// <para>u = (1/n) · Rh^(2/3) · S₀^(1/2)</para>
    /// Returns 0 if any input is non-positive.
    /// </summary>
    /// <param name="manningN">Manning's roughness coefficient n (dimensionless).
    /// Typical values: 0.020 (clean concrete), 0.035 (natural gravel), 0.060 (weedy).</param>
    /// <param name="hydraulicRadius">Hydraulic radius Rh (m) = A / P.</param>
    /// <param name="bedSlope">Longitudinal bed slope S₀ (m/m). Must be > 0.</param>
    public static double Velocity(double manningN, double hydraulicRadius, double bedSlope)
    {
        if (manningN <= 0 || hydraulicRadius <= 0 || bedSlope <= 0)
            return 0;

        return (1.0 / manningN) * Math.Pow(hydraulicRadius, 2.0 / 3.0) * Math.Sqrt(bedSlope);
    }

    /// <summary>
    /// Volumetric discharge Q (m³/s) = u · A.
    /// <para>Q = (A/n) · Rh^(2/3) · S₀^(1/2)</para>
    /// </summary>
    /// <param name="manningN">Manning's roughness coefficient n.</param>
    /// <param name="hydraulicRadius">Hydraulic radius Rh (m).</param>
    /// <param name="bedSlope">Bed slope S₀ (m/m).</param>
    /// <param name="crossSectionalArea">Flow cross-sectional area A (m²).</param>
    public static double Discharge(double manningN, double hydraulicRadius, double bedSlope, double crossSectionalArea)
    {
        return Velocity(manningN, hydraulicRadius, bedSlope) * crossSectionalArea;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Hydraulic radius helpers
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Generic hydraulic radius Rh = A / P (m).
    /// </summary>
    /// <param name="crossSectionalArea">Flow area A (m²).</param>
    /// <param name="wettedPerimeter">Wetted perimeter P (m).</param>
    public static double HydraulicRadius(double crossSectionalArea, double wettedPerimeter)
    {
        if (wettedPerimeter <= 0) return 0;
        return crossSectionalArea / wettedPerimeter;
    }

    /// <summary>
    /// Hydraulic radius for a rectangular channel: Rh = (W · d) / (W + 2d).
    /// </summary>
    /// <param name="width">Channel width W (m).</param>
    /// <param name="depth">Flow depth d (m).</param>
    public static double RectangularHydraulicRadius(double width, double depth)
    {
        if (width <= 0 || depth <= 0) return 0;
        double area = width * depth;
        double perimeter = width + 2.0 * depth;
        return area / perimeter;
    }

    /// <summary>
    /// Hydraulic radius for a trapezoidal channel:
    /// A = (b + z·d) · d,  P = b + 2d · √(1 + z²),  Rh = A / P.
    /// </summary>
    /// <param name="bottomWidth">Channel bottom width b (m).</param>
    /// <param name="depth">Flow depth d (m).</param>
    /// <param name="sideSlope">Side slope z (horizontal:vertical, e.g. 2 means 2H:1V).</param>
    public static double TrapezoidalHydraulicRadius(double bottomWidth, double depth, double sideSlope)
    {
        if (bottomWidth <= 0 || depth <= 0 || sideSlope < 0) return 0;
        double area = (bottomWidth + sideSlope * depth) * depth;
        double perimeter = bottomWidth + 2.0 * depth * Math.Sqrt(1.0 + sideSlope * sideSlope);
        if (perimeter <= 0) return 0;
        return area / perimeter;
    }
}
