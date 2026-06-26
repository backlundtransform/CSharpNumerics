using System;
using CSharpNumerics.Physics.Constants;

namespace CSharpNumerics.Physics.Relativity;

/// <summary>
/// Weak-field general-relativistic effects — the classical tests of GR:
/// perihelion precession, gravitational light deflection, and the Shapiro
/// time delay. Masses are in kg, lengths in metres, angles in radians.
/// </summary>
public static class GeneralRelativity
{
    private const double G = PhysicsConstants.GravitationalConstant;
    private const double C = PhysicsConstants.SpeedOfLight;
    private const double C2 = C * C;
    private const double C3 = C * C * C;

    private const double RadiansToArcseconds = 180.0 / Math.PI * 3600.0;

    /// <summary>
    /// Relativistic perihelion precession per orbit (radians):
    /// Δϖ = 6πGM / (c²·a·(1 − e²)).
    /// For Mercury this is ≈ 5.0×10⁻⁷ rad per orbit (≈ 43″ per century).
    /// </summary>
    /// <param name="semiMajorAxis">Orbital semi-major axis a (m).</param>
    /// <param name="eccentricity">Orbital eccentricity e, 0 ≤ e &lt; 1.</param>
    /// <param name="centralMass">Mass of the central body (kg).</param>
    public static double PerihelionPrecessionPerOrbit(double semiMajorAxis, double eccentricity, double centralMass)
    {
        if (semiMajorAxis <= 0)
            throw new ArgumentOutOfRangeException(nameof(semiMajorAxis), "Semi-major axis must be positive.");
        if (eccentricity < 0 || eccentricity >= 1)
            throw new ArgumentOutOfRangeException(nameof(eccentricity), "Eccentricity must be in [0, 1).");

        double latusFactor = semiMajorAxis * (1.0 - eccentricity * eccentricity);
        return 6.0 * Math.PI * G * centralMass / (C2 * latusFactor);
    }

    /// <summary>
    /// Relativistic perihelion precession in arcseconds per century, given the
    /// orbital period. Convenience wrapper around
    /// <see cref="PerihelionPrecessionPerOrbit"/>.
    /// </summary>
    /// <param name="semiMajorAxis">Semi-major axis a (m).</param>
    /// <param name="eccentricity">Eccentricity e, 0 ≤ e &lt; 1.</param>
    /// <param name="centralMass">Central mass (kg).</param>
    /// <param name="orbitalPeriodDays">Orbital period (days).</param>
    public static double PerihelionPrecessionArcsecPerCentury(
        double semiMajorAxis, double eccentricity, double centralMass, double orbitalPeriodDays)
    {
        if (orbitalPeriodDays <= 0)
            throw new ArgumentOutOfRangeException(nameof(orbitalPeriodDays), "Period must be positive.");

        double perOrbit = PerihelionPrecessionPerOrbit(semiMajorAxis, eccentricity, centralMass);
        double orbitsPerCentury = 36525.0 / orbitalPeriodDays;
        return perOrbit * orbitsPerCentury * RadiansToArcseconds;
    }

    /// <summary>
    /// Gravitational deflection of light passing a mass at impact parameter
    /// <paramref name="impactParameter"/>: α = 4GM / (c²·b) radians.
    /// For a ray grazing the Sun's limb this is ≈ 1.75″.
    /// </summary>
    /// <param name="impactParameter">Closest approach distance b (m).</param>
    /// <param name="centralMass">Deflecting mass (kg).</param>
    public static double LightDeflection(double impactParameter, double centralMass)
    {
        if (impactParameter <= 0)
            throw new ArgumentOutOfRangeException(nameof(impactParameter), "Impact parameter must be positive.");

        return 4.0 * G * centralMass / (C2 * impactParameter);
    }

    /// <summary>
    /// Shapiro time delay (seconds) — the extra one-way light-travel time for a
    /// signal passing a mass <paramref name="centralMass"/> at impact parameter
    /// <paramref name="impactParameter"/>, between a source at radial distance
    /// <paramref name="r1"/> and a receiver at <paramref name="r2"/>:
    /// Δt = (2GM/c³)·ln[ (r₁+√(r₁²−b²))(r₂+√(r₂²−b²)) / b² ].
    /// </summary>
    /// <param name="r1">Radial distance of the source from the mass (m).</param>
    /// <param name="r2">Radial distance of the receiver from the mass (m).</param>
    /// <param name="impactParameter">Impact parameter b of the signal path (m).</param>
    /// <param name="centralMass">Mass causing the delay (kg).</param>
    public static double ShapiroDelay(double r1, double r2, double impactParameter, double centralMass)
    {
        if (impactParameter <= 0)
            throw new ArgumentOutOfRangeException(nameof(impactParameter), "Impact parameter must be positive.");
        if (r1 < impactParameter || r2 < impactParameter)
            throw new ArgumentOutOfRangeException(nameof(r1),
                "Both radial distances must be at least the impact parameter.");

        double x1 = Math.Sqrt(r1 * r1 - impactParameter * impactParameter);
        double x2 = Math.Sqrt(r2 * r2 - impactParameter * impactParameter);
        double argument = (r1 + x1) * (r2 + x2) / (impactParameter * impactParameter);

        return 2.0 * G * centralMass / C3 * Math.Log(argument);
    }
}
