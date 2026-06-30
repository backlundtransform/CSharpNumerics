using System;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Astro.Enums;
using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Physics.OrbitalMechanics;

namespace CSharpNumerics.Physics.Astro;

/// <summary>
/// Approximate positions of the major planets from the low-precision Keplerian
/// element set published by E. M. Standish (JPL Solar System Dynamics). Valid
/// roughly 1800–2050 AD with an accuracy of a few arc-minutes.
///
/// The pipeline reuses the existing building blocks: Julian centuries from
/// <see cref="AstronomyExtensions"/>, Kepler's equation from
/// <see cref="KeplerOrbit"/>, and the perifocal→inertial rotation from
/// <see cref="ElementsToState"/>. Heliocentric coordinates are in the J2000
/// ecliptic frame; geocentric results are right ascension / declination in the
/// J2000 equatorial frame.
/// </summary>
public static class PlanetaryEphemeris
{
    private const double DegToRad = Math.PI / 180.0;
    private const double RadToDeg = 180.0 / Math.PI;

    // A non-zero gravitational parameter is required by ElementsToState; only the
    // returned position is used, which is independent of it.
    private static readonly double SunMu = PhysicsConstants.GravitationalConstant * PhysicsConstants.SolarMass;

    /// <summary>
    /// Keplerian elements and their per-century rates, J2000 epoch.
    /// Order: a (AU), e, I (deg), L mean longitude (deg), ϖ longitude of
    /// perihelion (deg), Ω longitude of ascending node (deg) — each followed by
    /// its rate (…/century). Source: Standish, "Keplerian Elements for
    /// Approximate Positions of the Major Planets" (1800–2050 AD set).
    /// </summary>
    private static readonly double[][] Elements =
    {
        // a, e, I, L, peri, node  /  rates
        new[] { 0.38709927, 0.20563593, 7.00497902, 252.25032350, 77.45779628, 48.33076593,
                0.00000037, 0.00001906, -0.00594749, 149472.67411175, 0.16047689, -0.12534081 }, // Mercury
        new[] { 0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255,
                0.00000390, -0.00004107, -0.00078890, 58517.81538729, 0.00268329, -0.27769418 }, // Venus
        new[] { 1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0,
                0.00000562, -0.00004392, -0.01294668, 35999.37244981, 0.32327364, 0.0 },          // Earth–Moon barycentre
        new[] { 1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891,
                0.00001847, 0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343 },   // Mars
        new[] { 5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983, 100.47390909,
                -0.00011607, -0.00013253, -0.00183714, 3034.74612775, 0.21252668, 0.20469106 },   // Jupiter
        new[] { 9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448,
                -0.00125060, -0.00050991, 0.00193609, 1222.49362201, -0.41897216, -0.28867794 },  // Saturn
        new[] { 19.18916464, 0.04725744, 0.77263783, 313.23810451, 170.95427630, 74.01692503,
                -0.00196176, -0.00004397, -0.00242939, 428.48202785, 0.40805281, 0.04240589 },    // Uranus
        new[] { 30.06992276, 0.00859048, 1.77004347, -55.12002969, 44.96476227, 131.78422574,
                0.00026291, 0.00005105, 0.00035372, 218.45945325, -0.32241464, -0.00508664 },     // Neptune
    };

    /// <summary>
    /// Heliocentric position of a planet in the J2000 ecliptic frame, in
    /// astronomical units.
    /// </summary>
    public static Vector HeliocentricEcliptic(Planet planet, DateTime utc)
    {
        double t = utc.JulianCenturiesSinceJ2000();
        double[] el = Elements[(int)planet];

        double a = el[0] + el[6] * t;            // AU
        double e = el[1] + el[7] * t;
        double inc = (el[2] + el[8] * t) * DegToRad;
        double meanLongitude = el[3] + el[9] * t;
        double longPerihelion = el[4] + el[10] * t;
        double longNode = el[5] + el[11] * t;

        double argPerihelion = (longPerihelion - longNode) * DegToRad;
        double meanAnomaly = NormalizeDegrees180(meanLongitude - longPerihelion) * DegToRad;

        double trueAnomaly = KeplerOrbit.TrueAnomaly(meanAnomaly, e);

        var elements = new OrbitalElements
        {
            SemiMajorAxis = a,
            Eccentricity = e,
            Inclination = inc,
            RAAN = longNode * DegToRad,
            ArgumentOfPeriapsis = argPerihelion,
            TrueAnomaly = trueAnomaly,
            Mu = SunMu
        };

        // Position is independent of Mu; velocity (unused here) is not.
        return ElementsToState.ToStateVector(elements).Position;
    }

    /// <summary>
    /// Geocentric apparent position of a planet: right ascension and declination
    /// (degrees, J2000 equatorial frame) together with the Earth–planet distance
    /// in astronomical units.
    /// </summary>
    /// <exception cref="ArgumentException">If <paramref name="planet"/> is Earth.</exception>
    public static (double RightAscensionDegrees, double DeclinationDegrees, double DistanceAU) GeocentricEquatorial(
        Planet planet, DateTime utc)
    {
        if (planet == Planet.Earth)
            throw new ArgumentException("Geocentric position is undefined for Earth.", nameof(planet));

        Vector planetHelio = HeliocentricEcliptic(planet, utc);
        Vector earthHelio = HeliocentricEcliptic(Planet.Earth, utc);

        // Geocentric ecliptic vector (AU).
        double gx = planetHelio.x - earthHelio.x;
        double gy = planetHelio.y - earthHelio.y;
        double gz = planetHelio.z - earthHelio.z;

        // Rotate ecliptic → equatorial about the x-axis by the obliquity.
        double t = utc.JulianCenturiesSinceJ2000();
        double obliquity = (23.43929 - 0.0130042 * t) * DegToRad;
        double cosE = Math.Cos(obliquity);
        double sinE = Math.Sin(obliquity);

        double xEq = gx;
        double yEq = gy * cosE - gz * sinE;
        double zEq = gy * sinE + gz * cosE;

        double distance = Math.Sqrt(gx * gx + gy * gy + gz * gz);
        double ra = NormalizeDegrees360(Math.Atan2(yEq, xEq) * RadToDeg);
        double dec = Math.Atan2(zEq, Math.Sqrt(xEq * xEq + yEq * yEq)) * RadToDeg;

        return (ra, dec, distance);
    }

    private static double NormalizeDegrees180(double deg)
    {
        deg %= 360.0;
        if (deg > 180.0) deg -= 360.0;
        if (deg < -180.0) deg += 360.0;
        return deg;
    }

    private static double NormalizeDegrees360(double deg)
    {
        deg %= 360.0;
        if (deg < 0) deg += 360.0;
        return deg;
    }
}
