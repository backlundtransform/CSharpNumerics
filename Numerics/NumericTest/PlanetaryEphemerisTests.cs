using System;
using CSharpNumerics.Physics.Astro;
using CSharpNumerics.Physics.Astro.Enums;

namespace NumericsTests;

[TestClass]
public class PlanetaryEphemerisTests
{
    private static DateTime J2000 => new DateTime(2000, 1, 1, 12, 0, 0, DateTimeKind.Utc);

    // Perihelion / aphelion bounds (AU) for each planet — a(1∓e).
    private static readonly (Planet Planet, double Min, double Max)[] DistanceBounds =
    {
        (Planet.Mercury, 0.307, 0.467),
        (Planet.Venus,   0.718, 0.729),
        (Planet.Earth,   0.983, 1.017),
        (Planet.Mars,    1.381, 1.666),
        (Planet.Jupiter, 4.950, 5.455),
        (Planet.Saturn,  9.024, 10.054),
        (Planet.Uranus,  18.28, 20.10),
        (Planet.Neptune, 29.81, 30.33),
    };

    [TestMethod]
    public void HeliocentricDistances_AtJ2000_AreWithinPerihelionAndAphelion()
    {
        foreach (var (planet, min, max) in DistanceBounds)
        {
            double r = HeliocentricEclipticDistance(planet);
            Assert.IsTrue(r >= min && r <= max, $"{planet}: r = {r} AU outside [{min}, {max}]");
        }
    }

    [TestMethod]
    public void Earth_AtJ2000_IsNearPerihelion()
    {
        // Earth reaches perihelion (~0.9833 AU) in early January; at J2000 it is close.
        double r = HeliocentricEclipticDistance(Planet.Earth);
        Assert.AreEqual(0.9833, r, 0.001);
    }

    [TestMethod]
    public void Earth_HeliocentricLongitude_AtJ2000_IsAbout100Degrees()
    {
        // Corresponds to the Sun's geocentric longitude ≈ 280.5° at J2000.
        var p = PlanetaryEphemeris.HeliocentricEcliptic(Planet.Earth, J2000);
        double lon = Math.Atan2(p.y, p.x) * 180.0 / Math.PI;
        if (lon < 0) lon += 360.0;
        Assert.AreEqual(100.4, lon, 0.5);
    }

    [TestMethod]
    public void HeliocentricLatitude_StaysNearTheEcliptic()
    {
        // All major planets orbit within a few degrees of the ecliptic plane.
        foreach (var (planet, _, _) in DistanceBounds)
        {
            var p = PlanetaryEphemeris.HeliocentricEcliptic(planet, J2000);
            double r = Math.Sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
            double latDeg = Math.Asin(p.z / r) * 180.0 / Math.PI;
            Assert.IsTrue(Math.Abs(latDeg) < 8.0, $"{planet}: latitude {latDeg}° too far from ecliptic");
        }
    }

    [TestMethod]
    public void GeocentricEquatorial_Mars_ReturnsValidCoordinates()
    {
        var (ra, dec, dist) = PlanetaryEphemeris.GeocentricEquatorial(Planet.Mars, J2000);

        Assert.IsTrue(ra >= 0.0 && ra < 360.0, $"RA {ra} out of range");
        Assert.IsTrue(dec >= -90.0 && dec <= 90.0, $"Dec {dec} out of range");
        Assert.IsTrue(dist > 0.37 && dist < 2.70, $"Mars distance {dist} AU out of range"); // min/max Earth–Mars
    }

    [TestMethod]
    public void GeocentricEquatorial_IsDeterministic()
    {
        var a = PlanetaryEphemeris.GeocentricEquatorial(Planet.Jupiter, J2000);
        var b = PlanetaryEphemeris.GeocentricEquatorial(Planet.Jupiter, J2000);
        Assert.AreEqual(a.RightAscensionDegrees, b.RightAscensionDegrees, 1e-12);
        Assert.AreEqual(a.DeclinationDegrees, b.DeclinationDegrees, 1e-12);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void GeocentricEquatorial_Earth_Throws()
    {
        PlanetaryEphemeris.GeocentricEquatorial(Planet.Earth, J2000);
    }

    private static double HeliocentricEclipticDistance(Planet planet)
    {
        var p = PlanetaryEphemeris.HeliocentricEcliptic(planet, J2000);
        return Math.Sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
    }
}
