using System;
using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Physics.Relativity;

namespace NumericsTests;

[TestClass]
public class GeneralRelativityTests
{
    private const double SolarMass = PhysicsConstants.SolarMass;
    private const double RadToArcsec = 180.0 / Math.PI * 3600.0;

    // Mercury's orbit
    private const double MercurySemiMajorAxis = 5.79e10;   // m
    private const double MercuryEccentricity = 0.2056;
    private const double MercuryPeriodDays = 87.969;

    [TestMethod]
    public void PerihelionPrecession_Mercury_PerOrbit()
    {
        // ≈ 5.0×10⁻⁷ rad per orbit
        double perOrbit = GeneralRelativity.PerihelionPrecessionPerOrbit(
            MercurySemiMajorAxis, MercuryEccentricity, SolarMass);
        Assert.AreEqual(5.02e-7, perOrbit, 0.1e-7);
    }

    [TestMethod]
    public void PerihelionPrecession_Mercury_IsAbout43ArcsecPerCentury()
    {
        double arcsec = GeneralRelativity.PerihelionPrecessionArcsecPerCentury(
            MercurySemiMajorAxis, MercuryEccentricity, SolarMass, MercuryPeriodDays);
        Assert.AreEqual(43.0, arcsec, 1.0);
    }

    [TestMethod]
    public void LightDeflection_AtSunLimb_IsAbout1Point75Arcsec()
    {
        double sunRadius = 6.957e8; // m
        double radians = GeneralRelativity.LightDeflection(sunRadius, SolarMass);
        Assert.AreEqual(1.75, radians * RadToArcsec, 0.02);
    }

    [TestMethod]
    public void LightDeflection_ScalesInverselyWithImpactParameter()
    {
        double a = GeneralRelativity.LightDeflection(1.0e9, SolarMass);
        double b = GeneralRelativity.LightDeflection(2.0e9, SolarMass);
        Assert.AreEqual(a / 2.0, b, a * 1e-12);
    }

    [TestMethod]
    public void ShapiroDelay_IsPositive_AndGrowsWithMass()
    {
        double au = PhysicsConstants.AstronomicalUnit;
        double sunRadius = 6.957e8;
        double delay = GeneralRelativity.ShapiroDelay(au, au, sunRadius, SolarMass);

        Assert.IsTrue(delay > 0.0);
        // Earth–graze–Earth one-way excess delay is on the order of 0.1 ms.
        Assert.AreEqual(1.2e-4, delay, 0.3e-4);

        double heavier = GeneralRelativity.ShapiroDelay(au, au, sunRadius, 2.0 * SolarMass);
        Assert.IsTrue(heavier > delay);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentOutOfRangeException))]
    public void ShapiroDelay_RadiusBelowImpactParameter_Throws()
    {
        GeneralRelativity.ShapiroDelay(1.0e8, 1.0e12, 5.0e8, SolarMass);
    }
}
