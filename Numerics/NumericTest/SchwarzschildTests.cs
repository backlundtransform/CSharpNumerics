using System;
using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Physics.Relativity;

namespace NumericsTests;

[TestClass]
public class SchwarzschildTests
{
    private const double C = PhysicsConstants.SpeedOfLight;
    private const double SolarMass = PhysicsConstants.SolarMass;

    [TestMethod]
    public void Radius_OfSun_IsAboutTwoPointNineKilometres()
    {
        // r_s(Sun) ≈ 2953 m
        Assert.AreEqual(2953.0, Schwarzschild.Radius(SolarMass), 5.0);
    }

    [TestMethod]
    public void PhotonSphere_And_Isco_AreMultiplesOfSchwarzschildRadius()
    {
        double rs = Schwarzschild.Radius(SolarMass);
        Assert.AreEqual(1.5 * rs, Schwarzschild.PhotonSphereRadius(SolarMass), rs * 1e-12);
        Assert.AreEqual(3.0 * rs, Schwarzschild.Isco(SolarMass), rs * 1e-12);
    }

    [TestMethod]
    public void TimeDilationFactor_FarFromMass_ApproachesOne()
    {
        double factor = Schwarzschild.TimeDilationFactor(1.0e12, SolarMass);
        Assert.IsTrue(factor < 1.0);
        Assert.AreEqual(1.0, factor, 1e-5);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentOutOfRangeException))]
    public void TimeDilationFactor_AtHorizon_Throws()
    {
        double rs = Schwarzschild.Radius(SolarMass);
        Schwarzschild.TimeDilationFactor(rs, SolarMass);
    }

    [TestMethod]
    public void Redshift_IsPositive_AndSmallFarFromMass()
    {
        double z = Schwarzschild.Redshift(1.0e9, SolarMass);
        Assert.IsTrue(z > 0.0);
        Assert.IsTrue(z < 1e-5);
    }

    [TestMethod]
    public void ClockRateRatio_LowerClockRunsSlow()
    {
        // Clock deeper in the well (smaller r) ticks slower than the higher one.
        double ratio = Schwarzschild.ClockRateRatio(1.0e7, 1.0e9, SolarMass);
        Assert.IsTrue(ratio < 1.0);
    }

    [TestMethod]
    public void EscapeVelocity_AtSchwarzschildRadius_EqualsSpeedOfLight()
    {
        double rs = Schwarzschild.Radius(SolarMass);
        Assert.AreEqual(C, Schwarzschild.EscapeVelocity(rs, SolarMass), 1.0);
    }
}
