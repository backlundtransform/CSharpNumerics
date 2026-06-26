using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Physics.Relativity;

namespace NumericsTests;

[TestClass]
public class RelativisticClockTests
{
    private const double EarthMass = PhysicsConstants.EarthMass;
    private const double SecondsPerDay = 86400.0;

    [TestMethod]
    public void GravitationalRate_HigherClockRunsFaster()
    {
        double rate = RelativisticClock.GravitationalFractionalRate(EarthMass, 6.371e6, 2.0e7);
        Assert.IsTrue(rate > 0); // upper clock ticks faster
    }

    [TestMethod]
    public void VelocityRate_IsNegative()
    {
        Assert.IsTrue(RelativisticClock.VelocityFractionalRate(7660.0) < 0); // moving clock runs slow
    }

    [TestMethod]
    public void GpsSatellite_NetClockGains_About38MicrosecondsPerDay()
    {
        // GPS: orbit radius ≈ 26 560 km, speed ≈ 3874 m/s
        double rate = RelativisticClock.OrbitingClockFractionalRate(
            EarthMass, groundRadius: 6.371e6, orbitRadius: 2.656e7, orbitSpeed: 3874.0);

        double microsecondsPerDay = rate * SecondsPerDay * 1e6;
        Assert.AreEqual(38.5, microsecondsPerDay, 1.5); // GPS clocks run ~38 µs/day fast
    }

    [TestMethod]
    public void IssOrbit_NetClockRunsSlow_VelocityDominates()
    {
        // Low, fast orbit (ISS ≈ 6 778 km, 7660 m/s): SR dilation beats GR blueshift
        double rate = RelativisticClock.OrbitingClockFractionalRate(
            EarthMass, groundRadius: 6.371e6, orbitRadius: 6.778e6, orbitSpeed: 7660.0);

        Assert.IsTrue(rate < 0);
    }
}
