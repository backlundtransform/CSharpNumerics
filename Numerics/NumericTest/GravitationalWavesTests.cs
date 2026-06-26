using System;
using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Physics.Relativity;

namespace NumericsTests;

[TestClass]
public class GravitationalWavesTests
{
    private const double SolarMass = PhysicsConstants.SolarMass;

    [TestMethod]
    public void ChirpMass_EqualMasses_IsMassOverFifthRootOfTwo()
    {
        // M_c = (m²)^(3/5)/(2m)^(1/5) = m / 2^(1/5)
        double m = 10.0;
        Assert.AreEqual(m / Math.Pow(2.0, 0.2), GravitationalWaves.ChirpMass(m, m), 1e-9);
    }

    [TestMethod]
    public void ChirpMass_IsSymmetric()
    {
        Assert.AreEqual(
            GravitationalWaves.ChirpMass(36 * SolarMass, 29 * SolarMass),
            GravitationalWaves.ChirpMass(29 * SolarMass, 36 * SolarMass), 1.0);
    }

    [TestMethod]
    public void GravitationalWaveFrequency_IsTwiceOrbitalFrequency()
    {
        double m1 = 30 * SolarMass, m2 = 30 * SolarMass, a = 1.0e8;
        Assert.AreEqual(
            2.0 * GravitationalWaves.OrbitalFrequency(m1, m2, a),
            GravitationalWaves.GravitationalWaveFrequency(m1, m2, a), 1e-9);
    }

    [TestMethod]
    public void MergerTime_TwoNeutronStars_MatchesPetersFormula()
    {
        // 1.4 + 1.4 M☉ at a₀ = 1×10⁹ m → ≈ 3.69×10¹⁵ s
        double m = 1.4 * SolarMass;
        double t = GravitationalWaves.MergerTime(m, m, 1.0e9);
        Assert.AreEqual(3.687e15, t, 0.05e15);
    }

    [TestMethod]
    public void MergerTime_ScalesWithSeparationToTheFourth()
    {
        double m = 1.4 * SolarMass;
        double t1 = GravitationalWaves.MergerTime(m, m, 1.0e9);
        double t2 = GravitationalWaves.MergerTime(m, m, 2.0e9);
        Assert.AreEqual(16.0 * t1, t2, t1 * 1e-9); // a⁴ scaling
    }

    [TestMethod]
    public void Luminosity_ScalesWithSeparationToMinusFifth()
    {
        double m1 = 30 * SolarMass, m2 = 30 * SolarMass;
        double l1 = GravitationalWaves.Luminosity(m1, m2, 1.0e8);
        double l2 = GravitationalWaves.Luminosity(m1, m2, 2.0e8);
        Assert.IsTrue(l1 > 0);
        Assert.AreEqual(l1 / 32.0, l2, l1 * 1e-9); // a⁻⁵ scaling
    }

    [TestMethod]
    public void StrainAmplitude_IsPositive_AndScalesInverselyWithDistance()
    {
        double mc = GravitationalWaves.ChirpMass(30 * SolarMass, 30 * SolarMass);
        double h1 = GravitationalWaves.StrainAmplitude(mc, distance: 1.0e24, gravitationalWaveFrequency: 100.0);
        double h2 = GravitationalWaves.StrainAmplitude(mc, distance: 2.0e24, gravitationalWaveFrequency: 100.0);

        Assert.IsTrue(h1 > 0);
        Assert.AreEqual(h1 / 2.0, h2, h1 * 1e-12);
    }
}
