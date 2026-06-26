using System;
using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Physics.Relativity;

namespace NumericsTests;

[TestClass]
public class SpecialRelativityTests
{
    private const double C = PhysicsConstants.SpeedOfLight;

    [TestMethod]
    public void LorentzFactor_AtRest_IsOne()
    {
        Assert.AreEqual(1.0, SpecialRelativity.LorentzFactor(0.0), 1e-15);
    }

    [TestMethod]
    public void LorentzFactor_At080c_IsFiveThirds()
    {
        // γ = 1/√(1 − 0.64) = 1/0.6 = 5/3
        Assert.AreEqual(5.0 / 3.0, SpecialRelativity.LorentzFactor(0.8 * C), 1e-9);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentOutOfRangeException))]
    public void LorentzFactor_AtLightSpeed_Throws()
    {
        SpecialRelativity.LorentzFactor(C);
    }

    [TestMethod]
    public void TimeDilation_MovingClockRunsSlow()
    {
        // 1 s of proper time at 0.8c → 5/3 s of coordinate time
        Assert.AreEqual(5.0 / 3.0, SpecialRelativity.TimeDilation(1.0, 0.8 * C), 1e-9);
    }

    [TestMethod]
    public void LengthContraction_At080c_IsThreeFifths()
    {
        Assert.AreEqual(0.6, SpecialRelativity.LengthContraction(1.0, 0.8 * C), 1e-9);
    }

    [TestMethod]
    public void KineticEnergy_LowSpeed_MatchesClassicalLimit()
    {
        // (γ − 1)mc² → ½mv² for v ≪ c
        double m = 2.0, v = 1000.0;
        double classical = 0.5 * m * v * v;
        Assert.AreEqual(classical, SpecialRelativity.KineticEnergy(m, v), 1.0);
    }

    [TestMethod]
    public void AddVelocities_HalfPlusHalf_IsBelowLight()
    {
        // (0.5c + 0.5c)/(1 + 0.25) = 0.8c
        Assert.AreEqual(0.8 * C, SpecialRelativity.AddVelocities(0.5 * C, 0.5 * C), 1.0);
    }

    [TestMethod]
    public void DopplerFactor_Receding_IsRedshift()
    {
        // β = 0.6 → √(0.4/1.6) = 0.5  (observed frequency halved)
        Assert.AreEqual(0.5, SpecialRelativity.DopplerFactor(0.6 * C), 1e-12);
    }

    [TestMethod]
    public void EnergyFromMomentum_MatchesTotalEnergy()
    {
        double m = 1.0, v = 0.9 * C;
        double p = SpecialRelativity.Momentum(m, v);
        double expected = SpecialRelativity.TotalEnergy(m, v);
        Assert.AreEqual(expected, SpecialRelativity.EnergyFromMomentum(m, p), expected * 1e-12);
    }

    [TestMethod]
    public void Rapidity_AddsLinearlyUnderVelocityComposition()
    {
        double u = 0.3 * C, v = 0.4 * C;
        double combined = SpecialRelativity.Rapidity(SpecialRelativity.AddVelocities(u, v));
        double sum = SpecialRelativity.Rapidity(u) + SpecialRelativity.Rapidity(v);
        Assert.AreEqual(sum, combined, 1e-9);
    }
}
