using System;
using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Physics.Cosmology;

namespace NumericsTests;

[TestClass]
public class FlrwModelTests
{
    private const double Mpc = PhysicsConstants.Megaparsec;
    private const double SecondsPerYear = 3.15576e7; // Julian year

    private static FlrwModel Planck() => FlrwModel.Planck2018();

    [TestMethod]
    public void Planck2018_HasFlatGeometry_AndExpectedParameters()
    {
        var m = Planck();
        Assert.AreEqual(67.4, m.HubbleConstantKmSMpc, 1e-9);
        Assert.AreEqual(0.315, m.OmegaMatter, 1e-9);
        Assert.AreEqual(0.685, m.OmegaLambda, 1e-9);
        Assert.AreEqual(0.0, m.OmegaCurvature, 1e-9); // flat
    }

    [TestMethod]
    public void DimensionlessHubble_AtPresent_IsOne()
    {
        Assert.AreEqual(1.0, Planck().DimensionlessHubble(0.0), 1e-12);
    }

    [TestMethod]
    public void HubbleParameter_AtPresent_EqualsHubbleConstant()
    {
        var m = Planck();
        Assert.AreEqual(m.HubbleConstant, m.HubbleParameter(0.0), m.HubbleConstant * 1e-12);
        Assert.IsTrue(m.HubbleParameter(1.0) > m.HubbleParameter(0.0)); // expansion was faster
    }

    [TestMethod]
    public void CriticalDensity_Today_IsAbout8Point5e_minus27()
    {
        // ρ_c,0 ≈ 8.5×10⁻²⁷ kg/m³ for H₀ = 67.4
        Assert.AreEqual(8.5e-27, Planck().CriticalDensity(), 0.3e-27);
    }

    [TestMethod]
    public void AgeOfUniverse_IsAbout13Point8Gyr()
    {
        double ageYears = Planck().AgeOfUniverse() / SecondsPerYear;
        Assert.AreEqual(13.8e9, ageYears, 0.3e9);
    }

    [TestMethod]
    public void ComovingDistance_AtZ1_IsAbout3400Mpc()
    {
        double dcMpc = Planck().ComovingDistance(1.0) / Mpc;
        Assert.AreEqual(3402.0, dcMpc, 30.0);
    }

    [TestMethod]
    public void ComovingDistance_AtLowRedshift_FollowsHubbleLaw()
    {
        // D_C(z) → D_H·z as z → 0
        var m = Planck();
        double dc = m.ComovingDistance(0.001);
        Assert.AreEqual(m.HubbleDistance * 0.001, dc, m.HubbleDistance * 0.001 * 0.005);
    }

    [TestMethod]
    public void DistanceDuality_LuminosityEqualsOnePlusZSquaredTimesAngular()
    {
        // Etherington relation: D_L = (1+z)²·D_A
        var m = Planck();
        double z = 2.0;
        double dl = m.LuminosityDistance(z);
        double da = m.AngularDiameterDistance(z);
        Assert.AreEqual((1.0 + z) * (1.0 + z) * da, dl, dl * 1e-9);
    }

    [TestMethod]
    public void DistanceModulus_AtZ1_IsAbout44()
    {
        Assert.AreEqual(44.16, Planck().DistanceModulus(1.0), 0.1);
    }

    [TestMethod]
    public void ScaleFactor_IsInverseOfOnePlusZ()
    {
        Assert.AreEqual(0.5, Planck().ScaleFactor(1.0), 1e-12);
    }

    [TestMethod]
    public void LookbackTimePlusAge_EqualsAgeOfUniverse()
    {
        var m = Planck();
        double total = m.AgeOfUniverse();
        double sum = m.LookbackTime(1.0) + m.Age(1.0);
        Assert.AreEqual(total, sum, total * 1e-3);
    }

    [TestMethod]
    public void OpenUniverse_TransverseDistanceExceedsComoving()
    {
        // Ω_k > 0 (open) → sinh form makes D_M > D_C
        var open = new FlrwModel(70.0, 0.3, 0.0); // Ω_k = 0.7
        Assert.IsTrue(open.OmegaCurvature > 0.0);
        Assert.IsTrue(open.TransverseComovingDistance(2.0) > open.ComovingDistance(2.0));
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentOutOfRangeException))]
    public void ComovingDistance_NegativeRedshift_Throws()
    {
        Planck().ComovingDistance(-0.5);
    }
}
