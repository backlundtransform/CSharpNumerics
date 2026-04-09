using CSharpNumerics.Engines.Exoplanet.Data;
using CSharpNumerics.Engines.Exoplanet.Enums;
using CSharpNumerics.Physics.Astro;
using CSharpNumerics.Physics.Constants;
using System;
using System.Linq;

namespace NumericTest;

[TestClass]
public class TransitPhysicsTests
{
    private const double SolarRadiusMeters = 6.957e8;

    #region TransitGeometry Tests

    [TestMethod]
    public void ImpactParameter_EdgeOn_ShouldBeZero()
    {
        // Edge-on: inclination = 90 degrees → b = 0
        double b = TransitGeometry.ImpactParameter(10.0, Math.PI / 2);
        Assert.AreEqual(0.0, b, 1e-10);
    }

    [TestMethod]
    public void ImpactParameter_Inclined_ShouldBePositive()
    {
        // a/R★ = 10, i = 85° → b = 10·cos(85°) ≈ 0.872
        double b = TransitGeometry.ImpactParameter(10.0, 85.0 * Math.PI / 180.0);
        Assert.AreEqual(10.0 * Math.Cos(85.0 * Math.PI / 180.0), b, 1e-10);
        Assert.IsTrue(b > 0);
    }

    [TestMethod]
    public void TransitProbability_ShouldScaleWithRadii()
    {
        // P = (R★ + Rp) / a
        double prob = TransitGeometry.TransitProbability(215.0, 1.0, 0.1);
        Assert.AreEqual(1.1 / 215.0, prob, 1e-10);
    }

    [TestMethod]
    public void TransitDuration_HotJupiter_ShouldBeReasonable()
    {
        // Hot Jupiter: P = 3 days, a/R★ ≈ 8.8, k = 0.1, i = 90°
        double duration = TransitGeometry.TransitDuration(3.0, 8.8, 0.1, Math.PI / 2);

        // Typical hot Jupiter transit ~2–4 hours = 0.08–0.17 days
        Assert.IsTrue(duration > 0.05 && duration < 0.25,
            $"Transit duration {duration} days outside expected range.");
    }

    [TestMethod]
    public void TransitDuration_NoTransit_ShouldBeZero()
    {
        // Very low inclination → no transit
        double duration = TransitGeometry.TransitDuration(3.0, 20.0, 0.1, 0.1);
        Assert.AreEqual(0.0, duration, 1e-10);
    }

    [TestMethod]
    public void IngressDuration_ShouldBeLessThanTotalDuration()
    {
        double t14 = TransitGeometry.TransitDuration(3.0, 8.8, 0.1, Math.PI / 2);
        double tIngress = TransitGeometry.IngressDuration(3.0, 8.8, 0.1, Math.PI / 2);

        Assert.IsTrue(tIngress > 0, "Ingress duration should be positive.");
        Assert.IsTrue(tIngress < t14 / 2.0, "Ingress should be less than half total duration.");
    }

    [TestMethod]
    public void ContactTimes_ShouldBeOrdered()
    {
        double epoch = 100.0;
        var (t1, t2, t3, t4) = TransitGeometry.ContactTimes(epoch, 3.0, 8.8, 0.1, Math.PI / 2);

        Assert.IsTrue(t1 < t2, $"T1={t1} should be < T2={t2}");
        Assert.IsTrue(t2 < t3, $"T2={t2} should be < T3={t3}");
        Assert.IsTrue(t3 < t4, $"T3={t3} should be < T4={t4}");

        // Symmetric around epoch
        Assert.AreEqual(epoch, (t1 + t4) / 2.0, 1e-10);
        Assert.AreEqual(epoch, (t2 + t3) / 2.0, 1e-10);
    }

    #endregion

    #region LimbDarkening Tests

    [TestMethod]
    public void Linear_AtCenter_ShouldBeOne()
    {
        // μ = 1 (center): I = 1 − u1·(1 − 1) = 1
        Assert.AreEqual(1.0, LimbDarkening.Linear(1.0, 0.5), 1e-10);
    }

    [TestMethod]
    public void Linear_AtLimb_ShouldBeDimmed()
    {
        // μ = 0 (limb): I = 1 − u1
        Assert.AreEqual(0.5, LimbDarkening.Linear(0.0, 0.5), 1e-10);
    }

    [TestMethod]
    public void Quadratic_AtCenter_ShouldBeOne()
    {
        Assert.AreEqual(1.0, LimbDarkening.Quadratic(1.0, 0.3, 0.2), 1e-10);
    }

    [TestMethod]
    public void Quadratic_AtLimb_ShouldEqual_1MinusU1MinusU2()
    {
        // μ = 0: I = 1 − u1 − u2
        double u1 = 0.3, u2 = 0.2;
        Assert.AreEqual(1.0 - u1 - u2, LimbDarkening.Quadratic(0.0, u1, u2), 1e-10);
    }

    [TestMethod]
    public void NonlinearFourParam_AtCenter_ShouldBeOne()
    {
        Assert.AreEqual(1.0, LimbDarkening.NonlinearFourParam(1.0, 0.1, 0.2, 0.1, 0.1), 1e-10);
    }

    [TestMethod]
    public void IntensityProfile_Uniform_ShouldBeAllOnes()
    {
        var mu = new double[] { 0.0, 0.25, 0.5, 0.75, 1.0 };
        var profile = LimbDarkening.IntensityProfile(LimbDarkeningModel.Uniform, new double[0], mu);

        foreach (var val in profile)
            Assert.AreEqual(1.0, val, 1e-10);
    }

    [TestMethod]
    public void IntensityProfile_Linear_ShouldDecreaseTowardLimb()
    {
        var mu = new double[] { 1.0, 0.75, 0.5, 0.25, 0.0 };
        var profile = LimbDarkening.IntensityProfile(LimbDarkeningModel.Linear, new[] { 0.6 }, mu);

        for (int i = 0; i < profile.Length - 1; i++)
        {
            Assert.IsTrue(profile[i] >= profile[i + 1],
                $"Intensity should decrease toward limb: I[{i}]={profile[i]} vs I[{i + 1}]={profile[i + 1]}");
        }
    }

    [TestMethod]
    public void IntensityProfile_Quadratic_ShouldMatchDirectCalculation()
    {
        double u1 = 0.4, u2 = 0.2;
        var mu = new double[] { 0.3, 0.7 };
        var profile = LimbDarkening.IntensityProfile(LimbDarkeningModel.Quadratic, new[] { u1, u2 }, mu);

        Assert.AreEqual(LimbDarkening.Quadratic(0.3, u1, u2), profile[0], 1e-10);
        Assert.AreEqual(LimbDarkening.Quadratic(0.7, u1, u2), profile[1], 1e-10);
    }

    #endregion

    #region KeplerOrbit Tests

    [TestMethod]
    public void TrueAnomaly_CircularOrbit_ShouldEqualMeanAnomaly()
    {
        // For e = 0, true anomaly = mean anomaly
        double M = 1.5;
        double nu = KeplerOrbit.TrueAnomaly(M, 0.0);
        Assert.AreEqual(M, nu, 1e-10);
    }

    [TestMethod]
    public void TrueAnomaly_EccentricOrbit_ShouldConverge()
    {
        // e = 0.5, M = π/2: known solution exists
        double nu = KeplerOrbit.TrueAnomaly(Math.PI / 2, 0.5);

        // ν should be > M for e > 0 when M is in first half
        Assert.IsTrue(nu > Math.PI / 2);
        Assert.IsTrue(nu < Math.PI);
    }

    [TestMethod]
    public void TrueAnomaly_AtPeriapsis_ShouldBeZero()
    {
        double nu = KeplerOrbit.TrueAnomaly(0.0, 0.3);
        Assert.AreEqual(0.0, nu, 1e-10);
    }

    [TestMethod]
    public void TrueAnomaly_AtApoapsis_ShouldBePi()
    {
        double nu = KeplerOrbit.TrueAnomaly(Math.PI, 0.3);
        Assert.AreEqual(Math.PI, nu, 1e-8);
    }

    [TestMethod]
    public void SemiMajorAxis_Earth_ShouldBeAbout1AU()
    {
        // Earth: P ≈ 365.25 days, M★ = 1 M☉
        double aAU = KeplerOrbit.SemiMajorAxisAU(365.25, 1.0);
        Assert.AreEqual(1.0, aAU, 0.01);
    }

    [TestMethod]
    public void SemiMajorAxis_HotJupiter_ShouldBeSmall()
    {
        // P = 3 days, M★ = 1 M☉ → a ≈ 0.04 AU
        double aAU = KeplerOrbit.SemiMajorAxisAU(3.0, 1.0);
        Assert.IsTrue(aAU > 0.03 && aAU < 0.06,
            $"Semi-major axis {aAU} AU outside expected range for hot Jupiter.");
    }

    [TestMethod]
    public void OrbitalVelocity_ShouldBe2PiAOverP()
    {
        double a = 1.0e11;
        double P = 1.0e7;
        double v = KeplerOrbit.OrbitalVelocity(a, P);
        Assert.AreEqual(2.0 * Math.PI * a / P, v, 1e-5);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentOutOfRangeException))]
    public void TrueAnomaly_EccentricityOutOfRange_ShouldThrow()
    {
        KeplerOrbit.TrueAnomaly(1.0, 1.0);
    }

    #endregion

    #region TransitModel Tests

    [TestMethod]
    public void TransitModel_Uniform_OutOfTransit_ShouldBeOne()
    {
        var model = new TransitModel();
        // Times well away from any transit (epoch=10 → transits at t=1,4,7,10,...)
        double[] times = { 2.0, 2.5, 3.0, 3.5 };
        double period = 3.0;
        double epoch = 10.0;
        double k = 0.1;
        double aOverRstar = 8.8;
        double inclination = Math.PI / 2;

        double[] flux = model.Evaluate(times, period, epoch, k, aOverRstar, inclination,
            LimbDarkeningModel.Uniform, new double[0]);

        foreach (var f in flux)
            Assert.AreEqual(1.0, f, 1e-6, $"Out-of-transit flux should be 1.0, got {f}");
    }

    [TestMethod]
    public void TransitModel_Uniform_AtMidTransit_ShouldShowDip()
    {
        var model = new TransitModel();
        double period = 3.0;
        double epoch = 0.0;
        double k = 0.1;
        double aOverRstar = 8.8;
        double inclination = Math.PI / 2;

        // At mid-transit (t = epoch), z = 0 for edge-on orbit
        double[] times = { 0.0 };

        double[] flux = model.Evaluate(times, period, epoch, k, aOverRstar, inclination,
            LimbDarkeningModel.Uniform, new double[0]);

        // Expected depth = k² = 0.01
        double expectedDepth = k * k;
        double actualDepth = 1.0 - flux[0];

        Assert.AreEqual(expectedDepth, actualDepth, 0.001,
            $"Uniform transit depth should be ~{expectedDepth}, got {actualDepth}");
    }

    [TestMethod]
    public void TransitModel_Uniform_SymmetricAroundEpoch()
    {
        var model = new TransitModel();
        double period = 3.0;
        double epoch = 0.0;
        double k = 0.1;
        double aOverRstar = 8.8;
        double inclination = Math.PI / 2;

        double[] times = { -0.05, -0.03, -0.01, 0.01, 0.03, 0.05 };

        double[] flux = model.Evaluate(times, period, epoch, k, aOverRstar, inclination,
            LimbDarkeningModel.Uniform, new double[0]);

        // flux[0] ≈ flux[5], flux[1] ≈ flux[4], flux[2] ≈ flux[3]
        Assert.AreEqual(flux[0], flux[5], 1e-6, "Transit should be symmetric.");
        Assert.AreEqual(flux[1], flux[4], 1e-6, "Transit should be symmetric.");
        Assert.AreEqual(flux[2], flux[3], 1e-6, "Transit should be symmetric.");
    }

    [TestMethod]
    public void TransitModel_Quadratic_ShouldHaveDeeperCenterThanUniform()
    {
        var model = new TransitModel();
        double period = 3.0;
        double epoch = 0.0;
        double k = 0.1;
        double aOverRstar = 8.8;
        double inclination = Math.PI / 2;
        double[] times = { 0.0 };

        double[] fluxUniform = model.Evaluate(times, period, epoch, k, aOverRstar, inclination,
            LimbDarkeningModel.Uniform, new double[0]);

        double[] fluxLD = model.Evaluate(times, period, epoch, k, aOverRstar, inclination,
            LimbDarkeningModel.Quadratic, new[] { 0.4, 0.2 });

        // LD transit at center should be deeper because center is brighter
        double depthUniform = 1.0 - fluxUniform[0];
        double depthLD = 1.0 - fluxLD[0];

        Assert.IsTrue(depthLD > depthUniform,
            $"LD depth ({depthLD}) should exceed uniform depth ({depthUniform}).");
    }

    [TestMethod]
    public void TransitModel_WithStellarProperties_ShouldProduceTransit()
    {
        var model = new TransitModel();
        var star = new StellarProperties(5778, 1.0, 1.0, 4.44, 0.0, SpectralType.G);
        var p = new TransitParameters
        {
            Period = 3.0,
            Epoch = 0.0,
            Depth = 0.01,
            Duration = 0.12,
            RadiusRatio = 0.1,
            ImpactParameter = 0.0,
            IngressDuration = 0.015
        };

        // Generate transit around epoch
        int n = 200;
        double[] times = new double[n];
        for (int i = 0; i < n; i++)
            times[i] = -0.1 + 0.001 * i; // ±0.1 days around epoch

        double[] flux = model.Evaluate(times, p, LimbDarkeningModel.Quadratic, new[] { 0.4, 0.2 }, star);

        // Verify transit
        double minFlux = flux.Min();
        double maxFlux = flux.Max();

        Assert.IsTrue(minFlux < 0.995, $"Should show transit dip, min flux = {minFlux}");
        Assert.AreEqual(1.0, maxFlux, 0.001, "Out-of-transit flux should be ~1.0");
    }

    [TestMethod]
    public void TransitModel_LargerPlanet_ShouldHaveDeeperTransit()
    {
        var model = new TransitModel();
        double period = 3.0;
        double epoch = 0.0;
        double aOverRstar = 8.8;
        double inclination = Math.PI / 2;
        double[] times = { 0.0 };

        double[] fluxSmall = model.Evaluate(times, period, epoch, 0.05, aOverRstar, inclination,
            LimbDarkeningModel.Uniform, new double[0]);

        double[] fluxLarge = model.Evaluate(times, period, epoch, 0.15, aOverRstar, inclination,
            LimbDarkeningModel.Uniform, new double[0]);

        Assert.IsTrue(fluxLarge[0] < fluxSmall[0],
            "Larger planet should produce deeper transit.");
    }

    [TestMethod]
    public void TransitModel_Periodic_ShouldRepeat()
    {
        var model = new TransitModel();
        double period = 3.0;
        double epoch = 0.0;
        double k = 0.1;
        double aOverRstar = 8.8;
        double inclination = Math.PI / 2;

        // At mid-transit and one period later
        double[] times = { 0.0, 3.0, 6.0 };

        double[] flux = model.Evaluate(times, period, epoch, k, aOverRstar, inclination,
            LimbDarkeningModel.Uniform, new double[0]);

        Assert.AreEqual(flux[0], flux[1], 1e-6, "Transit should repeat at period.");
        Assert.AreEqual(flux[0], flux[2], 1e-6, "Transit should repeat at 2×period.");
    }

    #endregion
}
