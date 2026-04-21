using CSharpNumerics.Physics.Environmental.Water;
using CSharpNumerics.Physics.Materials.Water;
using CSharpNumerics.Physics.Materials.Water.Enums;
using System;
using System.Linq;

namespace NumericsTests
{
    [TestClass]
    public class WaterContaminationPhysicsTests
    {
        // ════════════════════════════════════════════════════════════════
        //  ContaminantLibrary
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void ContaminantLibrary_GetByName_ReturnsCorrectContaminant()
        {
            var cs = ContaminantLibrary.Get("Cs137");
            Assert.AreEqual("Cs137", cs.Name);
            Assert.AreEqual(ContaminantType.Radioactive, cs.Type);
        }

        [TestMethod]
        public void ContaminantLibrary_CaseInsensitive()
        {
            var benzene = ContaminantLibrary.Get("benzene");
            Assert.AreEqual("Benzene", benzene.Name);
        }

        [TestMethod]
        [ExpectedException(typeof(System.Collections.Generic.KeyNotFoundException))]
        public void ContaminantLibrary_UnknownName_Throws()
        {
            ContaminantLibrary.Get("NonExistent");
        }

        [TestMethod]
        public void ContaminantLibrary_TryGet_ReturnsFalseForUnknown()
        {
            bool found = ContaminantLibrary.TryGet("NonExistent", out _);
            Assert.IsFalse(found);
        }

        [TestMethod]
        public void ContaminantLibrary_All_ReturnsPreloaded()
        {
            var all = ContaminantLibrary.All();
            // 3 radioactive + 6 chemical + 2 biological + 1 thermal = 12
            Assert.IsTrue(all.Count >= 12, $"Expected ≥12 pre-loaded contaminants, got {all.Count}");
        }

        [TestMethod]
        public void ContaminantLibrary_Register_CustomContaminant()
        {
            var custom = new AquaticContaminant("TestPollutant", ContaminantType.Chemical,
                halfLifeSeconds: 3600, partitionCoefficient: 10,
                toxicityThresholdMgL: 0.5, lethalThresholdMgL: 50);
            ContaminantLibrary.Register(custom);

            var retrieved = ContaminantLibrary.Get("TestPollutant");
            Assert.AreEqual("TestPollutant", retrieved.Name);
            Assert.AreEqual(3600, retrieved.HalfLifeSeconds);
        }

        [TestMethod]
        public void ContaminantLibrary_AllContaminantsHaveValidThresholds()
        {
            foreach (var c in ContaminantLibrary.All())
            {
                Assert.IsTrue(c.ToxicityThresholdMgL > 0,
                    $"{c.Name}: toxicity threshold must be > 0");
                Assert.IsTrue(c.LethalThresholdMgL >= c.ToxicityThresholdMgL,
                    $"{c.Name}: lethal threshold must be ≥ toxicity threshold");
            }
        }

        // ════════════════════════════════════════════════════════════════
        //  AquaticContaminant — decay constant
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void Cs137_DecayConstant_MatchesExpected()
        {
            var cs = AquaticContaminant.Cs137;
            // λ = ln(2) / (30.17 × 365.25 × 86400) ≈ 7.28 × 10⁻¹⁰ s⁻¹
            double expected = Math.Log(2) / (30.17 * 365.25 * 86400);
            Assert.AreEqual(expected, cs.DecayConstant, expected * 1e-6,
                "Cs-137 decay constant should match ln(2)/t½");
        }

        [TestMethod]
        public void ConservativeTracer_DecayConstantIsZero()
        {
            // Cyanide has halfLife = 0 → conservative
            var cn = AquaticContaminant.Cyanide;
            Assert.AreEqual(0.0, cn.DecayConstant);
            Assert.IsTrue(cn.IsConservative);
        }

        [TestMethod]
        public void Benzene_IsNotConservative()
        {
            Assert.IsFalse(AquaticContaminant.Benzene.IsConservative);
            Assert.IsTrue(AquaticContaminant.Benzene.DecayConstant > 0);
        }

        // ════════════════════════════════════════════════════════════════
        //  Manning's equation
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void Manning_RectangularChannel_VelocityMatchesHandbook()
        {
            // Rectangular channel: W=10 m, H=2 m, n=0.035, S=0.001
            double width = 10, depth = 2, n = 0.035, slope = 0.001;
            double Rh = ManningEquation.RectangularHydraulicRadius(width, depth);

            // Rh = (10×2)/(10+2×2) = 20/14 ≈ 1.4286 m
            Assert.AreEqual(20.0 / 14.0, Rh, 1e-6);

            double u = ManningEquation.Velocity(n, Rh, slope);
            // u = (1/0.035) × 1.4286^(2/3) × 0.001^(1/2)
            //   ≈ 28.571 × 1.267 × 0.03162 ≈ 1.145 m/s
            Assert.IsTrue(u > 1.0 && u < 1.3,
                $"Expected velocity ~1.14 m/s for standard channel, got {u:F3}");
        }

        [TestMethod]
        public void Manning_DischargeEqualsVelocityTimesArea()
        {
            double n = 0.035, slope = 0.001;
            double width = 10, depth = 2;
            double Rh = ManningEquation.RectangularHydraulicRadius(width, depth);
            double area = width * depth;

            double u = ManningEquation.Velocity(n, Rh, slope);
            double Q = ManningEquation.Discharge(n, Rh, slope, area);

            Assert.AreEqual(u * area, Q, 1e-10, "Q must equal u × A");
        }

        [TestMethod]
        public void Manning_ZeroSlope_ReturnsZero()
        {
            double u = ManningEquation.Velocity(0.035, 1.0, 0);
            Assert.AreEqual(0, u);
        }

        [TestMethod]
        public void Manning_ZeroManningN_ReturnsZero()
        {
            double u = ManningEquation.Velocity(0, 1.0, 0.001);
            Assert.AreEqual(0, u);
        }

        [TestMethod]
        public void Manning_TrapezoidalHydraulicRadius()
        {
            // Trapezoidal: b=5 m, d=2 m, z=1.5 (1.5H:1V)
            double Rh = ManningEquation.TrapezoidalHydraulicRadius(5, 2, 1.5);
            // A = (5 + 1.5×2)×2 = (5+3)×2 = 16 m²
            // P = 5 + 2×2×√(1+1.5²) = 5 + 4×√3.25 = 5 + 7.211 = 12.211 m
            // Rh = 16/12.211 ≈ 1.310 m
            double expectedA = (5 + 1.5 * 2) * 2;
            double expectedP = 5 + 2 * 2 * Math.Sqrt(1 + 1.5 * 1.5);
            double expected = expectedA / expectedP;
            Assert.AreEqual(expected, Rh, 1e-6);
        }

        [TestMethod]
        public void Manning_HydraulicRadius_Generic()
        {
            double Rh = ManningEquation.HydraulicRadius(20, 14);
            Assert.AreEqual(20.0 / 14.0, Rh, 1e-10);
        }

        // ════════════════════════════════════════════════════════════════
        //  Longitudinal dispersion
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void Fischer_WideShallowRiver_HigherEL()
        {
            double u = 1.0, uStar = 0.1;
            // Wide shallow: W=100 m, H=1 m
            double EL_wide = LongitudinalDispersion.FischerCoefficient(u, 100, 1, uStar);
            // Narrow deep: W=10 m, H=10 m
            double EL_narrow = LongitudinalDispersion.FischerCoefficient(u, 10, 10, uStar);

            Assert.IsTrue(EL_wide > EL_narrow,
                $"Wide shallow river should have higher EL ({EL_wide:F1}) than narrow deep ({EL_narrow:F1})");
        }

        [TestMethod]
        public void Fischer_ZeroInputs_ReturnsZero()
        {
            Assert.AreEqual(0, LongitudinalDispersion.FischerCoefficient(0, 10, 1, 0.1));
            Assert.AreEqual(0, LongitudinalDispersion.FischerCoefficient(1, 0, 1, 0.1));
            Assert.AreEqual(0, LongitudinalDispersion.FischerCoefficient(1, 10, 0, 0.1));
            Assert.AreEqual(0, LongitudinalDispersion.FischerCoefficient(1, 10, 1, 0));
        }

        [TestMethod]
        public void ShearVelocity_PositiveForValidInput()
        {
            double uStar = LongitudinalDispersion.ShearVelocity(1.4, 0.001);
            // u* = √(9.81 × 1.4 × 0.001) ≈ 0.117 m/s
            double expected = Math.Sqrt(9.80665 * 1.4 * 0.001);
            Assert.AreEqual(expected, uStar, 1e-6);
        }

        [TestMethod]
        public void ShearVelocity_ZeroSlope_ReturnsZero()
        {
            Assert.AreEqual(0, LongitudinalDispersion.ShearVelocity(1.0, 0));
        }

        [TestMethod]
        public void DecayConstant_Cs137_MatchesHalfLife()
        {
            double halfLife = 30.17 * 365.25 * 86400;
            double lambda = LongitudinalDispersion.DecayConstant(halfLife);
            Assert.AreEqual(Math.Log(2) / halfLife, lambda, 1e-15);
        }

        [TestMethod]
        public void DecayConstant_ZeroHalfLife_ReturnsZero()
        {
            Assert.AreEqual(0, LongitudinalDispersion.DecayConstant(0));
        }

        [TestMethod]
        public void DecayConstant_InfiniteHalfLife_ReturnsZero()
        {
            Assert.AreEqual(0, LongitudinalDispersion.DecayConstant(double.PositiveInfinity));
        }

        [TestMethod]
        public void RetardationFactor_NoAdsorption_ReturnsOne()
        {
            // Kd = 0 → Rf = 1
            double Rf = LongitudinalDispersion.RetardationFactor(1600, 0, 0.4);
            Assert.AreEqual(1.0, Rf, 1e-10);
        }

        [TestMethod]
        public void RetardationFactor_WithAdsorption_GreaterThanOne()
        {
            // ρs=1600 kg/m³, Kd=1000 L/kg, n=0.4
            // Rf = 1 + (1600/1000) × 1000 / 0.4 = 1 + 1.6 × 2500 = 4001
            double Rf = LongitudinalDispersion.RetardationFactor(1600, 1000, 0.4);
            Assert.IsTrue(Rf > 1, $"Rf should be > 1 with adsorption, got {Rf}");
            double expected = 1.0 + (1600.0 / 1000.0) * 1000.0 / 0.4;
            Assert.AreEqual(expected, Rf, 1e-6);
        }

        // ════════════════════════════════════════════════════════════════
        //  Mixing zone model
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void TributaryMixing_EqualFlows_AveragesConcentration()
        {
            // Main: 0 mg/L at 10 m³/s, Tributary: 100 mg/L at 10 m³/s → 50 mg/L
            double result = MixingZoneModel.TributaryMixing(0, 10, 100, 10);
            Assert.AreEqual(50.0, result, 1e-10);
        }

        [TestMethod]
        public void TributaryMixing_DominantMain_DilutesTributary()
        {
            // Main: 0 mg/L at 100 m³/s, Tributary: 1000 mg/L at 1 m³/s
            double result = MixingZoneModel.TributaryMixing(0, 100, 1000, 1);
            // Expected: 1000/101 ≈ 9.9 mg/L
            Assert.AreEqual(1000.0 / 101.0, result, 1e-10);
        }

        [TestMethod]
        public void TributaryMixing_ZeroDischarge_ReturnsZero()
        {
            double result = MixingZoneModel.TributaryMixing(50, 0, 100, 0);
            Assert.AreEqual(0, result);
        }

        [TestMethod]
        public void MixingLength_Positive()
        {
            // u=1 m/s, W=20 m, Et=0.1 m²/s
            double Lmix = MixingZoneModel.MixingLength(1, 20, 0.1);
            // Lmix = 0.4 × 1 × 400 / 0.1 = 1600 m
            Assert.AreEqual(1600.0, Lmix, 1e-6);
        }

        [TestMethod]
        public void MixingLength_ZeroVelocity_ReturnsZero()
        {
            Assert.AreEqual(0, MixingZoneModel.MixingLength(0, 20, 0.1));
        }

        [TestMethod]
        public void MixingLength_ZeroWidth_ReturnsZero()
        {
            Assert.AreEqual(0, MixingZoneModel.MixingLength(1, 0, 0.1));
        }
    }
}
