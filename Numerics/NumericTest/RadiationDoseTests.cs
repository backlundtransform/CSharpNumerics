using CSharpNumerics.Physics.Materials.Nuclear.Isotopes;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using DoseCalc = CSharpNumerics.Physics.Materials.Nuclear.RadiationDose.RadiationDose;

namespace NumericTest
{
    [TestClass]
    public class RadiationDoseTests
    {
        [TestMethod]
        public void DoseRate_InverseSquare()
        {
            double a = 1e9; // 1 GBq
            double d1 = DoseCalc.DoseRate(a, 1.0, Isotope.Cs137);
            double d2 = DoseCalc.DoseRate(a, 2.0, Isotope.Cs137);
            // At 2× distance, dose should be 1/4
            Assert.AreEqual(d1 / 4, d2, d1 * 0.001);
        }

        [TestMethod]
        public void DoseRate_PureBetaEmitter_IsZero()
        {
            // Sr-90 has no gamma
            double d = DoseCalc.DoseRate(1e9, 1.0, Isotope.Sr90);
            Assert.AreEqual(0.0, d);
        }

        [TestMethod]
        public void DoseRate_Cs137_ReasonableOrder()
        {
            // 1 TBq Cs-137 at 1 m: expect order of µSv/h to mSv/h
            double d = DoseCalc.DoseRate(1e12, 1.0, Isotope.Cs137);
            Assert.IsTrue(d > 1e-6 && d < 1, $"DoseRate: {d} Sv/h");
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void DoseRate_ZeroDistance_Throws()
        {
            DoseCalc.DoseRate(1e6, 0, Isotope.Cs137);
        }

        [TestMethod]
        public void GroundShineDose_Cs137()
        {
            // 1e6 Bq/m² for 1 hour
            double dose = DoseCalc.GroundShineDose(1e6, Isotope.Cs137, 3600);
            Assert.IsTrue(dose > 0);
            // Should be in the sub-µSv range for this level
            Assert.IsTrue(dose < 1e-3, $"Dose: {dose} Sv");
        }

        [TestMethod]
        public void GroundShineDose_Sr90_IsZero()
        {
            // Sr-90 has no ground-shine coefficient
            double dose = DoseCalc.GroundShineDose(1e6, Isotope.Sr90, 3600);
            Assert.AreEqual(0.0, dose);
        }

        [TestMethod]
        public void InhalationDose_Cs137()
        {
            // 1000 Bq/m³ for 1 hour at default breathing rate
            double dose = DoseCalc.InhalationDose(1000, 3600, Isotope.Cs137);
            // D = 1000 × 3.33e-4 × 3600 × 3.9e-8 ≈ 4.67e-5 Sv
            double expected = 1000 * DoseCalc.DefaultBreathingRate * 3600 * Isotope.Cs137.InhalationDoseCoeff;
            Assert.AreEqual(expected, dose, expected * 0.001);
        }

        [TestMethod]
        public void InhalationDose_CustomBreathingRate()
        {
            double customRate = 5e-4; // higher rate
            double dose = DoseCalc.InhalationDose(1000, customRate, 3600, Isotope.Cs137);
            double expected = 1000 * customRate * 3600 * Isotope.Cs137.InhalationDoseCoeff;
            Assert.AreEqual(expected, dose, expected * 0.001);
        }

        [TestMethod]
        public void InhalationDose_ZeroCoeff_IsZero()
        {
            // Ba-137m has no inhalation dose coefficient set
            double dose = DoseCalc.InhalationDose(1e6, 3600, Isotope.Ba137m);
            Assert.AreEqual(0.0, dose);
        }

        [TestMethod]
        public void InhalationDose_ProportionalToConcentration()
        {
            double d1 = DoseCalc.InhalationDose(100, 3600, Isotope.I131);
            double d2 = DoseCalc.InhalationDose(200, 3600, Isotope.I131);
            Assert.AreEqual(d1 * 2, d2, d1 * 0.001);
        }

        [TestMethod]
        public void InhalationDose_ProportionalToTime()
        {
            double d1 = DoseCalc.InhalationDose(100, 3600, Isotope.I131);
            double d2 = DoseCalc.InhalationDose(100, 7200, Isotope.I131);
            Assert.AreEqual(d1 * 2, d2, d1 * 0.001);
        }
    }
}
