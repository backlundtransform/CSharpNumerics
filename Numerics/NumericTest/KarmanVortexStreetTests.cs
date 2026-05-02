using CSharpNumerics.Physics.FluidDynamics;
using System;

namespace NumericTest
{
    [TestClass]
    public class KarmanVortexStreetTests
    {
        #region Shedding Frequency

        [TestMethod]
        public void VortexSheddingFrequency_StandardValues()
        {
            // St = 0.2, U = 10 m/s, D = 0.1 m → f = 0.2 * 10 / 0.1 = 20 Hz
            double f = 0.2.VortexSheddingFrequency(10.0, 0.1);
            Assert.AreEqual(20.0, f, 1e-10);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void VortexSheddingFrequency_ZeroDiameter_Throws()
        {
            0.2.VortexSheddingFrequency(10.0, 0.0);
        }

        [TestMethod]
        public void VortexSheddingFrequencyFromReynolds_Re1000()
        {
            // Re = 1000 → St ≈ 0.198 * (1 - 19.7/1000) ≈ 0.19410
            // f = St * U / D
            double U = 5.0, D = 0.05;
            double f = 1000.0.VortexSheddingFrequencyFromReynolds(U, D);
            double expectedSt = 0.198 * (1.0 - 19.7 / 1000.0);
            double expectedF = expectedSt * U / D;
            Assert.AreEqual(expectedF, f, 1e-10);
        }

        #endregion

        #region Roshko Strouhal Correlation

        [TestMethod]
        public void RoshkoStrouhalNumber_Re500()
        {
            // St = 0.198 * (1 - 19.7 / 500) = 0.198 * 0.9606 = 0.19020
            double st = 500.0.RoshkoStrouhalNumber();
            Assert.AreEqual(0.198 * (1.0 - 19.7 / 500.0), st, 1e-10);
        }

        [TestMethod]
        public void RoshkoStrouhalNumber_HighRe_ApproachesAsymptote()
        {
            // As Re → ∞, St → 0.198
            double st = 1e6.RoshkoStrouhalNumber();
            Assert.AreEqual(0.198, st, 0.005);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void RoshkoStrouhalNumber_NegativeRe_Throws()
        {
            (-100.0).RoshkoStrouhalNumber();
        }

        #endregion

        #region Roshko Number

        [TestMethod]
        public void RoshkoNumber_FromValues()
        {
            // f = 20 Hz, D = 0.05 m, ν = 1e-5 m²/s → Ro = 20 * 0.05² / 1e-5 = 5000
            double ro = 20.0.RoshkoNumber(0.05, 1e-5);
            Assert.AreEqual(5000.0, ro, 1e-6);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void RoshkoNumber_ZeroViscosity_Throws()
        {
            20.0.RoshkoNumber(0.05, 0.0);
        }

        #endregion

        #region Vortex Street Geometry

        [TestMethod]
        public void StableSpacingRatio_MatchesTheoretical()
        {
            // h/a = (1/π) * arccosh(√2) = (1/π) * ln(√2 + 1) ≈ 0.2806
            double expected = Math.Log(Math.Sqrt(2) + 1) / Math.PI;
            Assert.AreEqual(expected, KarmanVortexStreetExtensions.StableSpacingRatio, 1e-12);
            Assert.AreEqual(0.2806, KarmanVortexStreetExtensions.StableSpacingRatio, 0.001);
        }

        [TestMethod]
        public void VortexStreamwiseSpacing_Computed()
        {
            // U = 10 m/s, f = 20 Hz → a = 0.5 m
            double a = 10.0.VortexStreamwiseSpacing(20.0);
            Assert.AreEqual(0.5, a, 1e-10);
        }

        [TestMethod]
        public void VortexLateralSpacing_UsesStabilityRatio()
        {
            double a = 1.0;
            double h = a.VortexLateralSpacing();
            Assert.AreEqual(KarmanVortexStreetExtensions.StableSpacingRatio, h, 1e-12);
        }

        [TestMethod]
        public void VortexStreetWavelength_EqualsStreamwiseSpacing()
        {
            double lambda = 10.0.VortexStreetWavelength(20.0);
            double a = 10.0.VortexStreamwiseSpacing(20.0);
            Assert.AreEqual(a, lambda, 1e-12);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void VortexStreamwiseSpacing_ZeroFrequency_Throws()
        {
            10.0.VortexStreamwiseSpacing(0.0);
        }

        #endregion

        #region Convection Velocity

        [TestMethod]
        public void VortexConvectionVelocity_Default()
        {
            // Default ratio 0.875
            double Uv = 10.0.VortexConvectionVelocity();
            Assert.AreEqual(8.75, Uv, 1e-10);
        }

        [TestMethod]
        public void VortexConvectionVelocity_CustomRatio()
        {
            double Uv = 10.0.VortexConvectionVelocity(0.9);
            Assert.AreEqual(9.0, Uv, 1e-10);
        }

        #endregion

        #region Regime Classification

        [TestMethod]
        public void IsPeriodicSheddingRegime_WithinRange()
        {
            Assert.IsTrue(100.0.IsPeriodicSheddingRegime());
            Assert.IsTrue(1000.0.IsPeriodicSheddingRegime());
            Assert.IsTrue(50000.0.IsPeriodicSheddingRegime());
        }

        [TestMethod]
        public void IsPeriodicSheddingRegime_OutsideRange()
        {
            Assert.IsFalse(10.0.IsPeriodicSheddingRegime());
            Assert.IsFalse(47.0.IsPeriodicSheddingRegime());
            Assert.IsFalse(300000.0.IsPeriodicSheddingRegime());
        }

        [TestMethod]
        public void CylinderWakeRegime_AllRegimes()
        {
            Assert.AreEqual("Creeping", 1.0.CylinderWakeRegime());
            Assert.AreEqual("Steady separation", 20.0.CylinderWakeRegime());
            Assert.AreEqual("Laminar shedding", 100.0.CylinderWakeRegime());
            Assert.AreEqual("Turbulent transition", 5000.0.CylinderWakeRegime());
            Assert.AreEqual("Turbulent wake", 500000.0.CylinderWakeRegime());
        }

        #endregion

        #region Vortex Street Drag Coefficient

        [TestMethod]
        public void VortexStreetDragCoefficient_PositiveResult()
        {
            // Use typical values: a = 0.5 m, h = 0.281 * 0.5, U = 10, Uv = 8.75, D = 0.1
            double a = 0.5;
            double h = a.VortexLateralSpacing();
            double Cd = h.VortexStreetDragCoefficient(a, 10.0, 8.75, 0.1);
            Assert.IsTrue(Cd > 0, "Drag coefficient should be positive.");
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void VortexStreetDragCoefficient_ZeroDiameter_Throws()
        {
            0.14.VortexStreetDragCoefficient(0.5, 10.0, 8.75, 0.0);
        }

        #endregion
    }
}
