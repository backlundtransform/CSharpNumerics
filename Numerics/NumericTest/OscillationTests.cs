using CSharpNumerics.Physics.Oscillations;
using CSharpNumerics.Statistics.Data;
using System;
using System.Collections.Generic;
using System.Linq;

namespace NumericsTests
{
    [TestClass]
    public class OscillationTests
    {
        #region Constructor Validation

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Constructor_ZeroMass_Throws()
        {
            new SimpleHarmonicOscillator(0, 10);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Constructor_NegativeStiffness_Throws()
        {
            new SimpleHarmonicOscillator(1, -5);
        }

        #endregion

        #region Physical Properties

        [TestMethod]
        public void AngularFrequency_IsCorrect()
        {
            // ω₀ = √(k/m) = √(100/4) = 5 rad/s
            var sho = new SimpleHarmonicOscillator(4, 100);
            Assert.AreEqual(5.0, sho.AngularFrequency, 1e-12);
        }

        [TestMethod]
        public void Frequency_IsCorrect()
        {
            var sho = new SimpleHarmonicOscillator(4, 100);
            double expected = 5.0 / (2.0 * Math.PI);
            Assert.AreEqual(expected, sho.Frequency, 1e-12);
        }

        [TestMethod]
        public void Period_IsCorrect()
        {
            var sho = new SimpleHarmonicOscillator(4, 100);
            double expected = 2.0 * Math.PI / 5.0;
            Assert.AreEqual(expected, sho.Period, 1e-12);
        }

        [TestMethod]
        public void Amplitude_PositionOnly()
        {
            // x₀ = 3, v₀ = 0 → A = 3
            var sho = new SimpleHarmonicOscillator(1, 1, 3.0, 0.0);
            Assert.AreEqual(3.0, sho.Amplitude, 1e-12);
        }

        [TestMethod]
        public void Amplitude_WithVelocity()
        {
            // x₀ = 0, v₀ = 5, ω₀ = 1 → A = √(0 + 25) = 5
            var sho = new SimpleHarmonicOscillator(1, 1, 0.0, 5.0);
            Assert.AreEqual(5.0, sho.Amplitude, 1e-12);
        }

        [TestMethod]
        public void Amplitude_Mixed()
        {
            // x₀ = 3, v₀ = 4, ω₀ = 1 → A = √(9 + 16) = 5
            var sho = new SimpleHarmonicOscillator(1, 1, 3.0, 4.0);
            Assert.AreEqual(5.0, sho.Amplitude, 1e-12);
        }

        #endregion

        #region Analytic Solution

        [TestMethod]
        public void AnalyticPosition_AtTimeZero_EqualsInitialPosition()
        {
            var sho = new SimpleHarmonicOscillator(2, 50, 0.5, 0.0);
            Assert.AreEqual(0.5, sho.AnalyticPosition(0), 1e-12);
        }

        [TestMethod]
        public void AnalyticVelocity_AtTimeZero_EqualsInitialVelocity()
        {
            var sho = new SimpleHarmonicOscillator(2, 50, 0.5, 3.0);
            Assert.AreEqual(3.0, sho.AnalyticVelocity(0), 1e-10);
        }

        [TestMethod]
        public void AnalyticSolution_Periodic()
        {
            // After one full period, should return to initial position
            var sho = new SimpleHarmonicOscillator(1, 100, 1.0, 0.0);
            double T = sho.Period;
            Assert.AreEqual(1.0, sho.AnalyticPosition(T), 1e-10);
            Assert.AreEqual(0.0, sho.AnalyticVelocity(T), 1e-10);
        }

        #endregion

        #region Numerical vs Analytic

        [TestMethod]
        public void NumericalSolution_MatchesAnalytic()
        {
            // m = 1, k = 25, ω₀ = 5, x₀ = 1, v₀ = 0
            var sho = new SimpleHarmonicOscillator(1, 25, 1.0, 0.0);
            double dt = 0.001;
            double tEnd = sho.Period * 5; // 5 full periods

            for (double t = 0; t < tEnd; t += dt)
            {
                double analyticX = sho.AnalyticPosition(sho.Time);
                Assert.AreEqual(analyticX, sho.Position, 1e-4,
                    $"Position mismatch at t={sho.Time:F4}");
                sho.Step(dt);
            }
        }

        [TestMethod]
        public void NumericalVelocity_MatchesAnalytic()
        {
            var sho = new SimpleHarmonicOscillator(1, 25, 1.0, 0.0);
            double dt = 0.001;
            double tEnd = sho.Period * 3;

            for (double t = 0; t < tEnd; t += dt)
            {
                double analyticV = sho.AnalyticVelocity(sho.Time);
                Assert.AreEqual(analyticV, sho.Velocity, 1e-3,
                    $"Velocity mismatch at t={sho.Time:F4}");
                sho.Step(dt);
            }
        }

        [TestMethod]
        public void NumericalSolution_WithInitialVelocity_MatchesAnalytic()
        {
            var sho = new SimpleHarmonicOscillator(2, 8, 0.0, 3.0);
            double dt = 0.0005;
            double tEnd = sho.Period * 3;

            for (double t = 0; t < tEnd; t += dt)
            {
                double analyticX = sho.AnalyticPosition(sho.Time);
                Assert.AreEqual(analyticX, sho.Position, 1e-4,
                    $"Position mismatch at t={sho.Time:F4}");
                sho.Step(dt);
            }
        }

        #endregion

        #region Energy Conservation

        [TestMethod]
        public void EnergyIsConserved_OverManyPeriods()
        {
            var sho = new SimpleHarmonicOscillator(1, 25, 1.0, 0.0);
            double dt = 0.001;
            double tEnd = sho.Period * 10; // 10 full periods
            double E0 = sho.TotalEnergy;

            while (sho.Time < tEnd)
            {
                sho.Step(dt);
                double relativeError = Math.Abs(sho.TotalEnergy - E0) / E0;
                Assert.IsTrue(relativeError < 1e-5,
                    $"Energy drift {relativeError:E2} at t={sho.Time:F4}");
            }
        }

        [TestMethod]
        public void EnergyOverTime_ReturnsConsistentValues()
        {
            var sho = new SimpleHarmonicOscillator(1, 25, 1.0, 0.0);
            var energies = sho.EnergyOverTime(sho.Period * 2, 0.01);

            double E0 = energies.First().Value;
            foreach (var point in energies)
            {
                double relativeError = Math.Abs(point.Value - E0) / E0;
                Assert.IsTrue(relativeError < 1e-3,
                    $"Energy drift {relativeError:E2} at t={point.Index:F4}");
            }
        }

        [TestMethod]
        public void InitialEnergy_KineticPlusPotential()
        {
            // E = 0.5*k*x₀² + 0.5*m*v₀² = 0.5*25*1 + 0.5*1*0 = 12.5
            var sho = new SimpleHarmonicOscillator(1, 25, 1.0, 0.0);
            Assert.AreEqual(12.5, sho.TotalEnergy, 1e-12);
            Assert.AreEqual(0.0, sho.KineticEnergy, 1e-12);
            Assert.AreEqual(12.5, sho.PotentialEnergy, 1e-12);
        }

        #endregion

        #region Trajectory & PhasePortrait

        [TestMethod]
        public void Trajectory_ReturnsCorrectNumberOfPoints()
        {
            var sho = new SimpleHarmonicOscillator(1, 25, 1.0, 0.0);
            double dt = 0.01;
            double tEnd = 1.0;
            var traj = sho.Trajectory(tEnd, dt);

            // Expect roughly tEnd/dt + 1 points
            Assert.IsTrue(traj.Count >= (int)(tEnd / dt),
                $"Expected at least {(int)(tEnd / dt)} points, got {traj.Count}");
        }

        [TestMethod]
        public void Trajectory_StartsAtInitialCondition()
        {
            var sho = new SimpleHarmonicOscillator(1, 25, 2.0, 0.0);
            var traj = sho.Trajectory(1.0, 0.01);

            Assert.AreEqual(0.0, traj.First().Index, 1e-12);
            Assert.AreEqual(2.0, traj.First().Value, 1e-12);
        }

        [TestMethod]
        public void Trajectory_ResetsState()
        {
            var sho = new SimpleHarmonicOscillator(1, 25, 2.0, 0.0);
            sho.Trajectory(1.0, 0.01);

            // After Trajectory call, state should be reset
            Assert.AreEqual(0.0, sho.Time, 1e-12);
            Assert.AreEqual(2.0, sho.Position, 1e-12);
            Assert.AreEqual(0.0, sho.Velocity, 1e-12);
        }

        [TestMethod]
        public void PhasePortrait_IsEllipse()
        {
            // For SHM with x₀=A, v₀=0: phase portrait is an ellipse
            // x² / A² + v² / (Aω)² = 1
            var sho = new SimpleHarmonicOscillator(1, 25, 1.0, 0.0);
            double w = sho.AngularFrequency;
            double A = sho.Amplitude;
            var phase = sho.PhasePortrait(sho.Period * 2, 0.001);

            foreach (var point in phase)
            {
                double ellipseValue = (point.Index * point.Index) / (A * A)
                                    + (point.Value * point.Value) / (A * w * A * w);
                Assert.AreEqual(1.0, ellipseValue, 1e-3,
                    $"Point ({point.Index:F4}, {point.Value:F4}) off ellipse: {ellipseValue:F6}");
            }
        }

        #endregion

        #region Frequency Spectrum (FFT)

        [TestMethod]
        public void FrequencySpectrum_PeakAtNaturalFrequency()
        {
            var sho = new SimpleHarmonicOscillator(1, 25, 1.0, 0.0);
            // ω₀ = 5, f₀ ≈ 0.796 Hz
            double expectedFreq = sho.Frequency;

            // Long simulation for good frequency resolution
            double dt = 0.01;
            double tEnd = sho.Period * 20;
            var spectrum = sho.FrequencySpectrum(tEnd, dt);

            // Find the peak (skip DC at index 0)
            var peak = spectrum.Skip(1).Take(spectrum.Count / 2)
                               .OrderByDescending(s => s.Value)
                               .First();

            double freqResolution = 1.0 / tEnd;
            Assert.AreEqual(expectedFreq, peak.Index, freqResolution * 2,
                $"FFT peak at {peak.Index:F4} Hz, expected {expectedFreq:F4} Hz");
        }

        #endregion

        #region Reset

        [TestMethod]
        public void Reset_RestoresInitialState()
        {
            var sho = new SimpleHarmonicOscillator(1, 25, 2.0, 3.0);
            sho.Step(0.1);
            sho.Step(0.1);
            sho.Step(0.1);

            Assert.AreNotEqual(2.0, sho.Position);

            sho.Reset();

            Assert.AreEqual(0.0, sho.Time, 1e-12);
            Assert.AreEqual(2.0, sho.Position, 1e-12);
            Assert.AreEqual(3.0, sho.Velocity, 1e-12);
        }

        #endregion

        #region ToString

        [TestMethod]
        public void ToString_ContainsKeyParameters()
        {
            var sho = new SimpleHarmonicOscillator(2, 50);
            string s = sho.ToString();
            Assert.IsTrue(s.Contains("SHM"));
            Assert.IsTrue(s.Contains("kg"));
            Assert.IsTrue(s.Contains("N/m"));
        }

        #endregion
    }
}
