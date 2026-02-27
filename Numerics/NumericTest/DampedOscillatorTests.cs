using CSharpNumerics.Physics.Oscillations;
using CSharpNumerics.Statistics.Data;
using System;
using System.Collections.Generic;
using System.Linq;

namespace NumericsTests
{
    [TestClass]
    public class DampedOscillatorTests
    {
        #region Constructor Validation

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Constructor_ZeroMass_Throws()
        {
            new DampedOscillator(0, 10, 1);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Constructor_NegativeStiffness_Throws()
        {
            new DampedOscillator(1, -5, 1);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Constructor_NegativeDamping_Throws()
        {
            new DampedOscillator(1, 10, -1);
        }

        [TestMethod]
        public void Constructor_ZeroDamping_IsAllowed()
        {
            var osc = new DampedOscillator(1, 10, 0);
            Assert.AreEqual(0.0, osc.Gamma, 1e-12);
        }

        #endregion

        #region Damping Regime Detection

        [TestMethod]
        public void Regime_Underdamped_DetectedCorrectly()
        {
            // m=1, k=100 → ω₀=10. c=4 → γ=2 < 10
            var osc = new DampedOscillator(1, 100, 4);
            Assert.AreEqual(DampingRegime.Underdamped, osc.Regime);
        }

        [TestMethod]
        public void Regime_CriticallyDamped_DetectedCorrectly()
        {
            // m=1, k=100 → ω₀=10. c=20 → γ=10 = ω₀
            var osc = new DampedOscillator(1, 100, 20);
            Assert.AreEqual(DampingRegime.CriticallyDamped, osc.Regime);
        }

        [TestMethod]
        public void Regime_Overdamped_DetectedCorrectly()
        {
            // m=1, k=100 → ω₀=10. c=40 → γ=20 > 10
            var osc = new DampedOscillator(1, 100, 40);
            Assert.AreEqual(DampingRegime.Overdamped, osc.Regime);
        }

        #endregion

        #region Physical Properties

        [TestMethod]
        public void Gamma_IsCorrect()
        {
            // c=6, m=3 → γ = 6/(2·3) = 1
            var osc = new DampedOscillator(3, 27, 6);
            Assert.AreEqual(1.0, osc.Gamma, 1e-12);
        }

        [TestMethod]
        public void NaturalFrequency_IsCorrect()
        {
            // m=4, k=100 → ω₀ = √(100/4) = 5
            var osc = new DampedOscillator(4, 100, 2);
            Assert.AreEqual(5.0, osc.NaturalFrequency, 1e-12);
        }

        [TestMethod]
        public void DampedFrequency_Underdamped_IsCorrect()
        {
            // m=1, k=100 → ω₀=10. c=4 → γ=2. ω_d = √(100-4) = √96
            var osc = new DampedOscillator(1, 100, 4);
            double expected = Math.Sqrt(96);
            Assert.AreEqual(expected, osc.DampedFrequency, 1e-10);
        }

        [TestMethod]
        public void DampedFrequency_CriticallyDamped_IsZero()
        {
            var osc = new DampedOscillator(1, 100, 20);
            Assert.AreEqual(0.0, osc.DampedFrequency, 1e-10);
        }

        [TestMethod]
        public void DampedFrequency_Overdamped_IsZero()
        {
            var osc = new DampedOscillator(1, 100, 40);
            Assert.AreEqual(0.0, osc.DampedFrequency, 1e-10);
        }

        [TestMethod]
        public void QualityFactor_IsCorrect()
        {
            // m=1, k=100 → ω₀=10. c=4 → γ=2. Q = 10/(2·2) = 2.5
            var osc = new DampedOscillator(1, 100, 4);
            Assert.AreEqual(2.5, osc.QualityFactor, 1e-10);
        }

        [TestMethod]
        public void QualityFactor_ZeroDamping_ReturnsInfinity()
        {
            var osc = new DampedOscillator(1, 100, 0);
            Assert.AreEqual(double.PositiveInfinity, osc.QualityFactor);
        }

        [TestMethod]
        public void LogarithmicDecrement_Underdamped_IsCorrect()
        {
            // γ=2, ω_d=√96 → T_d = 2π/√96 → δ = 2 · 2π/√96
            var osc = new DampedOscillator(1, 100, 4);
            double wd = Math.Sqrt(96);
            double expected = 2.0 * 2.0 * Math.PI / wd;
            Assert.AreEqual(expected, osc.LogarithmicDecrement, 1e-10);
        }

        [TestMethod]
        public void LogarithmicDecrement_Overdamped_ReturnsNaN()
        {
            var osc = new DampedOscillator(1, 100, 40);
            Assert.IsTrue(double.IsNaN(osc.LogarithmicDecrement));
        }

        [TestMethod]
        public void DampedPeriod_Underdamped_IsCorrect()
        {
            var osc = new DampedOscillator(1, 100, 4);
            double wd = Math.Sqrt(96);
            double expected = 2.0 * Math.PI / wd;
            Assert.AreEqual(expected, osc.DampedPeriod, 1e-10);
        }

        [TestMethod]
        public void DampedPeriod_CriticallyDamped_ReturnsNaN()
        {
            var osc = new DampedOscillator(1, 100, 20);
            Assert.IsTrue(double.IsNaN(osc.DampedPeriod));
        }

        #endregion

        #region Analytic Solution — Underdamped

        [TestMethod]
        public void AnalyticPosition_Underdamped_AtTimeZero_EqualsInitial()
        {
            var osc = new DampedOscillator(1, 100, 4, 2.0, 0.0);
            Assert.AreEqual(2.0, osc.AnalyticPosition(0), 1e-12);
        }

        [TestMethod]
        public void AnalyticVelocity_Underdamped_AtTimeZero_EqualsInitial()
        {
            var osc = new DampedOscillator(1, 100, 4, 2.0, 3.0);
            Assert.AreEqual(3.0, osc.AnalyticVelocity(0), 1e-10);
        }

        [TestMethod]
        public void AnalyticPosition_Underdamped_DecaysToZero()
        {
            var osc = new DampedOscillator(1, 100, 4, 1.0, 0.0);
            // After many time constants, position should be near zero
            double t = 10.0 / osc.Gamma; // ~5 time constants
            Assert.AreEqual(0.0, osc.AnalyticPosition(t), 1e-4);
        }

        #endregion

        #region Analytic Solution — Critically Damped

        [TestMethod]
        public void AnalyticPosition_CriticallyDamped_AtTimeZero_EqualsInitial()
        {
            var osc = new DampedOscillator(1, 100, 20, 2.0, 0.0);
            Assert.AreEqual(2.0, osc.AnalyticPosition(0), 1e-12);
        }

        [TestMethod]
        public void AnalyticVelocity_CriticallyDamped_AtTimeZero_EqualsInitial()
        {
            var osc = new DampedOscillator(1, 100, 20, 2.0, -5.0);
            Assert.AreEqual(-5.0, osc.AnalyticVelocity(0), 1e-10);
        }

        [TestMethod]
        public void AnalyticPosition_CriticallyDamped_NoOscillation()
        {
            // Critically damped: should approach zero without crossing it (for x₀>0, v₀=0)
            var osc = new DampedOscillator(1, 100, 20, 1.0, 0.0);
            double dt = 0.001;
            int crossings = 0;
            double prev = osc.AnalyticPosition(0);
            for (double t = dt; t < 3.0; t += dt)
            {
                double curr = osc.AnalyticPosition(t);
                if (prev * curr < 0) crossings++;
                prev = curr;
            }
            Assert.AreEqual(0, crossings, "Critically damped should not oscillate with x₀>0, v₀=0.");
        }

        [TestMethod]
        public void AnalyticPosition_CriticallyDamped_DecaysToZero()
        {
            var osc = new DampedOscillator(1, 100, 20, 1.0, 0.0);
            Assert.AreEqual(0.0, osc.AnalyticPosition(5.0), 1e-6);
        }

        #endregion

        #region Analytic Solution — Overdamped

        [TestMethod]
        public void AnalyticPosition_Overdamped_AtTimeZero_EqualsInitial()
        {
            var osc = new DampedOscillator(1, 100, 40, 2.0, 0.0);
            Assert.AreEqual(2.0, osc.AnalyticPosition(0), 1e-12);
        }

        [TestMethod]
        public void AnalyticVelocity_Overdamped_AtTimeZero_EqualsInitial()
        {
            var osc = new DampedOscillator(1, 100, 40, 2.0, -3.0);
            Assert.AreEqual(-3.0, osc.AnalyticVelocity(0), 1e-10);
        }

        [TestMethod]
        public void AnalyticPosition_Overdamped_NoOscillation()
        {
            var osc = new DampedOscillator(1, 100, 40, 1.0, 0.0);
            double dt = 0.001;
            int crossings = 0;
            double prev = osc.AnalyticPosition(0);
            for (double t = dt; t < 5.0; t += dt)
            {
                double curr = osc.AnalyticPosition(t);
                if (prev * curr < 0) crossings++;
                prev = curr;
            }
            Assert.AreEqual(0, crossings, "Overdamped should not oscillate with x₀>0, v₀=0.");
        }

        [TestMethod]
        public void AnalyticPosition_Overdamped_DecaysToZero()
        {
            var osc = new DampedOscillator(1, 100, 40, 1.0, 0.0);
            Assert.AreEqual(0.0, osc.AnalyticPosition(5.0), 1e-5);
        }

        #endregion

        #region Numerical vs Analytic

        [TestMethod]
        public void Numerical_MatchesAnalytic_Underdamped()
        {
            // m=1, k=100 → ω₀=10. c=2 → γ=1 (lightly damped, Q=5)
            var osc = new DampedOscillator(1, 100, 2, 1.0, 0.0);
            double dt = 0.0005;
            double tEnd = 5.0; // ~8 damped periods

            for (double t = 0; t < tEnd; t += dt)
            {
                double analyticX = osc.AnalyticPosition(osc.Time);
                Assert.AreEqual(analyticX, osc.Position, 1e-4,
                    $"Position mismatch at t={osc.Time:F4}");
                osc.Step(dt);
            }
        }

        [TestMethod]
        public void Numerical_MatchesAnalytic_CriticallyDamped()
        {
            var osc = new DampedOscillator(1, 100, 20, 1.0, 0.0);
            double dt = 0.0005;
            double tEnd = 2.0;

            for (double t = 0; t < tEnd; t += dt)
            {
                double analyticX = osc.AnalyticPosition(osc.Time);
                Assert.AreEqual(analyticX, osc.Position, 1e-5,
                    $"Position mismatch at t={osc.Time:F4}");
                osc.Step(dt);
            }
        }

        [TestMethod]
        public void Numerical_MatchesAnalytic_Overdamped()
        {
            var osc = new DampedOscillator(1, 100, 40, 1.0, 0.0);
            double dt = 0.0005;
            double tEnd = 2.0;

            for (double t = 0; t < tEnd; t += dt)
            {
                double analyticX = osc.AnalyticPosition(osc.Time);
                Assert.AreEqual(analyticX, osc.Position, 1e-5,
                    $"Position mismatch at t={osc.Time:F4}");
                osc.Step(dt);
            }
        }

        [TestMethod]
        public void NumericalVelocity_MatchesAnalytic_Underdamped()
        {
            var osc = new DampedOscillator(1, 100, 2, 1.0, 0.0);
            double dt = 0.0005;
            double tEnd = 3.0;

            for (double t = 0; t < tEnd; t += dt)
            {
                double analyticV = osc.AnalyticVelocity(osc.Time);
                Assert.AreEqual(analyticV, osc.Velocity, 1e-3,
                    $"Velocity mismatch at t={osc.Time:F4}");
                osc.Step(dt);
            }
        }

        [TestMethod]
        public void Numerical_MatchesAnalytic_WithInitialVelocity()
        {
            var osc = new DampedOscillator(1, 100, 2, 0.0, 5.0);
            double dt = 0.0005;
            double tEnd = 3.0;

            for (double t = 0; t < tEnd; t += dt)
            {
                double analyticX = osc.AnalyticPosition(osc.Time);
                Assert.AreEqual(analyticX, osc.Position, 1e-4,
                    $"Position mismatch at t={osc.Time:F4}");
                osc.Step(dt);
            }
        }

        #endregion

        #region Energy — Monotonically Decreasing

        [TestMethod]
        public void Energy_MonotonicallyDecreasing_Underdamped()
        {
            var osc = new DampedOscillator(1, 100, 2, 1.0, 0.0);
            var energies = osc.EnergyOverTime(5.0, 0.001);

            for (int i = 1; i < energies.Count; i++)
            {
                Assert.IsTrue(energies[i].Value <= energies[i - 1].Value + 1e-10,
                    $"Energy increased at t={energies[i].Index:F4}: " +
                    $"E={energies[i].Value:E6} > E_prev={energies[i - 1].Value:E6}");
            }
        }

        [TestMethod]
        public void Energy_MonotonicallyDecreasing_Overdamped()
        {
            var osc = new DampedOscillator(1, 100, 40, 1.0, 0.0);
            var energies = osc.EnergyOverTime(3.0, 0.001);

            for (int i = 1; i < energies.Count; i++)
            {
                Assert.IsTrue(energies[i].Value <= energies[i - 1].Value + 1e-10,
                    $"Energy increased at t={energies[i].Index:F4}");
            }
        }

        [TestMethod]
        public void Energy_DecaysToZero()
        {
            var osc = new DampedOscillator(1, 100, 4, 1.0, 0.0);
            var energies = osc.EnergyOverTime(10.0, 0.01);
            double finalEnergy = energies.Last().Value;

            Assert.IsTrue(finalEnergy < 1e-6,
                $"Final energy {finalEnergy:E6} should be near zero.");
        }

        [TestMethod]
        public void InitialEnergy_EqualsExpected()
        {
            // E = 0.5*k*x₀² + 0.5*m*v₀² = 0.5*100*1 + 0 = 50
            var osc = new DampedOscillator(1, 100, 4, 1.0, 0.0);
            Assert.AreEqual(50.0, osc.TotalEnergy, 1e-12);
        }

        #endregion

        #region Energy Dissipation

        [TestMethod]
        public void EnergyDissipation_StartsAtZero()
        {
            var osc = new DampedOscillator(1, 100, 4, 1.0, 0.0);
            var dissipation = osc.EnergyDissipation(2.0, 0.01);
            Assert.AreEqual(0.0, dissipation.First().Value, 1e-12);
        }

        [TestMethod]
        public void EnergyDissipation_MonotonicallyIncreasing()
        {
            var osc = new DampedOscillator(1, 100, 4, 1.0, 0.0);
            var dissipation = osc.EnergyDissipation(5.0, 0.001);

            for (int i = 1; i < dissipation.Count; i++)
            {
                Assert.IsTrue(dissipation[i].Value >= dissipation[i - 1].Value - 1e-10,
                    $"Dissipated energy decreased at t={dissipation[i].Index:F4}");
            }
        }

        [TestMethod]
        public void EnergyDissipation_ApproachesInitialEnergy()
        {
            var osc = new DampedOscillator(1, 100, 4, 1.0, 0.0);
            double E0 = osc.TotalEnergy;
            var dissipation = osc.EnergyDissipation(10.0, 0.01);
            double finalDissipated = dissipation.Last().Value;

            Assert.AreEqual(E0, finalDissipated, 1e-3,
                $"Dissipated energy {finalDissipated:F6} should approach initial E₀={E0:F6}");
        }

        #endregion

        #region Envelope

        [TestMethod]
        public void Envelope_AtTimeZero_EqualsAmplitude()
        {
            // x₀=1, v₀=0, γ=2 → envelope A = |x₀| = 1 (since v₀+γx₀ / ωd matters)
            var osc = new DampedOscillator(1, 100, 4, 1.0, 0.0);
            double env0 = osc.Envelope(0);
            // A = √(x₀² + ((v₀+γx₀)/ωd)²) where γ=2, ωd=√96
            double wd = Math.Sqrt(96);
            double expectedA = Math.Sqrt(1.0 + (2.0 / wd) * (2.0 / wd));
            Assert.AreEqual(expectedA, env0, 1e-10);
        }

        [TestMethod]
        public void Envelope_Decays_Exponentially()
        {
            var osc = new DampedOscillator(1, 100, 4, 1.0, 0.0);
            double env1 = osc.Envelope(1.0);
            double env2 = osc.Envelope(2.0);
            // e^(-γ·1) / e^(-γ·2) = e^γ → env1/env2 = e^γ = e^2
            double ratio = env1 / env2;
            double expected = Math.Exp(osc.Gamma);
            Assert.AreEqual(expected, ratio, 1e-10);
        }

        [TestMethod]
        public void EnvelopeCurve_HasCorrectPointCount()
        {
            var osc = new DampedOscillator(1, 100, 4, 1.0, 0.0);
            var curve = osc.EnvelopeCurve(1.0, 0.1);
            // t = 0, 0.1, 0.2, … 1.0 → 11 points
            Assert.AreEqual(11, curve.Count);
        }

        [TestMethod]
        public void Envelope_BoundsAnalyticPosition_Underdamped()
        {
            var osc = new DampedOscillator(1, 100, 4, 1.0, 0.0);
            double dt = 0.01;
            for (double t = 0; t < 5.0; t += dt)
            {
                double pos = Math.Abs(osc.AnalyticPosition(t));
                double env = osc.Envelope(t);
                Assert.IsTrue(pos <= env + 1e-10,
                    $"|x({t:F3})| = {pos:E6} exceeds envelope {env:E6}");
            }
        }

        #endregion

        #region Trajectory & PhasePortrait

        [TestMethod]
        public void Trajectory_StartsAtInitialCondition()
        {
            var osc = new DampedOscillator(1, 100, 4, 2.0, 0.0);
            var traj = osc.Trajectory(1.0, 0.01);

            Assert.AreEqual(0.0, traj.First().Index, 1e-12);
            Assert.AreEqual(2.0, traj.First().Value, 1e-12);
        }

        [TestMethod]
        public void Trajectory_ResetsState()
        {
            var osc = new DampedOscillator(1, 100, 4, 2.0, 3.0);
            osc.Trajectory(1.0, 0.01);

            Assert.AreEqual(0.0, osc.Time, 1e-12);
            Assert.AreEqual(2.0, osc.Position, 1e-12);
            Assert.AreEqual(3.0, osc.Velocity, 1e-12);
        }

        [TestMethod]
        public void PhasePortrait_StartsAtInitialState()
        {
            var osc = new DampedOscillator(1, 100, 4, 2.0, 1.0);
            var phase = osc.PhasePortrait(1.0, 0.01);

            Assert.AreEqual(2.0, phase.First().Index, 1e-12);
            Assert.AreEqual(1.0, phase.First().Value, 1e-12);
        }

        [TestMethod]
        public void PhasePortrait_SpiralsToOrigin_Underdamped()
        {
            var osc = new DampedOscillator(1, 100, 4, 1.0, 0.0);
            var phase = osc.PhasePortrait(10.0, 0.01);

            // Last point should be near (0, 0)
            var last = phase.Last();
            Assert.AreEqual(0.0, last.Index, 0.01, $"Final position {last.Index:E4} not near zero.");
            Assert.AreEqual(0.0, last.Value, 0.05, $"Final velocity {last.Value:E4} not near zero.");
        }

        #endregion

        #region Frequency Spectrum

        [TestMethod]
        public void FrequencySpectrum_PeakNearDampedFrequency()
        {
            // Lightly damped so we get a clear peak
            var osc = new DampedOscillator(1, 100, 1, 1.0, 0.0);
            double expectedFreqHz = osc.DampedFrequency / (2 * Math.PI);

            double dt = 0.005;
            double tEnd = 20.0;
            var spectrum = osc.FrequencySpectrum(tEnd, dt);

            // Find the peak (skip DC)
            var peak = spectrum.Skip(1).Take(spectrum.Count / 2)
                               .OrderByDescending(s => s.Value)
                               .First();

            double freqResolution = 1.0 / tEnd;
            Assert.AreEqual(expectedFreqHz, peak.Index, freqResolution * 3,
                $"FFT peak at {peak.Index:F4} Hz, expected {expectedFreqHz:F4} Hz");
        }

        #endregion

        #region Q-Factor and Logarithmic Decrement — Known Formulas

        [TestMethod]
        public void QualityFactor_MatchesAlternativeFormula()
        {
            // The exact relation: δ = 2πγ/ωd, so π/δ = ωd/(2γ).
            // This equals Q = ω₀/(2γ) only for small damping.
            // Test the exact identity: π/δ = ωd/(2γ).
            var osc = new DampedOscillator(1, 100, 4, 1.0, 0.0);
            double delta = osc.LogarithmicDecrement;
            double piOverDelta = Math.PI / delta;
            double wdOver2Gamma = osc.DampedFrequency / (2.0 * osc.Gamma);

            Assert.AreEqual(wdOver2Gamma, piOverDelta, 1e-10,
                $"ωd/(2γ)={wdOver2Gamma:F6} ≠ π/δ={piOverDelta:F6}");
        }

        [TestMethod]
        public void LogarithmicDecrement_MeasuredFromTrajectory()
        {
            // Measure ratio of successive peaks in the numerical trajectory
            var osc = new DampedOscillator(1, 100, 2, 1.0, 0.0); // γ=1, ω₀=10, lightly damped
            double dt = 0.0001;
            double tEnd = 5.0;

            var traj = osc.Trajectory(tEnd, dt);

            // Find local maxima (peaks)
            var peaks = new List<double>();
            for (int i = 1; i < traj.Count - 1; i++)
            {
                if (traj[i].Value > traj[i - 1].Value && traj[i].Value > traj[i + 1].Value
                    && traj[i].Value > 0.001) // skip noise near zero
                {
                    peaks.Add(traj[i].Value);
                }
            }

            Assert.IsTrue(peaks.Count >= 3, $"Expected at least 3 peaks, found {peaks.Count}");

            // δ = ln(x_n / x_{n+1}) should be approximately constant
            double expectedDelta = osc.LogarithmicDecrement;
            for (int i = 0; i < peaks.Count - 1; i++)
            {
                double measuredDelta = Math.Log(peaks[i] / peaks[i + 1]);
                Assert.AreEqual(expectedDelta, measuredDelta, 0.05,
                    $"Measured δ={measuredDelta:F4} ≠ expected δ={expectedDelta:F4} at peak {i}");
            }
        }

        #endregion

        #region Reset

        [TestMethod]
        public void Reset_RestoresInitialState()
        {
            var osc = new DampedOscillator(1, 100, 4, 2.0, 3.0);
            osc.Step(0.1);
            osc.Step(0.1);

            Assert.AreNotEqual(2.0, osc.Position);

            osc.Reset();

            Assert.AreEqual(0.0, osc.Time, 1e-12);
            Assert.AreEqual(2.0, osc.Position, 1e-12);
            Assert.AreEqual(3.0, osc.Velocity, 1e-12);
        }

        #endregion

        #region ToString

        [TestMethod]
        public void ToString_ContainsKeyParameters()
        {
            var osc = new DampedOscillator(2, 50, 4);
            string s = osc.ToString();
            Assert.IsTrue(s.Contains("DampedOscillator"));
            Assert.IsTrue(s.Contains("kg"));
            Assert.IsTrue(s.Contains("N/m"));
            Assert.IsTrue(s.Contains("regime="));
        }

        #endregion

        #region Edge Cases — Zero Damping Behaves Like SHM

        [TestMethod]
        public void ZeroDamping_BehavesLikeSHM()
        {
            // With c=0, DampedOscillator (RK4) should closely match SimpleHarmonicOscillator (Verlet).
            // Different integrators, so we allow a small tolerance.
            var damped = new DampedOscillator(1, 25, 0, 1.0, 0.0);
            var sho = new SimpleHarmonicOscillator(1, 25, 1.0, 0.0);

            double dt = 0.001;
            for (int i = 0; i < 1000; i++)
            {
                Assert.AreEqual(sho.Position, damped.Position, 1e-5,
                    $"Position mismatch at step {i}");
                sho.Step(dt);
                damped.Step(dt);
            }
        }

        #endregion
    }
}
