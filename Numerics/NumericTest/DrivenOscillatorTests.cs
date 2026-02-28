using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Oscillations;
using CSharpNumerics.Statistics.Data;
using System;
using System.Linq;

namespace NumericsTests
{
    [TestClass]
    public class DrivenOscillatorTests
    {
        #region Constructor Validation

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Constructor_ZeroMass_Throws()
        {
            new DrivenOscillator(0, 10, 1, 5, 3);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Constructor_NegativeStiffness_Throws()
        {
            new DrivenOscillator(1, -5, 1, 5, 3);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Constructor_NegativeDamping_Throws()
        {
            new DrivenOscillator(1, 10, -1, 5, 3);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Constructor_NegativeDriveAmplitude_Throws()
        {
            new DrivenOscillator(1, 10, 1, -5, 3);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Constructor_NegativeDriveFrequency_Throws()
        {
            new DrivenOscillator(1, 10, 1, 5, -3);
        }

        [TestMethod]
        public void Constructor_ValidParameters_SetsInitialState()
        {
            var osc = new DrivenOscillator(2, 50, 4, 10, 7, 1.5, 0.5);
            Assert.AreEqual(1.5, osc.Position, 1e-12);
            Assert.AreEqual(0.5, osc.Velocity, 1e-12);
            Assert.AreEqual(0.0, osc.Time, 1e-12);
        }

        #endregion

        #region Physical Properties

        [TestMethod]
        public void NaturalFrequency_IsCorrect()
        {
            // m=4, k=100 → ω₀ = √(100/4) = 5
            var osc = new DrivenOscillator(4, 100, 2, 10, 3);
            Assert.AreEqual(5.0, osc.NaturalFrequency, 1e-12);
        }

        [TestMethod]
        public void Gamma_IsCorrect()
        {
            // c=6, m=3 → γ = 6/(2·3) = 1
            var osc = new DrivenOscillator(3, 27, 6, 10, 3);
            Assert.AreEqual(1.0, osc.Gamma, 1e-12);
        }

        [TestMethod]
        public void DampedFrequency_IsCorrect()
        {
            // m=1, k=100 → ω₀=10. c=4 → γ=2. ω_d = √(100-4) = √96
            var osc = new DrivenOscillator(1, 100, 4, 10, 5);
            Assert.AreEqual(Math.Sqrt(96), osc.DampedFrequency, 1e-10);
        }

        [TestMethod]
        public void QualityFactor_IsCorrect()
        {
            // m=1, k=100 → ω₀=10. c=4 → γ=2. Q = 10/(2·2) = 2.5
            var osc = new DrivenOscillator(1, 100, 4, 10, 5);
            Assert.AreEqual(2.5, osc.QualityFactor, 1e-10);
        }

        [TestMethod]
        public void ResonanceFrequency_IsCorrect()
        {
            // ω_r = √(ω₀² − 2γ²) = √(100 − 8) = √92
            var osc = new DrivenOscillator(1, 100, 4, 10, 5);
            Assert.AreEqual(Math.Sqrt(92), osc.ResonanceFrequency, 1e-10);
        }

        [TestMethod]
        public void ResonanceFrequency_HeavilyDamped_ReturnsZero()
        {
            // m=1, k=1 → ω₀=1. c=4 → γ=2. ω₀²=1, 2γ²=8 → disc<0
            var osc = new DrivenOscillator(1, 1, 4, 1, 0.5);
            Assert.AreEqual(0.0, osc.ResonanceFrequency, 1e-10);
        }

        [TestMethod]
        public void Bandwidth_IsCorrect()
        {
            // Δω = 2γ = 2·2 = 4
            var osc = new DrivenOscillator(1, 100, 4, 10, 5);
            Assert.AreEqual(4.0, osc.Bandwidth, 1e-10);
        }

        [TestMethod]
        public void Regime_DetectedCorrectly()
        {
            var under = new DrivenOscillator(1, 100, 4, 10, 5);
            Assert.AreEqual(DampingRegime.Underdamped, under.Regime);

            var crit = new DrivenOscillator(1, 100, 20, 10, 5);
            Assert.AreEqual(DampingRegime.CriticallyDamped, crit.Regime);

            var over = new DrivenOscillator(1, 100, 40, 10, 5);
            Assert.AreEqual(DampingRegime.Overdamped, over.Regime);
        }

        #endregion

        #region Steady-State Amplitude Formula

        [TestMethod]
        public void SteadyStateAmplitude_AtResonance_IsMaximum()
        {
            // m=1, k=100, c=2 → ω₀=10, γ=1
            // ω_r = √(100−2) = √98
            var osc = new DrivenOscillator(1, 100, 2, 10, Math.Sqrt(98));

            double wr = osc.ResonanceFrequency;
            double Amax = osc.SteadyStateAmplitude(wr);

            // Check that nearby frequencies give lower amplitudes
            double Abelow = osc.SteadyStateAmplitude(wr - 1);
            double Aabove = osc.SteadyStateAmplitude(wr + 1);

            Assert.IsTrue(Amax > Abelow, $"A(ωr)={Amax} should be > A(ωr-1)={Abelow}");
            Assert.IsTrue(Amax > Aabove, $"A(ωr)={Amax} should be > A(ωr+1)={Aabove}");
        }

        [TestMethod]
        public void SteadyStateAmplitude_AtZeroFrequency_EqualsStaticDeflection()
        {
            // At ω_d=0: A = (F₀/m) / ω₀² = F₀/k
            var osc = new DrivenOscillator(1, 100, 2, 10, 5);
            double expected = 10.0 / 100.0; // F₀/k = 0.1
            Assert.AreEqual(expected, osc.SteadyStateAmplitude(0), 1e-12);
        }

        [TestMethod]
        public void SteadyStateAmplitude_KnownValues()
        {
            // m=1, k=100, c=4 → ω₀=10, γ=2
            // A(ω_d=5) = (F₀/m) / √((ω₀²-ω_d²)² + (2γω_d)²)
            //          = 10 / √((100-25)² + (2·2·5)²) = 10 / √(5625 + 400) = 10 / √6025
            var osc = new DrivenOscillator(1, 100, 4, 10, 5);
            double expected = 10.0 / Math.Sqrt(6025);
            Assert.AreEqual(expected, osc.SteadyStateAmplitude(5), 1e-12);
        }

        [TestMethod]
        public void SteadyStateAmplitude_HighFrequency_ApproachesZero()
        {
            var osc = new DrivenOscillator(1, 100, 4, 10, 5);
            double A = osc.SteadyStateAmplitude(1000);
            Assert.IsTrue(A < 1e-4, $"Amplitude at ω=1000 should be tiny, got {A}");
        }

        #endregion

        #region Steady-State Phase

        [TestMethod]
        public void SteadyStatePhase_AtZeroFrequency_IsZero()
        {
            // At ω_d=0: φ = -atan2(0, ω₀²) = 0
            var osc = new DrivenOscillator(1, 100, 4, 10, 5);
            Assert.AreEqual(0.0, osc.SteadyStatePhase(0), 1e-12);
        }

        [TestMethod]
        public void SteadyStatePhase_AtNaturalFrequency_IsMinusPiOver2()
        {
            // At ω_d=ω₀: φ = -atan2(2γω₀, 0) = -π/2
            var osc = new DrivenOscillator(1, 100, 4, 10, 10);
            Assert.AreEqual(-Math.PI / 2.0, osc.SteadyStatePhase(10), 1e-10);
        }

        [TestMethod]
        public void SteadyStatePhase_HighFrequency_ApproachesMinusPi()
        {
            var osc = new DrivenOscillator(1, 100, 4, 10, 5);
            double phase = osc.SteadyStatePhase(1000);
            // Should approach -π from above
            Assert.AreEqual(-Math.PI, phase, 0.05);
        }

        #endregion

        #region Resonance & Phase Curves

        [TestMethod]
        public void ResonanceCurve_PeaksNearResonanceFrequency()
        {
            // m=1, k=100, c=2 → Q=5, ω_r=√98≈9.90
            var osc = new DrivenOscillator(1, 100, 2, 10, 5);
            var curve = osc.ResonanceCurve(0.1, 20, 500);

            // Find the peak
            var peak = curve.OrderByDescending(s => s.Value).First();
            double wr = osc.ResonanceFrequency;

            Assert.AreEqual(wr, peak.Index, 0.1,
                $"Resonance peak at ω={peak.Index:F3}, expected ω_r={wr:F3}");
        }

        [TestMethod]
        public void ResonanceCurve_CorrectStepCount()
        {
            var osc = new DrivenOscillator(1, 100, 2, 10, 5);
            var curve = osc.ResonanceCurve(0, 20, 100);
            Assert.AreEqual(100, curve.Count);
        }

        [TestMethod]
        public void PhaseResponse_CorrectStepCount()
        {
            var osc = new DrivenOscillator(1, 100, 2, 10, 5);
            var curve = osc.PhaseResponse(0.1, 20, 100);
            Assert.AreEqual(100, curve.Count);
        }

        [TestMethod]
        public void PhaseResponse_MonotonicallyDecreasing()
        {
            var osc = new DrivenOscillator(1, 100, 4, 10, 5);
            var curve = osc.PhaseResponse(0.1, 30, 200);

            for (int i = 1; i < curve.Count; i++)
            {
                Assert.IsTrue(curve[i].Value <= curve[i - 1].Value + 1e-10,
                    $"Phase increased at ω={curve[i].Index:F3}: " +
                    $"{curve[i].Value:F4} > {curve[i - 1].Value:F4}");
            }
        }

        #endregion

        #region Transfer Function

        [TestMethod]
        public void TransferFunction_AtImaginaryAxis_MatchesSteadyState()
        {
            // H(iω) evaluated on the imaginary axis should have magnitude = A(ω) / (F₀/m)
            var osc = new DrivenOscillator(1, 100, 4, 10, 5);
            double wd = 7.0;

            ComplexNumber H = osc.FrequencyResponse(wd);
            double magnitude = H.GetMagnitude();

            // |H(iω)| should equal A(ω) / (F₀/m)
            double expectedMag = osc.SteadyStateAmplitude(wd) / (10.0 / 1.0);
            Assert.AreEqual(expectedMag, magnitude, 1e-10);
        }

        [TestMethod]
        public void TransferFunction_PhaseMatchesSteadyStatePhase()
        {
            var osc = new DrivenOscillator(1, 100, 4, 10, 5);
            double wd = 7.0;

            ComplexNumber H = osc.FrequencyResponse(wd);
            double phaseH = H.GetArgument();
            double phaseSS = osc.SteadyStatePhase(wd);

            Assert.AreEqual(phaseSS, phaseH, 1e-10);
        }

        [TestMethod]
        public void TransferFunction_AtDCIsInverseOfStiffness()
        {
            // H(0) = 1/ω₀² = m/k
            var osc = new DrivenOscillator(2, 50, 4, 10, 5);
            ComplexNumber H = osc.TransferFunction(new ComplexNumber(0, 0));
            Assert.AreEqual(1.0 / 25.0, H.realPart, 1e-10); // 1/ω₀² = 1/(50/2) = 1/25
            Assert.AreEqual(0.0, H.imaginaryPart, 1e-10);
        }

        #endregion

        #region Numerical vs Steady-State Amplitude

        [TestMethod]
        public void NumericalAmplitude_MatchesSteadyState_AtResonance()
        {
            // m=1, k=100, c=4 → ω₀=10, γ=2
            // Drive at resonance ω_r = √(100−8) = √92
            double wr = Math.Sqrt(92);
            var osc = new DrivenOscillator(1, 100, 4, 10, wr);

            double theoreticalA = osc.SteadyStateAmplitudeAtDrive;
            double measuredA = osc.MeasuredSteadyStateAmplitude();

            Assert.AreEqual(theoreticalA, measuredA, theoreticalA * 0.02,
                $"Measured={measuredA:F6}, Expected={theoreticalA:F6}");
        }

        [TestMethod]
        public void NumericalAmplitude_MatchesSteadyState_OffResonance()
        {
            // Drive at ω_d = 5 (well below resonance)
            var osc = new DrivenOscillator(1, 100, 4, 10, 5);

            double theoreticalA = osc.SteadyStateAmplitudeAtDrive;
            double measuredA = osc.MeasuredSteadyStateAmplitude();

            Assert.AreEqual(theoreticalA, measuredA, theoreticalA * 0.02,
                $"Measured={measuredA:F6}, Expected={theoreticalA:F6}");
        }

        [TestMethod]
        public void NumericalAmplitude_MatchesSteadyState_AboveResonance()
        {
            // Drive at ω_d = 15 (above resonance)
            var osc = new DrivenOscillator(1, 100, 4, 10, 15);

            double theoreticalA = osc.SteadyStateAmplitudeAtDrive;
            double measuredA = osc.MeasuredSteadyStateAmplitude();

            Assert.AreEqual(theoreticalA, measuredA, theoreticalA * 0.02,
                $"Measured={measuredA:F6}, Expected={theoreticalA:F6}");
        }

        #endregion

        #region Without Driving — Behaves Like DampedOscillator

        [TestMethod]
        public void ZeroDrive_BehavesLikeDamped()
        {
            var driven = new DrivenOscillator(1, 100, 4, 0, 0, 1.0, 0.0);
            var damped = new DampedOscillator(1, 100, 4, 1.0, 0.0);

            double dt = 0.001;
            for (int i = 0; i < 2000; i++)
            {
                Assert.AreEqual(damped.Position, driven.Position, 1e-10,
                    $"Position mismatch at step {i}");
                Assert.AreEqual(damped.Velocity, driven.Velocity, 1e-10,
                    $"Velocity mismatch at step {i}");
                damped.Step(dt);
                driven.Step(dt);
            }
        }

        #endregion

        #region Trajectory & Reset

        [TestMethod]
        public void Trajectory_StartsAtInitialCondition()
        {
            var osc = new DrivenOscillator(1, 100, 4, 10, 5, 2.0, 0.0);
            var traj = osc.Trajectory(1.0, 0.01);

            Assert.AreEqual(0.0, traj.First().Index, 1e-12);
            Assert.AreEqual(2.0, traj.First().Value, 1e-12);
        }

        [TestMethod]
        public void Trajectory_ResetsState()
        {
            var osc = new DrivenOscillator(1, 100, 4, 10, 5, 2.0, 3.0);
            osc.Trajectory(1.0, 0.01);

            Assert.AreEqual(0.0, osc.Time, 1e-12);
            Assert.AreEqual(2.0, osc.Position, 1e-12);
            Assert.AreEqual(3.0, osc.Velocity, 1e-12);
        }

        [TestMethod]
        public void Reset_RestoresInitialState()
        {
            var osc = new DrivenOscillator(1, 100, 4, 10, 5, 2.0, 3.0);
            osc.Step(0.1);
            osc.Step(0.1);

            Assert.AreNotEqual(2.0, osc.Position);

            osc.Reset();

            Assert.AreEqual(0.0, osc.Time, 1e-12);
            Assert.AreEqual(2.0, osc.Position, 1e-12);
            Assert.AreEqual(3.0, osc.Velocity, 1e-12);
        }

        #endregion

        #region Frequency Spectrum

        [TestMethod]
        public void FrequencySpectrum_PeakAtDriveFrequency()
        {
            // m=1, k=100, c=4, F₀=10, ω_d=7 → f_d = 7/(2π) ≈ 1.114 Hz
            var osc = new DrivenOscillator(1, 100, 4, 10, 7, 0.0, 0.0);
            double expectedFreqHz = 7.0 / (2.0 * Math.PI);

            // FrequencySpectrum runs long enough for steady state
            var spectrum = osc.FrequencySpectrum(50.0, 0.005);

            // Find peak (skip DC)
            var peak = spectrum.Skip(1).Take(spectrum.Count / 2)
                               .OrderByDescending(s => s.Value)
                               .First();

            // Frequency resolution = fs / N
            double freqResolution = (1.0 / 0.005) / spectrum.Count;
            Assert.AreEqual(expectedFreqHz, peak.Index, freqResolution * 3,
                $"FFT peak at {peak.Index:F4} Hz, expected {expectedFreqHz:F4} Hz");
        }

        #endregion

        #region Energy & Power

        [TestMethod]
        public void Energy_IncreasesInitially_WhenDrivenFromRest()
        {
            // Start at rest, energy should increase as work is done
            var osc = new DrivenOscillator(1, 100, 4, 10, 10, 0.0, 0.0);
            var energies = osc.EnergyOverTime(1.0, 0.001);

            // After some time, energy should be above zero
            double laterEnergy = energies.Skip(energies.Count / 2).Average(e => e.Value);
            Assert.IsTrue(laterEnergy > 0.01, $"Expected energy to increase, got {laterEnergy}");
        }

        [TestMethod]
        public void PowerInput_HasCorrectPointCount()
        {
            var osc = new DrivenOscillator(1, 100, 4, 10, 5);
            var power = osc.PowerInput(1.0, 0.01);
            Assert.IsTrue(power.Count >= 100);
        }

        [TestMethod]
        public void AveragePower_InSteadyState_IsPositive()
        {
            // In steady state, net average power input = power dissipated > 0
            var osc = new DrivenOscillator(1, 100, 4, 10, 10, 0.0, 0.0);
            double g = osc.Gamma;
            double transient = 5.0 / g;

            osc.Reset();
            // Skip transient
            double dt = 0.001;
            while (osc.Time < transient) osc.Step(dt);

            // Measure average power over several drive periods
            double totalPower = 0;
            int steps = 0;
            double measureDuration = 10.0 * 2.0 * Math.PI / osc.DriveFrequency;
            double tStart = osc.Time;
            while (osc.Time - tStart < measureDuration)
            {
                double force = osc.DriveAmplitude * Math.Cos(osc.DriveFrequency * osc.Time);
                totalPower += force * osc.Velocity;
                osc.Step(dt);
                steps++;
            }

            double avgPower = totalPower / steps;
            Assert.IsTrue(avgPower > 0, $"Average power should be positive, got {avgPower}");
        }

        #endregion

        #region SteadyStateTrajectory

        [TestMethod]
        public void SteadyStateTrajectory_AmplitudeStabilises()
        {
            var osc = new DrivenOscillator(1, 100, 4, 10, 10, 0.0, 0.0);
            var steadyData = osc.SteadyStateTrajectory(5.0, 0.001);

            // The amplitude of the first and last half should be similar
            int halfCount = steadyData.Count / 2;
            double ampFirst = steadyData.Take(halfCount).Max(s => Math.Abs(s.Value));
            double ampLast = steadyData.Skip(halfCount).Max(s => Math.Abs(s.Value));

            double relativeDiff = Math.Abs(ampFirst - ampLast) / ampLast;
            Assert.IsTrue(relativeDiff < 0.05,
                $"Amplitude not stable: first_half={ampFirst:F4}, last_half={ampLast:F4}, diff={relativeDiff:P1}");
        }

        [TestMethod]
        public void SteadyStateTrajectory_ResetsState()
        {
            var osc = new DrivenOscillator(1, 100, 4, 10, 5, 2.0, 3.0);
            osc.SteadyStateTrajectory(1.0, 0.01);

            Assert.AreEqual(0.0, osc.Time, 1e-12);
            Assert.AreEqual(2.0, osc.Position, 1e-12);
        }

        #endregion

        #region Higher Q Gives Sharper Resonance

        [TestMethod]
        public void HigherQ_GivesSharperResonancePeak()
        {
            // Low damping → high Q → sharp peak
            var highQ = new DrivenOscillator(1, 100, 1, 10, 5); // Q=5
            var lowQ = new DrivenOscillator(1, 100, 8, 10, 5);  // Q=0.625

            double peakHigh = highQ.ResonanceCurve(0.1, 20, 500)
                .OrderByDescending(s => s.Value).First().Value;
            double peakLow = lowQ.ResonanceCurve(0.1, 20, 500)
                .OrderByDescending(s => s.Value).First().Value;

            Assert.IsTrue(peakHigh > peakLow,
                $"Higher Q peak ({peakHigh:F4}) should be greater than lower Q peak ({peakLow:F4})");
        }

        #endregion

        #region ToString

        [TestMethod]
        public void ToString_ContainsKeyParameters()
        {
            var osc = new DrivenOscillator(2, 50, 4, 10, 7);
            string s = osc.ToString();
            Assert.IsTrue(s.Contains("DrivenOscillator"));
            Assert.IsTrue(s.Contains("kg"));
            Assert.IsTrue(s.Contains("N/m"));
            Assert.IsTrue(s.Contains("F₀="));
            Assert.IsTrue(s.Contains("ω_d="));
        }

        #endregion
    }
}
