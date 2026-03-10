using CSharpNumerics.Numerics.FiniteDifference;
using CSharpNumerics.Numerics.FiniteDifference.TimeStepping;
using CSharpNumerics.Physics.Waves;
using CSharpNumerics.Numerics.SignalProcessing;
using CSharpNumerics.Statistics.Data;
using System;
using System.Collections.Generic;
using System.Linq;

namespace NumericsTests
{
    [TestClass]
    public class WaveTests
    {
        // ─────────────────────────── WaveEquation1D ───────────────────────────

        #region WaveEquation1D – Constructor Validation

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Wave1D_TooFewPoints_Throws()
        {
            new WaveEquation1D(2, 0.1, 1.0);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Wave1D_ZeroDx_Throws()
        {
            new WaveEquation1D(10, 0.0, 1.0);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Wave1D_NegativeSpeed_Throws()
        {
            new WaveEquation1D(10, 0.1, -1.0);
        }

        #endregion

        #region WaveEquation1D – Standing Wave Frequencies

        [TestMethod]
        public void Wave1D_StandingWaveFrequency_Mode1()
        {
            // f₁ = c / (2L) for fixed/fixed string
            double c = 340.0;
            int n = 101;
            double dx = 0.01;
            var wave = new WaveEquation1D(n, dx, c, BoundaryType.Fixed);
            double L = wave.Length; // 100 * 0.01 = 1.0

            double expected = c / (2.0 * L);
            Assert.AreEqual(expected, wave.StandingWaveFrequency(1), 1e-10);
        }

        [TestMethod]
        public void Wave1D_StandingWaveFrequency_Mode3()
        {
            double c = 100.0;
            int n = 51;
            double dx = 0.02; // L = 1.0
            var wave = new WaveEquation1D(n, dx, c);

            double expected = 3.0 * c / (2.0 * wave.Length);
            Assert.AreEqual(expected, wave.StandingWaveFrequency(3), 1e-10);
        }

        #endregion

        #region WaveEquation1D – Energy Conservation

        [TestMethod]
        public void Wave1D_EnergyConserved_OverManySteps()
        {
            // Undamped wave with Verlet → symplectic → energy should stay nearly constant
            int n = 101;
            double dx = 0.01;
            double c = 1.0;
            double dt = 0.005; // CFL = c*dt/dx = 0.5 < 1

            var wave = new WaveEquation1D(n, dx, c, BoundaryType.Fixed);
            double L = wave.Length;

            // Pluck the centre
            wave.SetInitialCondition(
                u0: x => Math.Sin(Math.PI * x / L),
                v0: null);

            double e0 = wave.TotalEnergy;
            Assert.IsTrue(e0 > 0, "Initial energy must be positive.");

            // Run 200 steps
            for (int i = 0; i < 200; i++)
                wave.Step(dt);

            double eFinal = wave.TotalEnergy;
            double relError = Math.Abs(eFinal - e0) / e0;
            Assert.IsTrue(relError < 0.02, $"Energy drifted by {relError:P2}");
        }

        #endregion

        #region WaveEquation1D – Snapshot

        [TestMethod]
        public void Wave1D_Snapshot_ReturnsCorrectCount()
        {
            int n = 50;
            var wave = new WaveEquation1D(n, 0.1, 1.0);
            wave.SetInitialCondition(x => Math.Sin(x));
            var snap = wave.Snapshot();
            Assert.AreEqual(n, snap.Count);
        }

        [TestMethod]
        public void Wave1D_Snapshot_XValuesMatchGrid()
        {
            int n = 20;
            double dx = 0.25;
            var wave = new WaveEquation1D(n, dx, 1.0);
            wave.SetInitialCondition(x => x);
            var snap = wave.Snapshot();

            for (int i = 0; i < n; i++)
                Assert.AreEqual(i * dx, snap[i].Index, 1e-12);
        }

        #endregion

        #region WaveEquation1D – Reset

        [TestMethod]
        public void Wave1D_Reset_RestoresInitialState()
        {
            int n = 51;
            double dx = 0.02;
            var wave = new WaveEquation1D(n, dx, 1.0);
            wave.SetInitialCondition(x => Math.Sin(Math.PI * x));

            var snapBefore = wave.Snapshot().Select(s => s.Value).ToArray();
            wave.Simulate(0.5, 0.005);

            wave.Reset();
            var snapAfter = wave.Snapshot().Select(s => s.Value).ToArray();

            for (int i = 0; i < n; i++)
                Assert.AreEqual(snapBefore[i], snapAfter[i], 1e-14);
            Assert.AreEqual(0.0, wave.Time, 1e-14);
        }

        #endregion

        #region WaveEquation1D – SpaceTimeField

        [TestMethod]
        public void Wave1D_SpaceTimeField_HasCorrectDimensions()
        {
            int n = 31;
            double dx = 0.1;
            double c = 1.0;
            double dt = 0.05;
            double tEnd = 0.5;

            var wave = new WaveEquation1D(n, dx, c);
            wave.SetInitialCondition(x => Math.Exp(-100 * (x - 1.5) * (x - 1.5)));

            var field = wave.SpaceTimeField(tEnd, dt);
            Assert.AreEqual(n, field.rowLength);
            Assert.IsTrue(field.columnLength >= 2, "Should have at least initial + 1 step.");
        }

        #endregion

        #region WaveEquation1D – CFL

        [TestMethod]
        public void Wave1D_CFL_CorrectValue()
        {
            double c = 2.0, dx = 0.1, dt = 0.03;
            var wave = new WaveEquation1D(10, dx, c);
            Assert.AreEqual(c * dt / dx, wave.CFL(dt), 1e-14);
        }

        #endregion

        // ─────────────────────────── WaveEquation2D ───────────────────────────

        #region WaveEquation2D – Constructor

        [TestMethod]
        [ExpectedException(typeof(ArgumentNullException))]
        public void Wave2D_NullGrid_Throws()
        {
            new WaveEquation2D(null, 1.0);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Wave2D_ZeroSpeed_Throws()
        {
            var grid = new Grid2D(10, 10, 0.1);
            new WaveEquation2D(grid, 0.0);
        }

        #endregion

        #region WaveEquation2D – Symmetric Initial Condition Stays Symmetric

        [TestMethod]
        public void Wave2D_SymmetricPulse_StaysSymmetric()
        {
            int nx = 21, ny = 21;
            double d = 0.05;
            var grid = new Grid2D(nx, ny, d);
            double c = 1.0;
            var wave = new WaveEquation2D(grid, c, BoundaryType.Fixed);

            double cx = (nx - 1) * d / 2.0;
            double cy = (ny - 1) * d / 2.0;

            wave.SetInitialCondition(
                u0: (x, y) => Math.Exp(-200.0 * ((x - cx) * (x - cx) + (y - cy) * (y - cy))));

            double dt = 0.02;
            for (int step = 0; step < 10; step++)
                wave.Step(dt);

            var arr = wave.SnapshotArray();

            // Check symmetry along both axes through centre
            int mid = nx / 2;
            for (int i = 1; i <= mid; i++)
            {
                // Left-right symmetry
                Assert.AreEqual(arr[mid - i, mid], arr[mid + i, mid], 1e-8,
                    $"X-symmetry broken at offset {i}");
                // Top-bottom symmetry
                Assert.AreEqual(arr[mid, mid - i], arr[mid, mid + i], 1e-8,
                    $"Y-symmetry broken at offset {i}");
            }
        }

        #endregion

        #region WaveEquation2D – Energy Conservation

        [TestMethod]
        public void Wave2D_EnergyConserved()
        {
            // Use a standing-wave mode sin(πx/Lx)·sin(πy/Ly) which satisfies
            // fixed boundary conditions exactly → minimal energy leakage.
            int nx = 21, ny = 21;
            double d = 0.05;
            var grid = new Grid2D(nx, ny, d);
            var wave = new WaveEquation2D(grid, 1.0, BoundaryType.Fixed);

            double Lx = (nx - 1) * d;
            double Ly = (ny - 1) * d;
            wave.SetInitialCondition(
                u0: (x, y) => Math.Sin(Math.PI * x / Lx) * Math.Sin(Math.PI * y / Ly));

            double e0 = wave.TotalEnergy;
            Assert.IsTrue(e0 > 0);

            wave.Simulate(0.1, 0.005);
            double eFinal = wave.TotalEnergy;

            double relError = Math.Abs(eFinal - e0) / e0;
            Assert.IsTrue(relError < 0.10, $"2D energy drifted by {relError:P2}");
        }

        #endregion

        // ────────────────────────── WaveSuperposition ─────────────────────────

        #region WaveSuperposition – Single Harmonic

        [TestMethod]
        public void Superposition_SingleHarmonic_MatchesSin()
        {
            var ws = new WaveSuperposition();
            double A = 2.0, w = 3.0, k = 1.5;
            ws.AddHarmonic(A, w, k);

            double x = 1.0, t = 0.5;
            double expected = A * Math.Sin(w * t - k * x);
            Assert.AreEqual(expected, ws.Evaluate(x, t), 1e-12);
        }

        #endregion

        #region WaveSuperposition – Standing Wave From Two Counter-Propagating

        [TestMethod]
        public void Superposition_StandingWave_NodesAtCorrectPositions()
        {
            // Two waves travelling in opposite directions with same A, ω, |k|
            // u = A sin(ωt − kx) + A sin(ωt + kx) = 2A cos(kx) sin(ωt)
            // Nodes where cos(kx) = 0 → kx = π/2, 3π/2, …
            var ws = new WaveSuperposition();
            double A = 1.0, w = 2.0, k = Math.PI; // k = π → nodes at x = 0.5, 1.5, …

            ws.AddHarmonic(A, w, k, 0);
            ws.AddHarmonic(A, w, -k, 0);

            double t = 0.25; // arbitrary non-zero time
            // At node x = 0.5: cos(π·0.5) = cos(π/2) = 0 → u ≈ 0
            Assert.AreEqual(0.0, ws.Evaluate(0.5, t), 1e-10, "Node at x=0.5 failed");
            Assert.AreEqual(0.0, ws.Evaluate(1.5, t), 1e-10, "Node at x=1.5 failed");
        }

        #endregion

        #region WaveSuperposition – Beat Frequency

        [TestMethod]
        public void Superposition_BeatFrequency_IsCorrect()
        {
            var ws = new WaveSuperposition();
            double w1 = 100.0, w2 = 106.0;
            ws.AddHarmonic(1, w1, 1);
            ws.AddHarmonic(1, w2, 1);

            double expected = Math.Abs(w1 - w2) / (2 * Math.PI);
            Assert.AreEqual(expected, ws.BeatFrequency(), 1e-12);
        }

        [TestMethod]
        [ExpectedException(typeof(InvalidOperationException))]
        public void Superposition_BeatFrequency_SingleHarmonic_Throws()
        {
            var ws = new WaveSuperposition();
            ws.AddHarmonic(1, 10, 1);
            ws.BeatFrequency();
        }

        #endregion

        #region WaveSuperposition – EvaluateRange

        [TestMethod]
        public void Superposition_EvaluateRange_MatchesPointwise()
        {
            var ws = new WaveSuperposition();
            ws.AddHarmonic(1, 5, 2);
            ws.AddHarmonic(0.5, 10, 4);

            double t = 0.3;
            var xVals = Enumerable.Range(0, 20).Select(i => i * 0.1).ToArray();
            var series = ws.EvaluateRange(xVals, t);

            for (int i = 0; i < xVals.Length; i++)
                Assert.AreEqual(ws.Evaluate(xVals[i], t), series[i].Value, 1e-14);
        }

        #endregion

        // ──────────────────────────── WavePacket ──────────────────────────────

        #region WavePacket – Constructor Validation

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void WavePacket_TooFewPoints_Throws()
        {
            new WavePacket(new double[] { 0, 1, 2 }, 5, 1, k => k);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void WavePacket_ZeroSigma_Throws()
        {
            var x = Enumerable.Range(0, 64).Select(i => i * 0.1).ToArray();
            new WavePacket(x, 5, 0, k => k);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentNullException))]
        public void WavePacket_NullDispersion_Throws()
        {
            var x = Enumerable.Range(0, 64).Select(i => i * 0.1).ToArray();
            new WavePacket(x, 5, 1, null);
        }

        #endregion

        #region WavePacket – Non-Dispersive Packet Maintains Shape

        [TestMethod]
        public void WavePacket_NonDispersive_WidthConstant()
        {
            // ω = c·k → no dispersion → packet doesn't spread
            double c = 2.0;
            int nPts = 256;
            double dx = 0.1;
            var x = Enumerable.Range(0, nPts).Select(i => i * dx).ToArray();

            var packet = new WavePacket(x, centerK: 5.0, sigma: 2.0,
                dispersion: k => c * k);

            double w0 = packet.Width(0);
            double w1 = packet.Width(1.0);

            // Width shouldn't change much for non-dispersive
            double change = Math.Abs(w1 - w0) / w0;
            Assert.IsTrue(change < 0.1, $"Non-dispersive packet width changed by {change:P1}");
        }

        #endregion

        #region WavePacket – Phase & Group Velocity

        [TestMethod]
        public void WavePacket_PhaseVelocity_LinearDispersion()
        {
            double c = 3.0;
            var x = Enumerable.Range(0, 64).Select(i => i * 0.1).ToArray();
            var packet = new WavePacket(x, centerK: 4.0, sigma: 1.0, dispersion: k => c * k);

            // v_p = ω(k₀)/k₀ = c
            Assert.AreEqual(c, packet.PhaseVelocity, 1e-8);
        }

        [TestMethod]
        public void WavePacket_GroupVelocity_LinearDispersion()
        {
            double c = 3.0;
            var x = Enumerable.Range(0, 64).Select(i => i * 0.1).ToArray();
            var packet = new WavePacket(x, centerK: 4.0, sigma: 1.0, dispersion: k => c * k);

            // v_g = dω/dk = c
            Assert.AreEqual(c, packet.GroupVelocity, 1e-4);
        }

        #endregion

        #region WavePacket – Dispersive Packet Spreads

        [TestMethod]
        public void WavePacket_Dispersive_WidthIncreases()
        {
            // ω = k² → dispersive (like Schrödinger)
            int nPts = 256;
            double dx = 0.1;
            var x = Enumerable.Range(0, nPts).Select(i => i * dx).ToArray();

            var packet = new WavePacket(x, centerK: 3.0, sigma: 2.0,
                dispersion: k => k * k);

            double w0 = packet.Width(0);
            double w2 = packet.Width(2.0);

            Assert.IsTrue(w2 > w0, "Dispersive packet should spread over time.");
        }

        #endregion

        #region WavePacket – Propagate Returns Correct Count

        [TestMethod]
        public void WavePacket_Propagate_CorrectPointCount()
        {
            int nPts = 128;
            var x = Enumerable.Range(0, nPts).Select(i => i * 0.1).ToArray();
            var packet = new WavePacket(x, 5, 1, k => 2 * k);

            var result = packet.Propagate(0.5);
            Assert.AreEqual(nPts, result.Count);
        }

        #endregion

        // ──────────────────────────── FourierSeries ───────────────────────────

        #region FourierSeries – Square Wave Coefficients

        [TestMethod]
        public void Fourier_SquareWave_SineCoefficients()
        {
            // Square wave on [0, 2π): +1 for t∈[0,π), -1 for t∈[π,2π)
            // bₙ = 4/(nπ) for odd n, 0 for even n
            var fs = new FourierSeries();
            double period = 2 * Math.PI;
            Func<double, double> square = t => (t % period) < period / 2 ? 1.0 : -1.0;

            fs.Analyze(square, period, 10, 4096);

            // a₀ ≈ 0 (zero mean)
            Assert.AreEqual(0.0, fs.A0, 0.01, "a₀ should be ~0 for symmetric square wave");

            // b₁ = 4/π ≈ 1.2732
            Assert.AreEqual(4.0 / Math.PI, fs.Bn[0], 0.01, "b₁ incorrect");

            // b₂ ≈ 0 (even)
            Assert.AreEqual(0.0, fs.Bn[1], 0.01, "b₂ should be ~0");

            // b₃ = 4/(3π) ≈ 0.4244
            Assert.AreEqual(4.0 / (3 * Math.PI), fs.Bn[2], 0.01, "b₃ incorrect");
        }

        #endregion

        #region FourierSeries – Pure Sine Analysis

        [TestMethod]
        public void Fourier_PureSine_OnlyB1NonZero()
        {
            var fs = new FourierSeries();
            double period = 1.0;
            Func<double, double> sine = t => Math.Sin(2 * Math.PI * t / period);

            fs.Analyze(sine, period, 5, 2048);

            Assert.AreEqual(0.0, fs.A0, 0.01);
            Assert.AreEqual(0.0, fs.An[0], 0.01, "a₁ should be 0 for pure sine");
            Assert.AreEqual(1.0, Math.Abs(fs.Bn[0]), 0.01, "|b₁| should be ≈ 1");

            // Higher harmonics should be negligible
            for (int n = 1; n < 5; n++)
            {
                Assert.AreEqual(0.0, fs.An[n], 0.01, $"a_{n + 1} should be ~0");
                Assert.AreEqual(0.0, fs.Bn[n], 0.01, $"b_{n + 1} should be ~0");
            }
        }

        #endregion

        #region FourierSeries – Parseval's Theorem

        [TestMethod]
        public void Fourier_Parseval_TimeEqualsFrquencyEnergy()
        {
            var fs = new FourierSeries();
            double period = 2.0;
            Func<double, double> sawtooth = t =>
            {
                double tp = t % period;
                return 2.0 * tp / period - 1.0; // linear ramp -1 to +1
            };

            fs.Analyze(sawtooth, period, 20, 4096);

            double eParseval = fs.ParsevalEnergy();
            double eTime = fs.TimeDomainEnergy(sawtooth, 4096);

            double relErr = Math.Abs(eParseval - eTime) / eTime;
            Assert.IsTrue(relErr < 0.05, $"Parseval error: {relErr:P2}");
        }

        #endregion

        #region FourierSeries – Synthesis Reconstructs Signal

        [TestMethod]
        public void Fourier_Synthesis_ApproximatesOriginal()
        {
            var fs = new FourierSeries();
            double period = 1.0;
            Func<double, double> f = t => Math.Sin(2 * Math.PI * t) + 0.5 * Math.Cos(4 * Math.PI * t);

            fs.Analyze(f, period, 5, 2048);

            // Check at several points
            for (int i = 0; i < 10; i++)
            {
                double t = i * period / 10.0;
                double expected = f(t);
                double actual = fs.Synthesize(t);
                Assert.AreEqual(expected, actual, 0.05, $"Synthesis mismatch at t={t:F2}");
            }
        }

        #endregion

        #region FourierSeries – SynthesizeRange

        [TestMethod]
        public void Fourier_SynthesizeRange_MatchesSynthesize()
        {
            var fs = new FourierSeries();
            fs.Analyze(t => Math.Sin(2 * Math.PI * t), 1.0, 3, 1024);

            var tVals = Enumerable.Range(0, 20).Select(i => i * 0.05).ToArray();
            var series = fs.SynthesizeRange(tVals);

            for (int i = 0; i < tVals.Length; i++)
                Assert.AreEqual(fs.Synthesize(tVals[i]), series[i].Value, 1e-12);
        }

        #endregion

        #region FourierSeries – Constructor Validation

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Fourier_ZeroPeriod_Throws()
        {
            new FourierSeries().Analyze(t => t, 0, 5);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Fourier_ZeroTerms_Throws()
        {
            new FourierSeries().Analyze(t => t, 1.0, 0);
        }

        #endregion

        // ───────────────────── DampedDrivenWaveEquation1D ─────────────────────

        #region Damped Wave – Energy Decays

        [TestMethod]
        public void DampedWave_EnergyDecays()
        {
            int n = 101;
            double dx = 0.01;
            double c = 1.0;
            double alpha = 5.0; // strong damping
            double dt = 0.005;

            var wave = new DampedDrivenWaveEquation1D(n, dx, c, alpha: alpha);
            double L = (n - 1) * dx;

            wave.SetInitialCondition(
                u0: x => Math.Sin(Math.PI * x / L));

            double e0 = wave.TotalEnergy;
            Assert.IsTrue(e0 > 0);

            wave.Simulate(0.5, dt);
            double eFinal = wave.TotalEnergy;

            Assert.IsTrue(eFinal < e0 * 0.5, "Damped energy should drop significantly.");
        }

        #endregion

        #region Damped Wave – Zero Damping ≈ Conservative

        [TestMethod]
        public void DampedWave_ZeroDamping_ConservesEnergy()
        {
            int n = 51;
            double dx = 0.02;
            double c = 1.0;
            double dt = 0.01;

            var wave = new DampedDrivenWaveEquation1D(n, dx, c, alpha: 0);
            double L = (n - 1) * dx;
            wave.SetInitialCondition(u0: x => Math.Sin(Math.PI * x / L));

            double e0 = wave.TotalEnergy;
            wave.Simulate(0.5, dt);
            double eFinal = wave.TotalEnergy;

            double relErr = Math.Abs(eFinal - e0) / e0;
            Assert.IsTrue(relErr < 0.02, $"Zero-damping energy drifted {relErr:P2}");
        }

        #endregion

        #region Driven Wave – Source Injects Energy

        [TestMethod]
        public void DrivenWave_SourceInjectsEnergy()
        {
            int n = 51;
            double dx = 0.02;
            double c = 1.0;
            double dt = 0.01;
            double L = (n - 1) * dx;

            // Continuous sinusoidal driving at the centre
            Func<double, double, double> source = (x, t) =>
            {
                double xc = L / 2.0;
                double env = Math.Exp(-100 * (x - xc) * (x - xc));
                return 10.0 * Math.Sin(20 * t) * env;
            };

            var wave = new DampedDrivenWaveEquation1D(n, dx, c, alpha: 0, source: source);
            wave.SetInitialCondition(); // start from rest

            double e0 = wave.TotalEnergy; // should be 0
            Assert.AreEqual(0.0, e0, 1e-14);

            wave.Simulate(1.0, dt);
            double eFinal = wave.TotalEnergy;

            Assert.IsTrue(eFinal > 0.01, "Source should inject energy into the system.");
        }

        #endregion

        #region Damped+Driven – Steady State Balance

        [TestMethod]
        public void DampedDriven_DampingAndSource_EnergyBounded()
        {
            int n = 51;
            double dx = 0.02;
            double c = 1.0;
            double alpha = 2.0;
            double dt = 0.005;
            double L = (n - 1) * dx;

            Func<double, double, double> source = (x, t) =>
            {
                double xc = L / 2.0;
                return 5.0 * Math.Sin(10 * t) * Math.Exp(-50 * (x - xc) * (x - xc));
            };

            var wave = new DampedDrivenWaveEquation1D(n, dx, c, alpha: alpha, source: source);
            wave.SetInitialCondition();

            // Run to approximate steady state
            wave.Simulate(5.0, dt);
            double eMid = wave.TotalEnergy;

            // Run more to check it doesn't blow up
            wave.Simulate(10.0, dt);
            double eLate = wave.TotalEnergy;

            // Energy should be bounded (not exponentially growing)
            Assert.IsTrue(eLate < 100 * eMid + 1.0,
                "Energy should not blow up in damped+driven system.");
        }

        #endregion

        #region DampedDriven – Constructor Validation

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void DampedDriven_NegativeAlpha_Throws()
        {
            new DampedDrivenWaveEquation1D(10, 0.1, 1.0, alpha: -1);
        }

        #endregion

        #region DampedDriven – Reset

        [TestMethod]
        public void DampedDriven_Reset_RestoresState()
        {
            var wave = new DampedDrivenWaveEquation1D(31, 0.1, 1.0, alpha: 1);
            wave.SetInitialCondition(u0: x => Math.Sin(x));

            var snap0 = wave.Snapshot().Select(s => s.Value).ToArray();
            wave.Simulate(0.5, 0.01);

            wave.Reset();
            var snapReset = wave.Snapshot().Select(s => s.Value).ToArray();

            for (int i = 0; i < snap0.Length; i++)
                Assert.AreEqual(snap0[i], snapReset[i], 1e-14);
            Assert.AreEqual(0.0, wave.Time, 1e-14);
        }

        #endregion
    }
}
