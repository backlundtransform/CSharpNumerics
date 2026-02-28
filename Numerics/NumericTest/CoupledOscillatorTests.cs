using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Oscillations;
using CSharpNumerics.Statistics.Data;
using System;
using System.Linq;

namespace NumericsTests
{
    [TestClass]
    public class CoupledOscillatorTests
    {
        #region Constructor Validation

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Constructor_EmptyMasses_Throws()
        {
            new CoupledOscillators(new double[0], new double[] { 1 });
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Constructor_NegativeMass_Throws()
        {
            new CoupledOscillators(new[] { -1.0 }, new[] { 1.0, 1.0 });
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Constructor_ZeroMass_Throws()
        {
            new CoupledOscillators(new[] { 0.0 }, new[] { 1.0, 1.0 });
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Constructor_WrongStiffnessCount_Throws()
        {
            new CoupledOscillators(new[] { 1.0, 1.0 }, new[] { 1.0, 1.0 });
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Constructor_NegativeStiffness_Throws()
        {
            new CoupledOscillators(new[] { 1.0 }, new[] { 1.0, -1.0 });
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Constructor_NegativeDamping_Throws()
        {
            new CoupledOscillators(new[] { 1.0 }, new[] { 1.0, 1.0 },
                dampings: new[] { -0.5 });
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Constructor_WrongDampingCount_Throws()
        {
            new CoupledOscillators(new[] { 1.0, 1.0 }, new[] { 1.0, 1.0, 1.0 },
                dampings: new[] { 0.1 });
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Constructor_WrongInitialPositionCount_Throws()
        {
            new CoupledOscillators(new[] { 1.0, 1.0 }, new[] { 1.0, 1.0, 1.0 },
                initialPositions: new[] { 1.0 });
        }

        [TestMethod]
        public void Constructor_ValidParameters_SetsState()
        {
            var osc = new CoupledOscillators(
                new[] { 2.0, 3.0 }, new[] { 10.0, 20.0, 30.0 },
                initialPositions: new[] { 0.5, -0.3 },
                initialVelocities: new[] { 0.1, 0.2 });

            Assert.AreEqual(2, osc.Count);
            Assert.AreEqual(0.0, osc.Time, 1e-12);
            Assert.AreEqual(0.5, osc.Position(0), 1e-12);
            Assert.AreEqual(-0.3, osc.Position(1), 1e-12);
            Assert.AreEqual(0.1, osc.Velocity(0), 1e-12);
            Assert.AreEqual(0.2, osc.Velocity(1), 1e-12);
        }

        [TestMethod]
        public void UniformConstructor_CreatesCorrectArrays()
        {
            var osc = new CoupledOscillators(3, mass: 2.0, stiffness: 5.0, damping: 0.1);

            Assert.AreEqual(3, osc.Count);
            CollectionAssert.AreEqual(new[] { 2.0, 2.0, 2.0 }, osc.Masses);
            CollectionAssert.AreEqual(new[] { 5.0, 5.0, 5.0, 5.0 }, osc.Stiffnesses);
            CollectionAssert.AreEqual(new[] { 0.1, 0.1, 0.1 }, osc.Dampings);
        }

        #endregion

        #region Matrices

        [TestMethod]
        public void StiffnessMatrix_2Mass_IsCorrect()
        {
            // wall-k-m-k-m-k-wall, all k equal
            var osc = new CoupledOscillators(2, mass: 1.0, stiffness: 5.0);
            var K = osc.StiffnessMatrix();

            Assert.AreEqual(10.0, K.values[0, 0], 1e-12); // k0+k1
            Assert.AreEqual(-5.0, K.values[0, 1], 1e-12);
            Assert.AreEqual(-5.0, K.values[1, 0], 1e-12);
            Assert.AreEqual(10.0, K.values[1, 1], 1e-12); // k1+k2
        }

        [TestMethod]
        public void MassMatrix_IsDiagonal()
        {
            var osc = new CoupledOscillators(
                new[] { 2.0, 3.0, 5.0 }, new[] { 1.0, 1.0, 1.0, 1.0 });
            var M = osc.MassMatrix();

            Assert.AreEqual(2.0, M.values[0, 0], 1e-12);
            Assert.AreEqual(3.0, M.values[1, 1], 1e-12);
            Assert.AreEqual(5.0, M.values[2, 2], 1e-12);
            Assert.AreEqual(0.0, M.values[0, 1], 1e-12);
            Assert.AreEqual(0.0, M.values[1, 2], 1e-12);
        }

        #endregion

        #region Normal Modes — 1-mass system

        [TestMethod]
        public void NormalModes_1Mass_IsCorrect()
        {
            // Single mass between two walls: ω = √((k0+k1)/m)
            double m = 2.0, k0 = 3.0, k1 = 5.0;
            var osc = new CoupledOscillators(new[] { m }, new[] { k0, k1 });

            var modes = osc.NormalModes();
            Assert.AreEqual(1, modes.Count);
            Assert.AreEqual(Math.Sqrt((k0 + k1) / m), modes[0], 1e-10);
        }

        #endregion

        #region Normal Modes — 2-mass uniform system

        [TestMethod]
        public void NormalModes_2MassUniform_MatchesAnalytical()
        {
            // Uniform: m=1, k=1
            // K = [[2, -1], [-1, 2]] → eigenvalues of K/m: 1 and 3
            // ω₁ = 1, ω₂ = √3
            var osc = new CoupledOscillators(2, mass: 1.0, stiffness: 1.0);
            var modes = osc.NormalModes();

            Assert.AreEqual(2, modes.Count);
            Assert.AreEqual(1.0, modes[0], 1e-8, $"ω₁ = {modes[0]:F8}, expected 1.0");
            Assert.AreEqual(Math.Sqrt(3.0), modes[1], 1e-8,
                $"ω₂ = {modes[1]:F8}, expected {Math.Sqrt(3.0):F8}");
        }

        [TestMethod]
        public void NormalModes_2MassUniform_AlternativeK()
        {
            // m=2, k=8 → K/m eigenvalues: 8/2=4 and 24/2=12
            // ω₁ = 2, ω₂ = 2√3
            var osc = new CoupledOscillators(2, mass: 2.0, stiffness: 8.0);
            var modes = osc.NormalModes();

            Assert.AreEqual(2.0, modes[0], 1e-8);
            Assert.AreEqual(2.0 * Math.Sqrt(3.0), modes[1], 1e-8);
        }

        #endregion

        #region Mode Shapes — 2-mass uniform

        [TestMethod]
        public void ModeShapes_2MassUniform_Symmetric()
        {
            var osc = new CoupledOscillators(2, mass: 1.0, stiffness: 1.0);
            var shapes = osc.ModeShapes();

            Assert.AreEqual(2, shapes.Count);

            // Mode 1 (lower freq): symmetric [1, 1]/√2
            double expectedComponent = 1.0 / Math.Sqrt(2.0);
            Assert.AreEqual(Math.Abs(shapes[0][0]), Math.Abs(shapes[0][1]), 1e-8,
                "Mode 1 should be symmetric: |φ₁[0]| = |φ₁[1]|");
            Assert.AreEqual(expectedComponent, Math.Abs(shapes[0][0]), 1e-8);
            // Both components should have the SAME sign
            Assert.IsTrue(shapes[0][0] * shapes[0][1] > 0,
                $"Mode 1 components should have same sign: {shapes[0][0]:F6}, {shapes[0][1]:F6}");

            // Mode 2 (higher freq): antisymmetric [1, -1]/√2
            Assert.AreEqual(Math.Abs(shapes[1][0]), Math.Abs(shapes[1][1]), 1e-8,
                "Mode 2 should be antisymmetric: |φ₂[0]| = |φ₂[1]|");
            Assert.IsTrue(shapes[1][0] * shapes[1][1] < 0,
                $"Mode 2 components should have opposite sign: {shapes[1][0]:F6}, {shapes[1][1]:F6}");
        }

        [TestMethod]
        public void ModeShapes_AreNormalised()
        {
            var osc = new CoupledOscillators(3, mass: 1.0, stiffness: 2.0);
            var shapes = osc.ModeShapes();

            foreach (var shape in shapes)
            {
                double norm = shape.Norm();
                Assert.AreEqual(1.0, norm, 1e-8, $"Mode shape norm = {norm}");
            }
        }

        #endregion

        #region Normal Modes — 3-mass uniform (analytical)

        [TestMethod]
        public void NormalModes_3MassUniform_MatchesAnalytical()
        {
            // For uniform chain with N masses, fixed walls:
            // ω_r = 2√(k/m) sin(rπ / (2(N+1)))   r = 1,2,...,N
            double m = 1.0, k = 4.0;
            int N = 3;
            var osc = new CoupledOscillators(N, mass: m, stiffness: k);

            var expected = Enumerable.Range(1, N)
                .Select(r => 2.0 * Math.Sqrt(k / m) * Math.Sin(r * Math.PI / (2.0 * (N + 1))))
                .ToList();

            var modes = osc.NormalModes();

            for (int r = 0; r < N; r++)
            {
                Assert.AreEqual(expected[r], modes[r], 1e-6,
                    $"Mode {r + 1}: expected ω={expected[r]:F6}, got {modes[r]:F6}");
            }
        }

        #endregion

        #region Normal Modes — 5-mass uniform

        [TestMethod]
        public void NormalModes_5MassUniform_MatchesAnalytical()
        {
            double m = 2.0, k = 10.0;
            int N = 5;
            var osc = new CoupledOscillators(N, mass: m, stiffness: k);

            var expected = Enumerable.Range(1, N)
                .Select(r => 2.0 * Math.Sqrt(k / m) * Math.Sin(r * Math.PI / (2.0 * (N + 1))))
                .ToList();

            var modes = osc.NormalModes();

            for (int r = 0; r < N; r++)
            {
                Assert.AreEqual(expected[r], modes[r], 1e-5,
                    $"Mode {r + 1}: expected ω={expected[r]:F6}, got {modes[r]:F6}");
            }
        }

        #endregion

        #region Energy Conservation (Undamped)

        [TestMethod]
        public void UndampedSystem_EnergyConserved()
        {
            var osc = new CoupledOscillators(3, mass: 1.0, stiffness: 5.0,
                initialPositions: new[] { 1.0, 0.0, -0.5 });

            double E0 = osc.TotalEnergy;
            double dt = 0.001;

            for (int i = 0; i < 10000; i++)
                osc.Step(dt);

            Assert.AreEqual(E0, osc.TotalEnergy, E0 * 1e-6,
                $"Energy at t={osc.Time}: {osc.TotalEnergy:F8}, initial: {E0:F8}");
        }

        [TestMethod]
        public void UndampedSystem_EnergyOverTime_IsConstant()
        {
            var osc = new CoupledOscillators(2, mass: 1.0, stiffness: 3.0,
                initialPositions: new[] { 1.0, -1.0 });

            var energy = osc.EnergyOverTime(5.0, 0.001);
            double E0 = energy.First().Value;
            double maxDeviation = energy.Max(e => Math.Abs(e.Value - E0));

            Assert.IsTrue(maxDeviation < E0 * 1e-6,
                $"Max energy deviation: {maxDeviation}, initial: {E0}");
        }

        #endregion

        #region Energy Dissipation (Damped)

        [TestMethod]
        public void DampedSystem_EnergyDecreases()
        {
            var osc = new CoupledOscillators(2, mass: 1.0, stiffness: 5.0, damping: 0.5,
                initialPositions: new[] { 1.0, 0.0 });

            var energy = osc.EnergyOverTime(5.0, 0.001);
            double firstEnergy = energy.First().Value;
            double lastEnergy = energy.Last().Value;

            Assert.IsTrue(lastEnergy < firstEnergy,
                $"Energy should decrease: first={firstEnergy:F6}, last={lastEnergy:F6}");
        }

        [TestMethod]
        public void DampedSystem_EnergyMonotonicallyDecreasing()
        {
            var osc = new CoupledOscillators(2, mass: 1.0, stiffness: 5.0, damping: 1.0,
                initialPositions: new[] { 1.0, 0.5 });

            var energy = osc.EnergyOverTime(3.0, 0.001);

            // Check that energy never increases (with small tolerance for RK4 rounding)
            for (int i = 1; i < energy.Count; i++)
            {
                Assert.IsTrue(energy[i].Value <= energy[i - 1].Value + 1e-10,
                    $"Energy increased at t={energy[i].Index:F4}: " +
                    $"{energy[i].Value:F8} > {energy[i - 1].Value:F8}");
            }
        }

        #endregion

        #region Symmetric Initial Conditions → Single Mode Excitation

        [TestMethod]
        public void SymmetricIC_ExcitesOnlySymmetricMode()
        {
            // [1, 1] initial condition: only the symmetric mode (ω₁=1) should be excited
            var osc = new CoupledOscillators(2, mass: 1.0, stiffness: 1.0,
                initialPositions: new[] { 1.0, 1.0 });

            double omega1 = 1.0; // symmetric mode frequency
            double period = 2.0 * Math.PI / omega1;

            var traj = osc.Trajectory(period, 0.001);

            // Both masses should oscillate in phase at ω₁
            // After one full period, both should return near initial positions
            var mass0End = traj[0].Last().Value;
            var mass1End = traj[1].Last().Value;

            Assert.AreEqual(1.0, mass0End, 0.01,
                $"Mass 0 after one period: {mass0End:F4}");
            Assert.AreEqual(1.0, mass1End, 0.01,
                $"Mass 1 after one period: {mass1End:F4}");
        }

        [TestMethod]
        public void AntisymmetricIC_ExcitesOnlyAntisymmetricMode()
        {
            // [1, -1] initial condition: only the antisymmetric mode (ω₂=√3) should be excited
            var osc = new CoupledOscillators(2, mass: 1.0, stiffness: 1.0,
                initialPositions: new[] { 1.0, -1.0 });

            double omega2 = Math.Sqrt(3.0);
            double period = 2.0 * Math.PI / omega2;

            var traj = osc.Trajectory(period, 0.001);

            // After one period of the antisymmetric mode, should return
            var mass0End = traj[0].Last().Value;
            var mass1End = traj[1].Last().Value;

            Assert.AreEqual(1.0, mass0End, 0.01,
                $"Mass 0 after one period: {mass0End:F4}");
            Assert.AreEqual(-1.0, mass1End, 0.01,
                $"Mass 1 after one period: {mass1End:F4}");
        }

        #endregion

        #region Modal Energy

        [TestMethod]
        public void ModalEnergy_SumApproximatesTotalEnergy()
        {
            var osc = new CoupledOscillators(3, mass: 1.0, stiffness: 2.0,
                initialPositions: new[] { 1.0, 0.5, -0.3 });

            double tEnd = 3.0, dt = 0.001;

            // Get total energy
            var totalEnergy = osc.EnergyOverTime(tEnd, dt);

            // Sum modal energies
            var modalSum = new double[totalEnergy.Count];
            for (int mode = 0; mode < 3; mode++)
            {
                var mEnergy = osc.ModalEnergy(mode, tEnd, dt);
                for (int i = 0; i < modalSum.Length && i < mEnergy.Count; i++)
                    modalSum[i] += mEnergy[i].Value;
            }

            // Compare at several time points
            for (int i = 0; i < totalEnergy.Count; i += totalEnergy.Count / 10)
            {
                Assert.AreEqual(totalEnergy[i].Value, modalSum[i],
                    totalEnergy[i].Value * 0.01,
                    $"Modal sum ≠ total at t={totalEnergy[i].Index:F3}: " +
                    $"sum={modalSum[i]:F6}, total={totalEnergy[i].Value:F6}");
            }
        }

        [TestMethod]
        public void ModalEnergy_SymmetricIC_OnlyMode1HasEnergy()
        {
            // Excite only the symmetric mode
            var osc = new CoupledOscillators(2, mass: 1.0, stiffness: 1.0,
                initialPositions: new[] { 1.0, 1.0 });

            double tEnd = 5.0, dt = 0.001;

            var mode0Energy = osc.ModalEnergy(0, tEnd, dt);
            var mode1Energy = osc.ModalEnergy(1, tEnd, dt);

            double avgE0 = mode0Energy.Average(e => e.Value);
            double avgE1 = mode1Energy.Average(e => e.Value);

            Assert.IsTrue(avgE0 > 0.1,
                $"Mode 0 should have significant energy: {avgE0:F6}");
            Assert.IsTrue(avgE1 < avgE0 * 0.01,
                $"Mode 1 should have negligible energy: {avgE1:F6} vs mode0={avgE0:F6}");
        }

        #endregion

        #region Dispersion Relation

        [TestMethod]
        public void DispersionRelation_UniformChain_MatchesFormula()
        {
            double m = 1.0, k = 4.0;
            var osc = new CoupledOscillators(5, mass: m, stiffness: k);

            var kValues = Enumerable.Range(0, 100)
                .Select(i => i * Math.PI / 100.0)
                .ToArray();

            var dispersion = osc.DispersionRelation(kValues);

            for (int i = 0; i < dispersion.Count; i++)
            {
                double expected = 2.0 * Math.Sqrt(k / m) * Math.Abs(Math.Sin(kValues[i] / 2.0));
                Assert.AreEqual(expected, dispersion[i].Value, 1e-12,
                    $"ω at k={kValues[i]:F3}");
            }
        }

        [TestMethod]
        public void DispersionRelation_AtZero_IsZero()
        {
            var osc = new CoupledOscillators(3, mass: 1.0, stiffness: 1.0);
            var disp = osc.DispersionRelation(new[] { 0.0 });
            Assert.AreEqual(0.0, disp[0].Value, 1e-12);
        }

        [TestMethod]
        public void DispersionRelation_AtBrillouinZoneEdge_IsMaximum()
        {
            double m = 1.0, k = 4.0;
            var osc = new CoupledOscillators(5, mass: m, stiffness: k);

            // ω_max at k = π/a (a=1) → ω = 2√(k/m)
            var disp = osc.DispersionRelation(new[] { Math.PI });
            Assert.AreEqual(2.0 * Math.Sqrt(k / m), disp[0].Value, 1e-12);
        }

        #endregion

        #region Phase & Group Velocity

        [TestMethod]
        public void PhaseVelocity_AtZero_EqualsLongWavelengthLimit()
        {
            double m = 2.0, k = 8.0, a = 1.0;
            var osc = new CoupledOscillators(5, mass: m, stiffness: k);

            // v_p(k→0) = a√(k/m) = √4 = 2
            Assert.AreEqual(a * Math.Sqrt(k / m), osc.PhaseVelocity(0, a), 1e-12);
        }

        [TestMethod]
        public void GroupVelocity_AtZero_EqualsPhaseVelocity()
        {
            double m = 2.0, k = 8.0, a = 1.0;
            var osc = new CoupledOscillators(5, mass: m, stiffness: k);

            // At k=0, vg = vp (no dispersion)
            double vp = osc.PhaseVelocity(0.001, a);
            double vg = osc.GroupVelocity(0.001, a);

            Assert.AreEqual(vp, vg, 0.01,
                $"At k≈0: vp={vp:F4}, vg={vg:F4}");
        }

        [TestMethod]
        public void GroupVelocity_AtBrillouinEdge_ApproachesZero()
        {
            double m = 1.0, k = 4.0, a = 1.0;
            var osc = new CoupledOscillators(5, mass: m, stiffness: k);

            // At k = π/a, dω/dk = a√(k/m)cos(π/2) = 0
            double vg = osc.GroupVelocity(Math.PI - 0.001, a);
            Assert.IsTrue(Math.Abs(vg) < 0.1,
                $"Group velocity near BZ edge should be small: {vg:F4}");
        }

        #endregion

        #region Trajectory & Reset

        [TestMethod]
        public void Trajectory_StartsAtInitialCondition()
        {
            var osc = new CoupledOscillators(2, mass: 1.0, stiffness: 2.0,
                initialPositions: new[] { 0.5, -0.3 });

            var traj = osc.Trajectory(1.0, 0.01);

            Assert.AreEqual(2, traj.Count);
            Assert.AreEqual(0.0, traj[0].First().Index, 1e-12);
            Assert.AreEqual(0.5, traj[0].First().Value, 1e-12);
            Assert.AreEqual(-0.3, traj[1].First().Value, 1e-12);
        }

        [TestMethod]
        public void Trajectory_CorrectPointCount()
        {
            var osc = new CoupledOscillators(2, mass: 1.0, stiffness: 1.0,
                initialPositions: new[] { 1.0, 0.0 });

            var traj = osc.Trajectory(1.0, 0.01);

            Assert.AreEqual(101, traj[0].Count);
            Assert.AreEqual(101, traj[1].Count);
        }

        [TestMethod]
        public void Trajectory_ResetsState()
        {
            var osc = new CoupledOscillators(2, mass: 1.0, stiffness: 1.0,
                initialPositions: new[] { 1.0, -0.5 },
                initialVelocities: new[] { 0.1, 0.2 });

            osc.Trajectory(1.0, 0.01);

            Assert.AreEqual(0.0, osc.Time, 1e-12);
            Assert.AreEqual(1.0, osc.Position(0), 1e-12);
            Assert.AreEqual(-0.5, osc.Position(1), 1e-12);
        }

        [TestMethod]
        public void Reset_RestoresInitialState()
        {
            var osc = new CoupledOscillators(2, mass: 1.0, stiffness: 1.0,
                initialPositions: new[] { 1.0, 0.0 });

            osc.Step(0.1);
            osc.Step(0.1);
            Assert.AreNotEqual(1.0, osc.Position(0));

            osc.Reset();
            Assert.AreEqual(0.0, osc.Time, 1e-12);
            Assert.AreEqual(1.0, osc.Position(0), 1e-12);
            Assert.AreEqual(0.0, osc.Position(1), 1e-12);
        }

        #endregion

        #region Phase Portrait

        [TestMethod]
        public void PhasePortrait_CorrectStructure()
        {
            var osc = new CoupledOscillators(2, mass: 1.0, stiffness: 1.0,
                initialPositions: new[] { 1.0, 0.0 });

            var pp = osc.PhasePortrait(0, 1.0, 0.01);

            Assert.IsTrue(pp.Count >= 100);
            // First point: position=1, velocity=0
            Assert.AreEqual(1.0, pp.First().Index, 1e-12);
            Assert.AreEqual(0.0, pp.First().Value, 1e-12);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void PhasePortrait_InvalidIndex_Throws()
        {
            var osc = new CoupledOscillators(2, mass: 1.0, stiffness: 1.0);
            osc.PhasePortrait(5, 1.0, 0.01);
        }

        #endregion

        #region FFT — Peaks at Normal Mode Frequencies

        [TestMethod]
        public void FrequencySpectrum_2Mass_PeaksAtNormalModes()
        {
            // m=1, k=1: ω₁=1, ω₂=√3
            // f₁ = 1/(2π) ≈ 0.159 Hz, f₂ = √3/(2π) ≈ 0.276 Hz
            // Excite both modes with [1, 0] initial condition
            var osc = new CoupledOscillators(2, mass: 1.0, stiffness: 1.0,
                initialPositions: new[] { 1.0, 0.0 });

            double f1 = 1.0 / (2.0 * Math.PI);
            double f2 = Math.Sqrt(3.0) / (2.0 * Math.PI);

            var spectrum = osc.FrequencySpectrum(0, 50.0, 0.01);

            // Find the two highest peaks (skip DC)
            var peaks = spectrum.Skip(1).Take(spectrum.Count / 2)
                .OrderByDescending(s => s.Value)
                .Take(2)
                .OrderBy(s => s.Index)
                .ToList();

            double freqRes = (1.0 / 0.01) / spectrum.Count;

            Assert.AreEqual(f1, peaks[0].Index, freqRes * 3,
                $"Peak 1 at {peaks[0].Index:F4}, expected {f1:F4}");
            Assert.AreEqual(f2, peaks[1].Index, freqRes * 3,
                $"Peak 2 at {peaks[1].Index:F4}, expected {f2:F4}");
        }

        #endregion

        #region Numerical vs Analytical — 2-mass beat pattern

        [TestMethod]
        public void TwoMass_BeatPattern_CorrectFrequencies()
        {
            // [1, 0] initial condition decomposes as 0.5*[1,1] + 0.5*[1,-1]
            // Mass 0: x₀(t) = 0.5*cos(ω₁t) + 0.5*cos(ω₂t)
            // Mass 1: x₁(t) = 0.5*cos(ω₁t) - 0.5*cos(ω₂t)
            // Beat envelope of x₀ has period T_beat = 2π/(ω₂−ω₁)
            // Energy oscillates between the two masses
            var osc = new CoupledOscillators(2, mass: 1.0, stiffness: 1.0,
                initialPositions: new[] { 1.0, 0.0 });

            double omega1 = 1.0, omega2 = Math.Sqrt(3.0);

            // Verify analytical trajectory of mass 0 at several time points
            double dt = 0.0001;
            double tEnd = 2.0 * Math.PI / (omega2 - omega1); // one full beat period
            int steps = (int)(tEnd / dt);

            double maxDeviation = 0;
            for (int i = 0; i < steps; i++)
            {
                double t = osc.Time;
                double expected = 0.5 * Math.Cos(omega1 * t) + 0.5 * Math.Cos(omega2 * t);
                double deviation = Math.Abs(osc.Position(0) - expected);
                if (deviation > maxDeviation) maxDeviation = deviation;
                osc.Step(dt);
            }

            Assert.IsTrue(maxDeviation < 0.01,
                $"Max deviation from analytical: {maxDeviation:E3}");
        }

        #endregion

        #region Non-uniform masses

        [TestMethod]
        public void NonUniformMasses_EnergyConserved()
        {
            var osc = new CoupledOscillators(
                new[] { 1.0, 2.0, 3.0 },
                new[] { 5.0, 10.0, 8.0, 3.0 },
                initialPositions: new[] { 0.5, 0.0, -0.3 });

            double E0 = osc.TotalEnergy;
            double dt = 0.001;

            for (int i = 0; i < 5000; i++) osc.Step(dt);

            Assert.AreEqual(E0, osc.TotalEnergy, E0 * 1e-5,
                $"Energy at t={osc.Time}: {osc.TotalEnergy:F8}, initial: {E0:F8}");
        }

        [TestMethod]
        public void NonUniformMasses_NormalModesCount()
        {
            var osc = new CoupledOscillators(
                new[] { 1.0, 3.0, 2.0 },
                new[] { 5.0, 10.0, 8.0, 3.0 });

            var modes = osc.NormalModes();
            Assert.AreEqual(3, modes.Count);

            // All frequencies should be positive and in ascending order
            Assert.IsTrue(modes[0] > 0);
            Assert.IsTrue(modes[1] > modes[0]);
            Assert.IsTrue(modes[2] > modes[1]);
        }

        #endregion

        #region ToString

        [TestMethod]
        public void ToString_ContainsKeyInfo()
        {
            var osc = new CoupledOscillators(3, mass: 2.0, stiffness: 5.0);
            string s = osc.ToString();
            Assert.IsTrue(s.Contains("CoupledOscillators"));
            Assert.IsTrue(s.Contains("N=3"));
            Assert.IsTrue(s.Contains("kg"));
            Assert.IsTrue(s.Contains("N/m"));
        }

        [TestMethod]
        public void ToString_DampedSystem_ShowsDamping()
        {
            var osc = new CoupledOscillators(2, mass: 1.0, stiffness: 3.0, damping: 0.5);
            string s = osc.ToString();
            Assert.IsTrue(s.Contains("Ns/m"), $"Missing damping in: {s}");
        }

        #endregion
    }
}
