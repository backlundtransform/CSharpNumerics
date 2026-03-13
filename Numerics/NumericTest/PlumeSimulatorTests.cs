using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Simulation;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics;
using CSharpNumerics.Physics.Enums;
using System;
using System.Linq;

namespace NumericsTests
{
    [TestClass]
    public class PlumeSimulatorTests
    {
        // ════════════════════════════════════════════════════════════════
        //  GaussianPuff (transient) — EnvironmentalExtensions
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void GaussianPuff_PeakMovesDownwind()
        {
            double Q = 5.0;
            double release = 10.0;
            double windSpeed = 10.0;
            double H = 50;
            var src = new Vector(0, 0, 50);
            var windDir = new Vector(1, 0, 0);

            // Evaluate at two times: puff centre should move with wind
            var field10 = Q.GaussianPuff(release, windSpeed, H, src, windDir, 10, StabilityClass.D);
            var field60 = Q.GaussianPuff(release, windSpeed, H, src, windDir, 60, StabilityClass.D);

            // At t=10s, puff centre is at x = u*t = 100m downwind
            double c10_at100 = field10.Evaluate(new Vector(100, 0, H));
            double c10_at600 = field10.Evaluate(new Vector(600, 0, H));
            Assert.IsTrue(c10_at100 > c10_at600, "At t=10 concentration should be higher at 100m than 600m");

            // At t=60s, puff centre is at x = 600m downwind
            double c60_at600 = field60.Evaluate(new Vector(600, 0, H));
            double c60_at100 = field60.Evaluate(new Vector(100, 0, H));
            Assert.IsTrue(c60_at600 > c60_at100, "At t=60 concentration should be higher at 600m than 100m");
        }

        [TestMethod]
        public void GaussianPuff_ZeroTime_ReturnsZero()
        {
            var field = 5.0.GaussianPuff(10, 10, 50,
                new Vector(0, 0, 50), new Vector(1, 0, 0), 0, StabilityClass.D);

            Assert.AreEqual(0, field.Evaluate(new Vector(100, 0, 50)));
        }

        [TestMethod]
        public void GaussianPuff_IsPositive()
        {
            var field = 5.0.GaussianPuff(10, 10, 50,
                new Vector(0, 0, 50), new Vector(1, 0, 0), 30, StabilityClass.D);

            // Near the puff centre at x = 10*30 = 300m, z = H = 50
            double c = field.Evaluate(new Vector(300, 0, 50));
            Assert.IsTrue(c > 0, "Concentration at puff centre should be positive.");
        }

        // ════════════════════════════════════════════════════════════════
        //  PlumeSimulator — steady state
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void PlumeSimulator_SteadyState_ProducesSnapshots()
        {
            var sim = new PlumeSimulator(
                emissionRate: 5.0,
                windSpeed: 10,
                windDirection: new Vector(1, 0, 0),
                stackHeight: 50,
                sourcePosition: new Vector(0, 0, 50),
                stability: StabilityClass.D,
                mode: PlumeMode.SteadyState);

            var grid = new GeoGrid(-200, 200, -200, 200, 0, 0, 50);
            var tf = new TimeFrame(0, 300, 60);

            var snaps = sim.Run(grid, tf);

            Assert.AreEqual(tf.Count, snaps.Count);

            // Steady-state → all snapshots should have identical values
            var first = snaps[0].GetValues();
            for (int i = 1; i < snaps.Count; i++)
            {
                var vals = snaps[i].GetValues();
                for (int j = 0; j < first.Length; j++)
                    Assert.AreEqual(first[j], vals[j], 1e-20,
                        $"Snapshot {i}, cell {j} differs from snapshot 0");
            }
        }

        [TestMethod]
        public void PlumeSimulator_SteadyState_ConcentrationHigherDownwind()
        {
            var sim = new PlumeSimulator(
                emissionRate: 5.0,
                windSpeed: 10,
                windDirection: new Vector(1, 0, 0),
                stackHeight: 0,
                sourcePosition: new Vector(0, 0, 0),
                stability: StabilityClass.D,
                mode: PlumeMode.SteadyState);

            // Small grid along wind axis at ground level
            var grid = new GeoGrid(-100, 400, -50, 50, 0, 0, 50);
            var tf = new TimeFrame(0, 0, 1);

            var snap = sim.Run(grid, tf)[0];

            // Upwind cell (x=-100) should be 0
            double upwind = snap.ValueAt(new Vector(-100, 0, 0));
            Assert.AreEqual(0, upwind, "Upwind concentration should be zero.");

            // Downwind cell should be positive
            double downwind = snap.ValueAt(new Vector(200, 0, 0));
            Assert.IsTrue(downwind > 0, "Downwind concentration should be positive.");
        }

        [TestMethod]
        public void PlumeSimulator_SteadyState_LateralDecay()
        {
            var sim = new PlumeSimulator(
                emissionRate: 5.0,
                windSpeed: 10,
                windDirection: new Vector(1, 0, 0),
                stackHeight: 0,
                sourcePosition: new Vector(0, 0, 0),
                stability: StabilityClass.D,
                mode: PlumeMode.SteadyState);

            var grid = new GeoGrid(50, 150, -200, 200, 0, 0, 50);
            var tf = new TimeFrame(0, 0, 1);
            var snap = sim.Run(grid, tf)[0];

            // Centreline (y=0) should be higher than off-axis (y=200)
            double centre = snap.ValueAt(new Vector(100, 0, 0));
            double offAxis = snap.ValueAt(new Vector(100, 200, 0));
            Assert.IsTrue(centre > offAxis,
                $"Centreline={centre:E3} should exceed off-axis={offAxis:E3}");
        }

        // ════════════════════════════════════════════════════════════════
        //  PlumeSimulator — transient
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void PlumeSimulator_Transient_SnapshotsDiffer()
        {
            var sim = new PlumeSimulator(
                emissionRate: 5.0,
                windSpeed: 10,
                windDirection: new Vector(1, 0, 0),
                stackHeight: 50,
                sourcePosition: new Vector(0, 0, 50),
                stability: StabilityClass.D,
                mode: PlumeMode.Transient);
            sim.ReleaseSeconds = 10;

            var grid = new GeoGrid(-200, 800, -200, 200, 0, 100, 50);
            var tf = new TimeFrame(10, 60, 10);

            var snaps = sim.Run(grid, tf);
            Assert.AreEqual(tf.Count, snaps.Count);

            // Snapshots should differ (puff moves/disperses)
            var v0 = snaps[0].GetValues();
            var vLast = snaps[snaps.Count - 1].GetValues();
            bool anyDiff = false;
            for (int j = 0; j < v0.Length; j++)
                if (Math.Abs(v0[j] - vLast[j]) > 1e-30) { anyDiff = true; break; }

            Assert.IsTrue(anyDiff, "Transient snapshots at different times should differ.");
        }

        [TestMethod]
        public void PlumeSimulator_RunSingle_MatchesRunAtSameTime()
        {
            var sim = new PlumeSimulator(
                emissionRate: 5.0,
                windSpeed: 10,
                windDirection: new Vector(1, 0, 0),
                stackHeight: 0,
                sourcePosition: new Vector(0, 0, 0),
                stability: StabilityClass.D,
                mode: PlumeMode.SteadyState);

            var grid = new GeoGrid(0, 100, 0, 100, 0, 0, 50);
            var tf = new TimeFrame(0, 0, 1);

            var fromRun = sim.Run(grid, tf)[0];
            var single = sim.RunSingle(grid, 0);

            var a = fromRun.GetValues();
            var b = single.GetValues();
            for (int i = 0; i < a.Length; i++)
                Assert.AreEqual(a[i], b[i], 1e-20);
        }

        // ════════════════════════════════════════════════════════════════
        //  ScenarioVariation
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void ScenarioVariation_Default_HasNoVariation()
        {
            var sv = new ScenarioVariation();
            Assert.IsFalse(sv.HasVariation);
        }

        [TestMethod]
        public void ScenarioVariation_WindSpeed_SetsRange()
        {
            var sv = new ScenarioVariation().WindSpeed(5, 15);
            Assert.IsTrue(sv.HasVariation);
            Assert.AreEqual(5, sv.WindSpeedMin);
            Assert.AreEqual(15, sv.WindSpeedMax);
        }

        [TestMethod]
        public void ScenarioVariation_Fluent_Chains()
        {
            var sv = new ScenarioVariation()
                .WindSpeed(5, 15)
                .WindDirectionJitter(10)
                .EmissionRate(3, 7)
                .SetStabilityWeights(d: 0.6, c: 0.2, e: 0.2);

            Assert.IsTrue(sv.HasVariation);
            Assert.AreEqual(10, sv.WindDirectionJitterDeg);
            Assert.AreEqual(3, sv.EmissionRateMin);
            Assert.AreEqual(7, sv.EmissionRateMax);
            Assert.IsNotNull(sv.StabilityWeights);
            Assert.AreEqual(0.6, sv.StabilityWeights[StabilityClass.D]);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void ScenarioVariation_WindSpeed_InvertedRange_Throws()
        {
            new ScenarioVariation().WindSpeed(20, 5);
        }

        // ════════════════════════════════════════════════════════════════
        //  PlumeSimulator — validation
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void PlumeSimulator_ZeroWind_Throws()
        {
            new PlumeSimulator(5, 0, new Vector(1, 0, 0), 50, new Vector(0, 0, 50));
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void PlumeSimulator_NegativeEmission_Throws()
        {
            new PlumeSimulator(-1, 10, new Vector(1, 0, 0), 50, new Vector(0, 0, 50));
        }
    }
}
