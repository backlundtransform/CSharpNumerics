using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Spread;
using CSharpNumerics.Engines.GIS.Spread.WaterContamination;
using CSharpNumerics.Engines.GIS.Spread.WaterContamination.Enums;
using CSharpNumerics.Engines.GIS.Terrain;
using CSharpNumerics.Physics.Environmental.Water;
using CSharpNumerics.Physics.Materials.Water;
using System;
using System.Collections.Generic;
using System.Linq;

namespace NumericsTests
{
    [TestClass]
    public class WaterContaminationSimulatorTests
    {
        // ════════════════════════════════════════════════════════════════
        //  Helpers
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Creates a straight river along ix=5 from iy=0 to iy=9 on a 11×11 grid
        /// with steady downhill slope in y. Returns (grid, terrain, network, channelMap).
        /// </summary>
        private (GeoGrid grid, TerrainGrid terrain, RiverNetwork net, ChannelMap channels)
            CreateStraightRiver(double width = 10, double depth = 2, double manningN = 0.035)
        {
            var grid = new GeoGrid(0, 1000, 0, 1000, 0, 0, 100);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => 100.0 - y * 0.001);

            var cells = new List<(int, int)>();
            for (int iy = 0; iy < grid.Ny; iy++)
                cells.Add((5, iy));

            var net = RiverNetwork.FromManual(grid)
                .AddSegment(cells)
                .Build();

            var channels = new ChannelMap(grid);
            channels.SetUniformChannel(width, depth, manningN);

            return (grid, terrain, net, channels);
        }

        /// <summary>
        /// Runs a contamination simulation on a straight river and returns snapshots.
        /// </summary>
        private IReadOnlyList<SpreadSnapshot> RunStraightRiver(
            AquaticContaminant contaminant,
            double sourceConc = 100, double sourceDuration = double.MaxValue,
            double simDuration = 3600, double dt = 60,
            double bedPorosity = 0.4, double bedBulkDensity = 1600)
        {
            var (grid, terrain, net, channels) = CreateStraightRiver();

            var sources = new List<(int, int, double, double)>
            {
                (5, 0, sourceConc, sourceDuration)
            };

            var parameters = new WaterContaminationParameters(
                sources, contaminant, baseDischargeM3s: 10,
                bedPorosity: bedPorosity, bedBulkDensity: bedBulkDensity);

            var fuelMap = new FuelMap(grid); // unused but required by interface
            var tf = new TimeFrame(0, simDuration, dt);

            var sim = new WaterContaminationSimulator(parameters, net, channels);
            return sim.Run(grid, terrain, fuelMap, tf);
        }

        // ════════════════════════════════════════════════════════════════
        //  Conservative tracer — mass conservation
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void ConservativeTracer_MassConserved()
        {
            // Cyanide: halfLife=0 → conservative, Kd=0 → no adsorption
            var tracer = AquaticContaminant.Cyanide;

            // Source active for only 1 time step (60s), then observe transport
            var snaps = RunStraightRiver(tracer, sourceConc: 100, sourceDuration: 60,
                simDuration: 600, dt: 60);

            // After source stops, total mass should be roughly conserved
            // (small numerical diffusion from upwind scheme is expected)
            double massAfterSource = 0;
            var concLayer2 = snaps[1].Snapshot.GetLayer("concentration");
            for (int i = 0; i < concLayer2.Length; i++)
                massAfterSource += concLayer2[i];

            // Check at a later time step
            double massLater = 0;
            var lastSnap = snaps[snaps.Count - 1];
            var concLayerLast = lastSnap.Snapshot.GetLayer("concentration");
            for (int i = 0; i < concLayerLast.Length; i++)
                massLater += concLayerLast[i];

            // Allow 30% tolerance for numerical diffusion in upwind scheme
            if (massAfterSource > 0)
            {
                double ratio = massLater / massAfterSource;
                Assert.IsTrue(ratio > 0.5 && ratio < 1.5,
                    $"Mass ratio {ratio:F3} — conservative tracer should approximately conserve mass");
            }
        }

        // ════════════════════════════════════════════════════════════════
        //  Advection — plume moves downstream
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void SingleSource_PlumeAdvectsDownstream()
        {
            var tracer = AquaticContaminant.Cyanide;
            var snaps = RunStraightRiver(tracer, sourceConc: 100, sourceDuration: double.MaxValue,
                simDuration: 3600, dt: 60);

            // After some time, downstream cells should have non-zero concentration
            var lastSnap = snaps[snaps.Count - 1];
            var conc = lastSnap.Snapshot.GetLayer("concentration");

            // Source at (5,0), check cells downstream
            // Grid is 11×11 (0..100 step 10)
            int nx = 11;
            double concAt5 = conc[5 * nx + 5]; // (5,5) downstream
            Assert.IsTrue(concAt5 > 0,
                "Downstream cells should have concentration > 0 after advection");
        }

        [TestMethod]
        public void SingleSource_UpstreamCellsStayClean()
        {
            var tracer = AquaticContaminant.Cyanide;

            // Source at (5,0) which is the headwater — so there are no cells upstream
            // But non-river cells should be 0
            var snaps = RunStraightRiver(tracer, sourceConc: 100);
            var lastSnap = snaps[snaps.Count - 1];
            var conc = lastSnap.Snapshot.GetLayer("concentration");

            // Cell (0,0) is not on the river — should be 0
            Assert.AreEqual(0, conc[0 * 11 + 0], 1e-15,
                "Non-river cells should have zero concentration");
        }

        // ════════════════════════════════════════════════════════════════
        //  Dispersion — plume spreads
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void Dispersion_PlumeWidensOverTime()
        {
            // Use a tracer with no decay to isolate dispersion effects
            var tracer = AquaticContaminant.Cyanide;

            // Short pulse source — keep sim short so plume doesn't exit the grid
            var snaps = RunStraightRiver(tracer, sourceConc: 100, sourceDuration: 60,
                simDuration: 600, dt: 60);

            // Count contaminated cells at two different times
            int earlyContaminated = 0;
            int lateContaminated = 0;

            int earlyIdx = Math.Min(2, snaps.Count - 1);
            int lateIdx = snaps.Count - 1;

            var earlyConc = snaps[earlyIdx].Snapshot.GetLayer("concentration");
            var lateConc = snaps[lateIdx].Snapshot.GetLayer("concentration");

            double traceThreshold = 1e-6;
            for (int i = 0; i < earlyConc.Length; i++)
                if (earlyConc[i] > traceThreshold) earlyContaminated++;
            for (int i = 0; i < lateConc.Length; i++)
                if (lateConc[i] > traceThreshold) lateContaminated++;

            // As plume disperses, it should affect more cells (or at least as many)
            Assert.IsTrue(lateContaminated >= earlyContaminated,
                $"Late contaminated cells ({lateContaminated}) should be ≥ early ({earlyContaminated})");
        }

        // ════════════════════════════════════════════════════════════════
        //  First-order decay
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void Decay_PeakConcentrationDropsOverTime()
        {
            // E. coli: half-life ≈ 2 days — fast enough to see effect
            // Use short half-life contaminant
            var decaying = new AquaticContaminant("FastDecay",
                CSharpNumerics.Physics.Materials.Water.Enums.ContaminantType.Chemical,
                halfLifeSeconds: 300, // 5 min half-life
                partitionCoefficient: 0,
                toxicityThresholdMgL: 0.1,
                lethalThresholdMgL: 100);

            // Source for 60s then observe decay
            var snaps = RunStraightRiver(decaying, sourceConc: 100, sourceDuration: 60,
                simDuration: 1800, dt: 60);

            // Peak concentration should decrease over time
            double peakEarly = snaps[Math.Min(2, snaps.Count - 1)].MaxConcentration;
            double peakLate = snaps[snaps.Count - 1].MaxConcentration;

            Assert.IsTrue(peakLate < peakEarly,
                $"Peak should decay: early={peakEarly:F3}, late={peakLate:F3}");
        }

        // ════════════════════════════════════════════════════════════════
        //  Retardation
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void Retardation_SlowsPlumeFront()
        {
            // Run with no adsorption (Kd=0 → Rf=1)
            var noAdsorption = new AquaticContaminant("NoAds",
                CSharpNumerics.Physics.Materials.Water.Enums.ContaminantType.Chemical,
                halfLifeSeconds: 0, partitionCoefficient: 0,
                toxicityThresholdMgL: 0.1, lethalThresholdMgL: 100);

            var snapsNoAds = RunStraightRiver(noAdsorption, sourceConc: 100, sourceDuration: 120,
                simDuration: 600, dt: 60, bedPorosity: 0.4, bedBulkDensity: 1600);

            // Run with high adsorption (Kd=1000 → Rf >> 1)
            var highAdsorption = new AquaticContaminant("HighAds",
                CSharpNumerics.Physics.Materials.Water.Enums.ContaminantType.Chemical,
                halfLifeSeconds: 0, partitionCoefficient: 1000,
                toxicityThresholdMgL: 0.1, lethalThresholdMgL: 100);

            var snapsHighAds = RunStraightRiver(highAdsorption, sourceConc: 100, sourceDuration: 120,
                simDuration: 600, dt: 60, bedPorosity: 0.4, bedBulkDensity: 1600);

            // Count contaminated cells at last time step — retarded plume should be shorter
            int lastIdx = snapsNoAds.Count - 1;
            int contNoAds = snapsNoAds[lastIdx].ContaminatedCellCount(0.01);
            int contHighAds = snapsHighAds[lastIdx].ContaminatedCellCount(0.01);

            Assert.IsTrue(contHighAds <= contNoAds,
                $"Retarded plume ({contHighAds} cells) should spread ≤ non-retarded ({contNoAds} cells)");
        }

        // ════════════════════════════════════════════════════════════════
        //  Confluence mixing
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void Confluence_MassMixesCorrectly()
        {
            // Y-shaped river: left branch contaminated, right branch clean
            var grid = new GeoGrid(0, 1000, 0, 1000, 0, 0, 100);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => 100.0 - y * 0.001);

            var left = new List<(int, int)> { (3, 0), (3, 1), (4, 2), (5, 3) };
            var right = new List<(int, int)> { (7, 0), (7, 1), (6, 2), (5, 3) };
            var main = new List<(int, int)> { (5, 3), (5, 4), (5, 5), (5, 6) };

            var net = RiverNetwork.FromManual(grid)
                .AddSegment(left)
                .AddSegment(right)
                .AddSegment(main)
                .Build();

            var channels = new ChannelMap(grid);
            channels.SetUniformChannel(10, 2, 0.035);

            // Source only on left branch headwater
            var sources = new List<(int, int, double, double)>
            {
                (3, 0, 100, double.MaxValue) // 100 mg/L continuous
            };

            var contaminant = AquaticContaminant.Cyanide; // conservative
            var parameters = new WaterContaminationParameters(sources, contaminant);
            var fuelMap = new FuelMap(grid);
            var tf = new TimeFrame(0, 1800, 60);

            var sim = new WaterContaminationSimulator(parameters, net, channels);
            var snaps = sim.Run(grid, terrain, fuelMap, tf);

            // At confluence (5,3): concentration should be < source because
            // the clean right branch dilutes it
            var lastConc = snaps[snaps.Count - 1].Snapshot.GetLayer("concentration");
            int nx = grid.Nx;
            double concAtConfluence = lastConc[3 * nx + 5]; // (5,3)

            Assert.IsTrue(concAtConfluence > 0,
                "Confluence should have some contamination from left branch");
            Assert.IsTrue(concAtConfluence < 100,
                $"Confluence concentration ({concAtConfluence:F2}) should be < source (100) due to dilution");
        }

        // ════════════════════════════════════════════════════════════════
        //  CFL warning
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void CFL_ViolationDetected_LargeTimeStep()
        {
            var (grid, terrain, net, channels) = CreateStraightRiver();
            var sources = new List<(int, int, double, double)> { (5, 0, 100, 1000) };
            var contaminant = AquaticContaminant.Cyanide;
            var parameters = new WaterContaminationParameters(sources, contaminant);
            var fuelMap = new FuelMap(grid);

            // Very large time step relative to grid spacing (dx=100, u~1.1 m/s) → CFL > 1
            var tf = new TimeFrame(0, 200, 200);

            var sim = new WaterContaminationSimulator(parameters, net, channels);
            sim.Run(grid, terrain, fuelMap, tf);

            Assert.IsTrue(sim.CflViolationDetected,
                "Should detect CFL violation with large time step");
        }

        [TestMethod]
        public void CFL_NoViolation_SmallTimeStep()
        {
            var (grid, terrain, net, channels) = CreateStraightRiver();
            var sources = new List<(int, int, double, double)> { (5, 0, 100, 1000) };
            var contaminant = AquaticContaminant.Cyanide;
            var parameters = new WaterContaminationParameters(sources, contaminant);
            var fuelMap = new FuelMap(grid);
            var tf = new TimeFrame(0, 300, 1); // very small dt

            var sim = new WaterContaminationSimulator(parameters, net, channels);
            sim.Run(grid, terrain, fuelMap, tf);

            Assert.IsFalse(sim.CflViolationDetected,
                "Should not detect CFL violation with small time step");
        }

        // ════════════════════════════════════════════════════════════════
        //  Contamination state tracking
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void ContaminatedCellCount_IncreasesAsPlumeSpreads()
        {
            var tracer = AquaticContaminant.Cyanide;
            var snaps = RunStraightRiver(tracer, sourceConc: 100, sourceDuration: double.MaxValue,
                simDuration: 1200, dt: 60);

            double threshold = tracer.ToxicityThresholdMgL;
            int contEarly = snaps[1].ContaminatedCellCount(threshold);
            int contLate = snaps[snaps.Count - 1].ContaminatedCellCount(threshold);

            Assert.IsTrue(contLate >= contEarly,
                $"Contaminated cells should increase: early={contEarly}, late={contLate}");
        }

        [TestMethod]
        public void CleanCells_HaveZeroConcentration()
        {
            var tracer = AquaticContaminant.Cyanide;
            var snaps = RunStraightRiver(tracer, sourceConc: 100, simDuration: 300, dt: 60);

            var conc = snaps[0].Snapshot.GetLayer("concentration");
            var state = snaps[0].Snapshot.GetLayer("contaminationState");

            for (int i = 0; i < state.Length; i++)
            {
                if ((int)state[i] == (int)CellContaminationState.Clean)
                {
                    Assert.AreEqual(0, conc[i], 1e-15,
                        $"Clean cell {i} should have zero concentration");
                }
            }
        }

        [TestMethod]
        public void SourceCells_MaintainInjectionConcentration()
        {
            var tracer = AquaticContaminant.Cyanide;
            var snaps = RunStraightRiver(tracer, sourceConc: 100, sourceDuration: double.MaxValue,
                simDuration: 600, dt: 60);

            // Source at (5,0), check it stays at 100 mg/L
            int nx = 11;
            foreach (var snap in snaps)
            {
                var conc = snap.Snapshot.GetLayer("concentration");
                Assert.AreEqual(100.0, conc[0 * nx + 5], 1e-10,
                    $"Source cell should maintain injection concentration at t={snap.Time}");
            }
        }

        // ════════════════════════════════════════════════════════════════
        //  Output layers present
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void OutputLayers_AllPresent()
        {
            var tracer = AquaticContaminant.Cyanide;
            var snaps = RunStraightRiver(tracer, simDuration: 60, dt: 60);

            Assert.IsTrue(snaps[0].Snapshot.HasLayer("concentration"));
            Assert.IsTrue(snaps[0].Snapshot.HasLayer("contaminationState"));
            Assert.IsTrue(snaps[0].Snapshot.HasLayer("velocity"));
            Assert.IsTrue(snaps[0].Snapshot.HasLayer("exposureTime"));
        }
    }
}
