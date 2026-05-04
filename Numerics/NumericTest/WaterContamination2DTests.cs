using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Spread;
using CSharpNumerics.Engines.GIS.Spread.WaterContamination2D;
using CSharpNumerics.Engines.GIS.Terrain;
using CSharpNumerics.Physics.Materials.Water;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace NumericTest;

[TestClass]
public class WaterContamination2DTests
{
    // ═══════════════════════════════════════════════════════════════
    //  Helpers
    // ═══════════════════════════════════════════════════════════════

    private static GeoGrid MakeGrid(int n = 20, double step = 10.0)
        => new GeoGrid(0, n * step, 0, n * step, 0, 0, step);

    private static AquaticContaminant NoDecay
        => new AquaticContaminant("Tracer", CSharpNumerics.Physics.Materials.Water.Enums.ContaminantType.Chemical,
                double.PositiveInfinity, 0, 0.01, 100);

    private static AquaticContaminant FastDecay
        => new AquaticContaminant("Fast", CSharpNumerics.Physics.Materials.Water.Enums.ContaminantType.Chemical,
                10.0, 0, 0.01, 100); // 10 s half-life

    private static TerrainGrid FlatTerrain(GeoGrid grid)
        => TerrainGrid.FromFunction(grid, (x, y) => 0);

    // ═══════════════════════════════════════════════════════════════
    //  Pure diffusion
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void PureDiffusion_SpreadsRadially()
    {
        var grid = MakeGrid();
        int cx = 10, cy = 10;
        var pars = new WaterContamination2DParameters(
            new[] { (cx, cy, 50.0, double.MaxValue) },
            NoDecay,
            diffusionX: 1.0);

        var sim = new WaterContamination2DSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 200, 1));

        // After 200 steps the plume should have spread beyond the source
        var last = snaps[snaps.Count - 1];
        double concSource = last.Snapshot.GetLayer("concentration")[cy * grid.Nx + cx];
        Assert.IsTrue(concSource > 0, "Source cell should have concentration.");

        // Check that neighbours have non-zero concentration
        int idxNeighbour = cy * grid.Nx + (cx + 3);
        double concNeighbour = last.Snapshot.GetLayer("concentration")[idxNeighbour];
        Assert.IsTrue(concNeighbour > 0, "Diffusion should spread to neighbours.");
    }

    [TestMethod]
    public void PureDiffusion_IsSymmetric()
    {
        var grid = MakeGrid();
        int cx = 10, cy = 10;
        var pars = new WaterContamination2DParameters(
            new[] { (cx, cy, 50.0, 1.0) }, // source active only for 1 second
            NoDecay,
            diffusionX: 2.0);

        var sim = new WaterContamination2DSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 100, 0.5));

        var conc = snaps[snaps.Count - 1].Snapshot.GetLayer("concentration");
        // Check symmetry: left == right, up == down
        double left = conc[cy * grid.Nx + (cx - 2)];
        double right = conc[cy * grid.Nx + (cx + 2)];
        double up = conc[(cy - 2) * grid.Nx + cx];
        double down = conc[(cy + 2) * grid.Nx + cx];

        Assert.AreEqual(left, right, 1e-8, "Diffusion should be symmetric in X.");
        Assert.AreEqual(up, down, 1e-8, "Diffusion should be symmetric in Y.");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Advection
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Advection_ShiftsPlumeDownstream()
    {
        var grid = MakeGrid(30, 10.0);
        int cx = 5, cy = 15;
        double vx = 0.5; // rightward 0.5 m/s

        var pars = new WaterContamination2DParameters(
            new[] { (cx, cy, 100.0, 5.0) }, // source active 5 seconds
            NoDecay,
            diffusionX: 0.1,
            velocityField: (ix, iy) => (vx, 0));

        var sim = new WaterContamination2DSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 200, 0.5));

        var conc = snaps[snaps.Count - 1].Snapshot.GetLayer("concentration");

        // Compute concentration centre-of-mass in x
        double sumCx = 0, sumC = 0;
        for (int iy = 0; iy < grid.Ny; iy++)
            for (int ix = 0; ix < grid.Nx; ix++)
            {
                double c = conc[iy * grid.Nx + ix];
                sumCx += c * ix;
                sumC += c;
            }

        double comX = sumCx / sumC;
        Assert.IsTrue(comX > cx + 1,
            $"Centre of mass ({comX:F1}) should be downstream of source ({cx}).");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Decay
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Decay_ReducesConcentration()
    {
        var grid = MakeGrid(10, 10.0);
        int cx = 5, cy = 5;

        // Source active only for 2 seconds, then observe decay
        var pars = new WaterContamination2DParameters(
            new[] { (cx, cy, 100.0, 2.0) },
            FastDecay,
            diffusionX: 0.5);

        var sim = new WaterContamination2DSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 100, 0.5));

        // Peak concentration at t=2s vs final should be lower
        double peakAtSource = 0;
        foreach (var s in snaps)
        {
            double c = s.Snapshot.GetLayer("concentration")[cy * grid.Nx + cx];
            if (c > peakAtSource) peakAtSource = c;
        }

        var final = snaps[snaps.Count - 1].Snapshot.GetLayer("concentration");
        double totalFinal = 0;
        for (int i = 0; i < final.Length; i++) totalFinal += final[i];

        Assert.IsTrue(totalFinal < peakAtSource,
            "Total final concentration should be less than peak source concentration due to decay.");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Land mask
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void LandMask_BlocksSpread()
    {
        var grid = MakeGrid(10, 10.0);
        int cx = 2, cy = 5;

        // Create a wall at ix = 5
        var mask = new bool[grid.CellCount];
        for (int iy = 0; iy < grid.Ny; iy++)
            mask[iy * grid.Nx + 5] = true;

        var pars = new WaterContamination2DParameters(
            new[] { (cx, cy, 100.0, double.MaxValue) },
            NoDecay,
            diffusionX: 1.0,
            landMask: mask);

        var sim = new WaterContamination2DSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 200, 1));

        var conc = snaps[snaps.Count - 1].Snapshot.GetLayer("concentration");

        // Verify no concentration past the wall
        for (int iy = 0; iy < grid.Ny; iy++)
            for (int ix = 6; ix < grid.Nx; ix++)
            {
                double c = conc[iy * grid.Nx + ix];
                Assert.AreEqual(0, c, 1e-12,
                    $"Cell ({ix},{iy}) beyond land wall should be zero.");
            }

        // Verify concentration on source side exists
        double sourceConc = conc[cy * grid.Nx + cx];
        Assert.IsTrue(sourceConc > 0, "Source should have concentration.");
    }

    // ═══════════════════════════════════════════════════════════════
    //  CFL detection
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void CflViolation_DetectedForHighVelocity()
    {
        var grid = MakeGrid(10, 10.0);

        // Very high velocity => should violate CFL
        var pars = new WaterContamination2DParameters(
            new[] { (5, 5, 100.0, double.MaxValue) },
            NoDecay,
            diffusionX: 0.1,
            velocityField: (ix, iy) => (100.0, 100.0)); // way too fast

        var sim = new WaterContamination2DSimulator(pars);
        sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 10, 1));

        Assert.IsTrue(sim.CflViolationDetected, "CFL violation should be detected.");
    }

    [TestMethod]
    public void CflViolation_DetectedForHighDiffusion()
    {
        var grid = MakeGrid(10, 10.0);

        // Very high diffusion => should violate stability
        var pars = new WaterContamination2DParameters(
            new[] { (5, 5, 100.0, double.MaxValue) },
            NoDecay,
            diffusionX: 100.0); // dt/(dx²) * D >> 0.5

        var sim = new WaterContamination2DSimulator(pars);
        sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 10, 1));

        Assert.IsTrue(sim.CflViolationDetected, "Diffusion stability violation should be detected.");
    }

    [TestMethod]
    public void NoCflViolation_ForSmallTimeStep()
    {
        var grid = MakeGrid(10, 10.0);

        var pars = new WaterContamination2DParameters(
            new[] { (5, 5, 100.0, double.MaxValue) },
            NoDecay,
            diffusionX: 0.5,
            velocityField: (ix, iy) => (0.1, 0.05));

        var sim = new WaterContamination2DSimulator(pars);
        sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 10, 0.5));

        Assert.IsFalse(sim.CflViolationDetected, "Should be stable with small dt.");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Source injection
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Source_InjectsFixedConcentration()
    {
        var grid = MakeGrid(10, 10.0);
        int cx = 5, cy = 5;
        double sourceConc = 75.0;

        var pars = new WaterContamination2DParameters(
            new[] { (cx, cy, sourceConc, double.MaxValue) },
            NoDecay,
            diffusionX: 0.5);

        var sim = new WaterContamination2DSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 50, 1));

        // Source cell should maintain fixed concentration at every step
        for (int i = 0; i < snaps.Count; i++)
        {
            double c = snaps[i].Snapshot.GetLayer("concentration")[cy * grid.Nx + cx];
            Assert.AreEqual(sourceConc, c, 1e-10,
                $"Source cell at step {i} should have fixed concentration.");
        }
    }

    [TestMethod]
    public void Source_StopsAfterDuration()
    {
        var grid = MakeGrid(10, 10.0);
        int cx = 5, cy = 5;
        double sourceConc = 100.0;
        double duration = 10.0; // source off after 10s

        var pars = new WaterContamination2DParameters(
            new[] { (cx, cy, sourceConc, duration) },
            NoDecay,
            diffusionX: 0.5);

        var sim = new WaterContamination2DSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 60, 1));

        // After source deactivates, concentration at source cell should drop
        var lastConc = snaps[snaps.Count - 1].Snapshot.GetLayer("concentration")[cy * grid.Nx + cx];
        Assert.IsTrue(lastConc < sourceConc,
            "Concentration at source should drop after injection ends.");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Result metrics
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Result_MaxConcentration_IsPositive()
    {
        var grid = MakeGrid(10, 10.0);

        var pars = new WaterContamination2DParameters(
            new[] { (5, 5, 100.0, double.MaxValue) },
            NoDecay,
            diffusionX: 1.0);

        var sim = new WaterContamination2DSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 20, 1));

        var result = new WaterContamination2DResult(snaps, grid, 0.01);
        Assert.IsTrue(result.MaxConcentration > 0, "Max concentration should be positive.");
        Assert.AreEqual(100.0, result.MaxConcentration, 1e-8,
            "Max concentration should equal the source injection.");
    }

    [TestMethod]
    public void Result_AffectedAreaM2_Grows()
    {
        var grid = MakeGrid(20, 10.0);

        var pars = new WaterContamination2DParameters(
            new[] { (10, 10, 100.0, double.MaxValue) },
            NoDecay,
            diffusionX: 2.0);

        var sim = new WaterContamination2DSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 100, 1));

        var earlyResult = new WaterContamination2DResult(
            new[] { snaps[5] }, grid, 0.01);
        var lateResult = new WaterContamination2DResult(snaps, grid, 0.01);

        Assert.IsTrue(lateResult.AffectedAreaM2 >= earlyResult.AffectedAreaM2,
            "Affected area should grow over time.");
    }

    [TestMethod]
    public void Result_ExceedanceDuration_IsPositive()
    {
        var grid = MakeGrid(10, 10.0);

        var pars = new WaterContamination2DParameters(
            new[] { (5, 5, 100.0, double.MaxValue) },
            NoDecay,
            diffusionX: 1.0);

        var sim = new WaterContamination2DSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 50, 1));

        var result = new WaterContamination2DResult(snaps, grid, 0.01);
        double exceedance = result.ExceedanceDurationSeconds(0.01);
        Assert.IsTrue(exceedance > 0, "Exceedance duration should be positive.");
    }

    [TestMethod]
    public void Result_PeakArrivalTime_IsNonNegative()
    {
        var grid = MakeGrid(10, 10.0);

        var pars = new WaterContamination2DParameters(
            new[] { (5, 5, 100.0, double.MaxValue) },
            NoDecay,
            diffusionX: 1.0);

        var sim = new WaterContamination2DSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 20, 1));

        var result = new WaterContamination2DResult(snaps, grid, 0.01);
        Assert.IsTrue(result.PeakArrivalTimeSeconds >= 0);
    }

    [TestMethod]
    public void Result_GenerateContaminationExtent_ReturnsPolygon()
    {
        var grid = MakeGrid(10, 10.0);

        var pars = new WaterContamination2DParameters(
            new[] { (5, 5, 100.0, double.MaxValue) },
            NoDecay,
            diffusionX: 1.0);

        var sim = new WaterContamination2DSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 50, 1));

        var result = new WaterContamination2DResult(snaps, grid, 0.01);
        var polygon = result.GenerateContaminationExtent(snaps.Count - 1, 0.01);

        Assert.IsNotNull(polygon, "Polygon should not be null.");
        Assert.IsTrue(polygon.AreaSquareMetres > 0,
            "Contamination extent area should be positive.");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Scenario builder (fluent API)
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void ScenarioBuilder_RunSingle_ProducesResult()
    {
        var grid = MakeGrid(10, 10.0);

        var result = RiskScenario
            .ForWaterContamination2D()
            .WithSource(5, 5, 100.0, double.MaxValue)
            .WithContaminant(ContaminantLibrary.Get("Benzene"))
            .WithDiffusion(1.0)
            .OverGrid(grid)
            .OverTime(0, 50, 1)
            .RunSingle();

        Assert.IsNotNull(result);
        Assert.IsTrue(result.MaxConcentration > 0);
        Assert.IsTrue(result.Snapshots.Count > 0);
    }

    [TestMethod]
    public void ScenarioBuilder_WithCurrent_AppliesAdvection()
    {
        var grid = MakeGrid(20, 10.0);

        var result = RiskScenario
            .ForWaterContamination2D()
            .WithSource(5, 10, 100.0, 5.0)
            .WithContaminant(ContaminantLibrary.Get("Benzene"))
            .WithDiffusion(0.1)
            .WithCurrent(0.5, 0)
            .OverGrid(grid)
            .OverTime(0, 100, 0.5)
            .RunSingle();

        // Plume should have been advected rightward
        var conc = result.Snapshots[result.Snapshots.Count - 1]
            .Snapshot.GetLayer("concentration");
        double sumRight = 0, sumLeft = 0;
        for (int iy = 0; iy < grid.Ny; iy++)
            for (int ix = 0; ix < grid.Nx; ix++)
            {
                double c = conc[iy * grid.Nx + ix];
                if (ix > 10) sumRight += c;
                else if (ix < 5) sumLeft += c;
            }

        Assert.IsTrue(sumRight > sumLeft,
            "More contamination should be downstream (right) of the source.");
    }

    [TestMethod]
    public void ScenarioBuilder_WithLandMask_Blocks()
    {
        var grid = MakeGrid(10, 10.0);
        var mask = new bool[grid.CellCount];
        // Wall at ix = 7
        for (int iy = 0; iy < grid.Ny; iy++)
            mask[iy * grid.Nx + 7] = true;

        var result = RiskScenario
            .ForWaterContamination2D()
            .WithSource(5, 5, 100.0, double.MaxValue)
            .WithContaminant(ContaminantLibrary.Get("Benzene"))
            .WithDiffusion(1.0)
            .WithLandMask(mask)
            .OverGrid(grid)
            .OverTime(0, 100, 1)
            .RunSingle();

        var conc = result.Snapshots[result.Snapshots.Count - 1]
            .Snapshot.GetLayer("concentration");

        for (int iy = 0; iy < grid.Ny; iy++)
        {
            double c = conc[iy * grid.Nx + 9];
            Assert.AreEqual(0, c, 1e-12,
                $"Cell (9,{iy}) past land wall should be zero.");
        }
    }

    // ═══════════════════════════════════════════════════════════════
    //  Anisotropic diffusion
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void AnisotropicDiffusion_SpreadsMoreInXThanY()
    {
        var grid = MakeGrid(20, 10.0);
        int cx = 10, cy = 10;

        var pars = new WaterContamination2DParameters(
            new[] { (cx, cy, 50.0, 1.0) },
            NoDecay,
            diffusionX: 5.0,
            diffusionY: 0.5);

        var sim = new WaterContamination2DSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 50, 0.05));

        var conc = snaps[snaps.Count - 1].Snapshot.GetLayer("concentration");

        // Concentration at (cx+3, cy) should be larger than at (cx, cy+3)
        double spreadX = conc[cy * grid.Nx + (cx + 3)];
        double spreadY = conc[(cy + 3) * grid.Nx + cx];

        Assert.IsTrue(spreadX > spreadY,
            $"Spread in X ({spreadX:E3}) should exceed spread in Y ({spreadY:E3}) for Dx >> Dy.");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Multiple sources
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void MultipleSources_BothContribute()
    {
        var grid = MakeGrid(20, 10.0);

        var pars = new WaterContamination2DParameters(
            new[] { (5, 10, 50.0, double.MaxValue), (15, 10, 50.0, double.MaxValue) },
            NoDecay,
            diffusionX: 1.0);

        var sim = new WaterContamination2DSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 100, 1));

        var conc = snaps[snaps.Count - 1].Snapshot.GetLayer("concentration");

        // Both sources should have concentration
        Assert.IsTrue(conc[10 * grid.Nx + 5] > 0, "Source 1 should contribute.");
        Assert.IsTrue(conc[10 * grid.Nx + 15] > 0, "Source 2 should contribute.");

        // Midpoint (10, 10) should have concentration from both
        double mid = conc[10 * grid.Nx + 10];
        Assert.IsTrue(mid > 0, "Midpoint between sources should have concentration.");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Contamination state layers
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Snapshots_ContainExpectedLayers()
    {
        var grid = MakeGrid(10, 10.0);

        var pars = new WaterContamination2DParameters(
            new[] { (5, 5, 100.0, double.MaxValue) },
            NoDecay,
            diffusionX: 1.0);

        var sim = new WaterContamination2DSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 10, 1));

        var gs = snaps[5].Snapshot;
        Assert.IsTrue(gs.HasLayer("concentration"), "Should have concentration layer.");
        Assert.IsTrue(gs.HasLayer("contaminationState"), "Should have contaminationState layer.");
        Assert.IsTrue(gs.HasLayer("velocity"), "Should have velocity layer.");
        Assert.IsTrue(gs.HasLayer("exposureTime"), "Should have exposureTime layer.");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Conservation (no-decay tracer in closed domain)
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void ConservativeTracer_MassConserved()
    {
        // A source injecting for 5 seconds at 100 mg/L in a single cell,
        // then pure diffusion with Neumann BCs. Total mass should remain
        // constant after source stops.
        var grid = MakeGrid(20, 10.0);
        int cx = 10, cy = 10;

        var pars = new WaterContamination2DParameters(
            new[] { (cx, cy, 100.0, 5.0) }, // 5 s injection
            NoDecay,
            diffusionX: 0.5);

        var sim = new WaterContamination2DSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 100, 0.5));

        // Find the snapshot right after source stops (t >= 5)
        int postSourceIdx = -1;
        for (int i = 0; i < snaps.Count; i++)
        {
            if (snaps[i].Time >= 6.0) { postSourceIdx = i; break; }
        }
        Assert.IsTrue(postSourceIdx > 0, "Should find a post-source snapshot.");

        double massPost = 0;
        var concPost = snaps[postSourceIdx].Snapshot.GetLayer("concentration");
        for (int i = 0; i < concPost.Length; i++) massPost += concPost[i];

        double massFinal = 0;
        var concFinal = snaps[snaps.Count - 1].Snapshot.GetLayer("concentration");
        for (int i = 0; i < concFinal.Length; i++) massFinal += concFinal[i];

        // Should be conserved within numerical tolerance
        Assert.AreEqual(massPost, massFinal, massPost * 0.05,
            "Mass should be conserved after source stops (within 5%).");
    }
}
