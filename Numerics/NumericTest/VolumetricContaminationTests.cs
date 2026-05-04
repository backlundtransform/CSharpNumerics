using CSharpNumerics.Engines.GIS.Export;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Spread.VolumetricContamination;
using CSharpNumerics.Engines.GIS.Terrain;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace NumericTest;

[TestClass]
public class VolumetricContaminationTests
{
    // ═══════════════════════════════════════════════════════════════
    //  Helpers
    // ═══════════════════════════════════════════════════════════════

    private static GeoGrid MakeGrid3D(int n = 10, int nz = 5, double step = 10.0)
        => new GeoGrid(0, n * step, 0, n * step, 0, nz * step, step);

    private static TerrainGrid FlatTerrain(GeoGrid grid)
        => TerrainGrid.FromFunction(grid, (x, y) => 0);

    // ═══════════════════════════════════════════════════════════════
    //  Phase 5 Tests
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void VolumetricContamination_PureDiffusion_SphericalSpread()
    {
        var grid = MakeGrid3D(12, 6, 10);
        int cx = 6, cy = 6, cz = 3;

        var pars = new VolumetricContaminationParameters(
            new[] { (cx, cy, cz, 100.0, double.MaxValue) },
            horizontalDiffusivity: 1.0,
            verticalDiffusivity: 1.0); // isotropic for symmetry check

        var sim = new VolumetricContaminationSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 100, 1));

        var last = snaps[snaps.Count - 1];
        var conc = last.Snapshot.GetLayer("concentration");

        // Source should have concentration
        int srcIdx = grid.Index(cx, cy, cz);
        Assert.IsTrue(conc[srcIdx] > 0, "Source should have concentration.");

        // Neighbours in all 3 axes should have some spread
        Assert.IsTrue(conc[grid.Index(cx + 2, cy, cz)] > 0, "Should spread in +x.");
        Assert.IsTrue(conc[grid.Index(cx, cy + 2, cz)] > 0, "Should spread in +y.");
        if (cz + 2 < grid.Nz)
            Assert.IsTrue(conc[grid.Index(cx, cy, cz + 2)] > 0, "Should spread in +z.");

        // With isotropic diffusion, x/y spread should be symmetric
        double xSpread = conc[grid.Index(cx + 2, cy, cz)];
        double ySpread = conc[grid.Index(cx, cy + 2, cz)];
        Assert.AreEqual(xSpread, ySpread, xSpread * 0.1,
            "Isotropic diffusion should produce symmetric x/y spread.");
    }

    [TestMethod]
    public void VolumetricContamination_VerticalStratification_LimitedDepthSpread()
    {
        // Dv << Dh: plume should spread more horizontally than vertically
        var grid = MakeGrid3D(12, 8, 10);
        int cx = 6, cy = 6, cz = 4;

        var pars = new VolumetricContaminationParameters(
            new[] { (cx, cy, cz, 50.0, 1.0) }, // short injection
            horizontalDiffusivity: 5.0,
            verticalDiffusivity: 0.05); // 100x smaller

        var sim = new VolumetricContaminationSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 100, 0.1));

        var conc = snaps[snaps.Count - 1].Snapshot.GetLayer("concentration");

        double hSpread = conc[grid.Index(cx + 3, cy, cz)];
        double vSpread = cz + 3 < grid.Nz ? conc[grid.Index(cx, cy, cz + 3)] : 0;

        Assert.IsTrue(hSpread > vSpread,
            $"Horizontal spread ({hSpread:E3}) should exceed vertical ({vSpread:E3}) when Dv << Dh.");
    }

    [TestMethod]
    public void VolumetricContamination_HorizontalCurrent_TranslatesPlume()
    {
        var grid = MakeGrid3D(20, 4, 10);
        int cx = 5, cy = 10, cz = 2;

        var pars = new VolumetricContaminationParameters(
            new[] { (cx, cy, cz, 100.0, 5.0) },
            horizontalDiffusivity: 0.1,
            verticalDiffusivity: 0.01,
            velocityField: (ix, iy, iz) => (0.5, 0, 0)); // rightward current

        var sim = new VolumetricContaminationSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 100, 0.5));

        var conc = snaps[snaps.Count - 1].Snapshot.GetLayer("concentration");

        // Compute centre of mass in x at source depth
        double sumCx = 0, sumC = 0;
        int iz = cz;
        for (int iy = 0; iy < grid.Ny; iy++)
            for (int ix = 0; ix < grid.Nx; ix++)
            {
                double c = conc[grid.Index(ix, iy, iz)];
                sumCx += c * ix;
                sumC += c;
            }

        double comX = sumCx / sumC;
        Assert.IsTrue(comX > cx + 1,
            $"Centre of mass ({comX:F1}) should be downstream of source ({cx}).");
    }

    [TestMethod]
    public void VolumetricContamination_Decay_ReducesConcentration()
    {
        var grid = MakeGrid3D(8, 4, 10);
        int cx = 4, cy = 4, cz = 2;
        double lambda = Math.Log(2) / 10.0; // 10 s half-life

        var pars = new VolumetricContaminationParameters(
            new[] { (cx, cy, cz, 100.0, 2.0) }, // short injection
            horizontalDiffusivity: 0.5,
            verticalDiffusivity: 0.1,
            decayRate: lambda);

        var sim = new VolumetricContaminationSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 60, 0.5));

        // Find peak mass and compare to final
        double peakMass = 0;
        foreach (var s in snaps)
        {
            var c = s.Snapshot.GetLayer("concentration");
            double mass = 0;
            for (int i = 0; i < c.Length; i++) mass += c[i];
            if (mass > peakMass) peakMass = mass;
        }

        var finalConc = snaps[snaps.Count - 1].Snapshot.GetLayer("concentration");
        double finalMass = 0;
        for (int i = 0; i < finalConc.Length; i++) finalMass += finalConc[i];

        Assert.IsTrue(finalMass < peakMass * 0.5,
            "Decay should significantly reduce total mass.");
    }

    [TestMethod]
    public void VolumetricContamination_LandMask_BlocksSpread()
    {
        var grid = MakeGrid3D(10, 4, 10);
        int cx = 2, cy = 5, cz = 2;

        // Wall at ix = 5 for all depths
        var mask = new bool[grid.CellCount];
        for (int iz = 0; iz < grid.Nz; iz++)
            for (int iy = 0; iy < grid.Ny; iy++)
                mask[grid.Index(5, iy, iz)] = true;

        var pars = new VolumetricContaminationParameters(
            new[] { (cx, cy, cz, 100.0, double.MaxValue) },
            horizontalDiffusivity: 1.0,
            verticalDiffusivity: 0.5,
            landMask: mask);

        var sim = new VolumetricContaminationSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 100, 1));

        var conc = snaps[snaps.Count - 1].Snapshot.GetLayer("concentration");

        // Nothing should pass the wall
        for (int iz = 0; iz < grid.Nz; iz++)
            for (int iy = 0; iy < grid.Ny; iy++)
                for (int ix = 6; ix < grid.Nx; ix++)
                {
                    double c = conc[grid.Index(ix, iy, iz)];
                    Assert.AreEqual(0, c, 1e-12,
                        $"Cell ({ix},{iy},{iz}) past land wall should be zero.");
                }
    }

    [TestMethod]
    public void VolumetricContamination_DepthProfile_ReturnsCorrectShape()
    {
        var grid = MakeGrid3D(10, 6, 10);
        int cx = 5, cy = 5, cz = 3;

        var pars = new VolumetricContaminationParameters(
            new[] { (cx, cy, cz, 100.0, double.MaxValue) },
            horizontalDiffusivity: 0.5,
            verticalDiffusivity: 0.5);

        var sim = new VolumetricContaminationSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 50, 1));

        var result = new VolumetricContaminationResult(snaps, grid);
        var profile = result.GetDepthProfile(cx, cy);

        Assert.AreEqual(grid.Nz, profile.Depths.Length, "Depth profile should have Nz entries.");
        Assert.IsTrue(profile.Concentrations[cz] > 0, "Source depth should have concentration.");
        Assert.AreEqual(grid.ZMin + cz * grid.Step, profile.MaxConcentrationDepth, 1e-8,
            "Peak should be at source depth.");
    }

    [TestMethod]
    public void VolumetricContamination_SurfaceSlice_MatchesTopLayer()
    {
        var grid = MakeGrid3D(8, 4, 10);
        // Source at surface
        int surfIz = grid.Nz - 1;

        var pars = new VolumetricContaminationParameters(
            new[] { (4, 4, surfIz, 50.0, double.MaxValue) },
            horizontalDiffusivity: 1.0,
            verticalDiffusivity: 0.1);

        var sim = new VolumetricContaminationSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 30, 1));

        var result = new VolumetricContaminationResult(snaps, grid);
        var surface = result.SurfaceSlice(snaps.Count - 1);

        // Surface slice should match concentration at iz = Nz-1
        var conc = snaps[snaps.Count - 1].Snapshot.GetLayer("concentration");
        for (int iy = 0; iy < grid.Ny; iy++)
            for (int ix = 0; ix < grid.Nx; ix++)
            {
                double expected = conc[grid.Index(ix, iy, surfIz)];
                double actual = surface[iy * grid.Nx + ix];
                Assert.AreEqual(expected, actual, 1e-12,
                    $"Surface slice at ({ix},{iy}) should match concentration.");
            }
    }

    [TestMethod]
    public void VolumetricContamination_MaxAffectedDepth_IncreasesOverTime()
    {
        var grid = MakeGrid3D(10, 8, 10);
        int cz = 4; // mid-depth source

        var pars = new VolumetricContaminationParameters(
            new[] { (5, 5, cz, 100.0, double.MaxValue) },
            horizontalDiffusivity: 0.5,
            verticalDiffusivity: 1.0); // strong vertical diffusion

        var sim = new VolumetricContaminationSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 100, 1));

        // Early result vs late result: depth penetration should grow
        var earlyResult = new VolumetricContaminationResult(
            new[] { snaps[5] }, grid, 0.01);
        var lateResult = new VolumetricContaminationResult(snaps, grid, 0.01);

        double earlyDepth = earlyResult.MaxAffectedDepth(0.01);
        double lateDepth = lateResult.MaxAffectedDepth(0.01);

        Assert.IsTrue(lateDepth >= earlyDepth,
            $"Late depth ({lateDepth}) should be >= early depth ({earlyDepth}).");
    }

    [TestMethod]
    public void VolumetricContamination_MassConservation_NeumannBCs()
    {
        var grid = MakeGrid3D(12, 6, 10);
        int cx = 6, cy = 6, cz = 3;

        var pars = new VolumetricContaminationParameters(
            new[] { (cx, cy, cz, 100.0, 5.0) }, // 5 s injection
            horizontalDiffusivity: 0.5,
            verticalDiffusivity: 0.5,
            decayRate: 0); // no decay — conservative tracer

        var sim = new VolumetricContaminationSimulator(pars);
        var snaps = sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 80, 0.5));

        // Find first post-source snapshot
        int postIdx = -1;
        for (int i = 0; i < snaps.Count; i++)
            if (snaps[i].Time >= 6.0) { postIdx = i; break; }
        Assert.IsTrue(postIdx > 0);

        double massPost = 0;
        var concPost = snaps[postIdx].Snapshot.GetLayer("concentration");
        for (int i = 0; i < concPost.Length; i++) massPost += concPost[i];

        double massFinal = 0;
        var concFinal = snaps[snaps.Count - 1].Snapshot.GetLayer("concentration");
        for (int i = 0; i < concFinal.Length; i++) massFinal += concFinal[i];

        Assert.AreEqual(massPost, massFinal, massPost * 0.05,
            "Mass should be conserved within 5% after source stops.");
    }

    [TestMethod]
    public void VolumetricContamination_CFL_Warning()
    {
        var grid = MakeGrid3D(6, 4, 10);

        var pars = new VolumetricContaminationParameters(
            new[] { (3, 3, 2, 100.0, double.MaxValue) },
            horizontalDiffusivity: 0.1,
            verticalDiffusivity: 0.01,
            velocityField: (ix, iy, iz) => (100.0, 100.0, 100.0)); // way too fast

        var sim = new VolumetricContaminationSimulator(pars);
        sim.Run(grid, FlatTerrain(grid), new FuelMap(grid),
            new TimeFrame(0, 5, 1));

        Assert.IsTrue(sim.CflViolationDetected, "CFL violation should be detected.");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Phase 6 Tests — Integration & Export
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void ScenarioBuilder_RunSingle_ReturnsResult()
    {
        var grid = MakeGrid3D(8, 4, 10);

        var result = RiskScenario
            .ForVolumetricContamination()
            .WithSource(4, 4, 2, 100.0, double.MaxValue)
            .WithDiffusivity(1.0, 0.1)
            .OverGrid(grid)
            .OverTime(0, 30, 1)
            .RunSingle();

        Assert.IsNotNull(result);
        Assert.IsTrue(result.MaxConcentration > 0);
        Assert.IsTrue(result.Snapshots.Count > 0);
    }

    [TestMethod]
    public void ScenarioBuilder_WithCurrent_AppliesAdvection()
    {
        var grid = MakeGrid3D(15, 4, 10);

        var result = RiskScenario
            .ForVolumetricContamination()
            .WithSource(3, 7, 2, 100.0, 5.0)
            .WithDiffusivity(0.1, 0.01)
            .WithCurrent(0.5, 0, 0) // rightward
            .OverGrid(grid)
            .OverTime(0, 60, 0.5)
            .RunSingle();

        var conc = result.Snapshots[result.Snapshots.Count - 1]
            .Snapshot.GetLayer("concentration");
        double sumRight = 0, sumLeft = 0;
        for (int iz = 0; iz < grid.Nz; iz++)
            for (int iy = 0; iy < grid.Ny; iy++)
                for (int ix = 0; ix < grid.Nx; ix++)
                {
                    double c = conc[grid.Index(ix, iy, iz)];
                    if (ix > 8) sumRight += c;
                    else if (ix < 3) sumLeft += c;
                }

        Assert.IsTrue(sumRight > sumLeft, "More contamination should be downstream.");
    }

    [TestMethod]
    public void GenerateContaminationExtent_ReturnsPolygon()
    {
        var grid = MakeGrid3D(10, 4, 10);

        var result = RiskScenario
            .ForVolumetricContamination()
            .WithSource(5, 5, 3, 100.0, double.MaxValue) // source at top layer
            .WithDiffusivity(1.0, 0.5)
            .OverGrid(grid)
            .OverTime(0, 30, 1)
            .RunSingle();

        var polygon = result.GenerateContaminationExtent(result.Snapshots.Count - 1, 0.01);
        Assert.IsNotNull(polygon);
        Assert.IsTrue(polygon.AreaSquareMetres > 0,
            "Surface contamination extent should have positive area.");
    }

    [TestMethod]
    public void GeoJsonExport_HorizontalSlice_ProducesValidJson()
    {
        var grid = MakeGrid3D(6, 3, 10);

        var result = RiskScenario
            .ForVolumetricContamination()
            .WithSource(3, 3, 1, 50.0, double.MaxValue)
            .WithDiffusivity(0.5, 0.1)
            .OverGrid(grid)
            .OverTime(0, 10, 1)
            .RunSingle();

        string json = GeoJsonExporter.VolumetricContaminationToGeoJson(result, 1);
        Assert.IsTrue(json.Contains("\"type\":\"FeatureCollection\""));
        Assert.IsTrue(json.Contains("\"concentration\""));
        Assert.IsTrue(json.Contains("\"depth\""));
        Assert.IsTrue(json.Contains("\"isLand\""));
    }

    [TestMethod]
    public void GeoJsonExport_DepthProfile_ProducesValidJson()
    {
        var grid = MakeGrid3D(6, 4, 10);

        var result = RiskScenario
            .ForVolumetricContamination()
            .WithSource(3, 3, 2, 50.0, double.MaxValue)
            .WithDiffusivity(0.5, 0.5)
            .OverGrid(grid)
            .OverTime(0, 20, 1)
            .RunSingle();

        string json = GeoJsonExporter.DepthProfileToGeoJson(result, 3, 3);
        Assert.IsTrue(json.Contains("\"type\":\"FeatureCollection\""));
        Assert.IsTrue(json.Contains("\"depth\""));
        Assert.IsTrue(json.Contains("\"concentration\""));
    }

    [TestMethod]
    public void DepthProfileCsv_HasCorrectFormat()
    {
        var grid = MakeGrid3D(6, 4, 10);

        var result = RiskScenario
            .ForVolumetricContamination()
            .WithSource(3, 3, 2, 50.0, double.MaxValue)
            .WithDiffusivity(0.5, 0.5)
            .OverGrid(grid)
            .OverTime(0, 20, 1)
            .RunSingle();

        string csv = DepthProfileCsvExporter.ToCsv(result, 3, 3);
        var lines = csv.Split('\n', StringSplitOptions.RemoveEmptyEntries);
        Assert.AreEqual("Depth_m,Concentration_mgL", lines[0].Trim());
        // Should have header + Nz data rows
        Assert.AreEqual(grid.Nz + 1, lines.Length, "CSV should have header + Nz rows.");
    }

    [TestMethod]
    public void PeakConcentration3D_ReturnsCorrectLength()
    {
        var grid = MakeGrid3D(6, 3, 10);

        var result = RiskScenario
            .ForVolumetricContamination()
            .WithSource(3, 3, 1, 100.0, double.MaxValue)
            .WithDiffusivity(1.0, 0.5)
            .OverGrid(grid)
            .OverTime(0, 10, 1)
            .RunSingle();

        var peak = result.PeakConcentration3D();
        Assert.AreEqual(grid.CellCount, peak.Length);

        // Source cell should have peak = 100
        int srcIdx = grid.Index(3, 3, 1);
        Assert.AreEqual(100.0, peak[srcIdx], 1e-8);
    }
}
