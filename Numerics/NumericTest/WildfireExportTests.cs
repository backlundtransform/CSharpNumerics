using CSharpNumerics.Engines.GIS.Export;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Spread;
using CSharpNumerics.Engines.GIS.Spread.Wildfire;
using CSharpNumerics.Engines.GIS.Spread.Wildfire.Enums;
using CSharpNumerics.Engines.GIS.Terrain;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Materials.Fire.Enums;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.IO;

namespace NumericTest
{
    [TestClass]
    public class WildfireExportTests
    {
        // ═══════════════════════════════════════════════════════════════
        //  Helpers
        // ═══════════════════════════════════════════════════════════════

        private static (GeoGrid grid, IReadOnlyList<SpreadSnapshot> snapshots) RunSmallFire()
        {
            var grid = new GeoGrid(0, 200, 0, 200, 0, 0, 25);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => 0);
            var fuelMap = new FuelMap(grid);
            fuelMap.SetUniformFuel(FuelModelType.ShortGrass);
            fuelMap.SetUniformMoisture(0.05);

            var result = RiskScenario
                .ForWildfire()
                .WithTerrain(terrain)
                .WithFuel(fuelMap)
                .WithIgnition(grid.Nx / 2, grid.Ny / 2)
                .WithWind(3.0, new Vector(1, 0, 0))
                .OverGrid(grid)
                .OverTime(0, 300, 60)
                .RunSingle();

            return (grid, result.Snapshots);
        }

        private static WildfireMonteCarloResult RunSmallMC()
        {
            var grid = new GeoGrid(0, 200, 0, 200, 0, 0, 25);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => 0);
            var fuelMap = new FuelMap(grid);
            fuelMap.SetUniformFuel(FuelModelType.ShortGrass);
            fuelMap.SetUniformMoisture(0.05);

            return RiskScenario
                .ForWildfire()
                .WithTerrain(terrain)
                .WithFuel(fuelMap)
                .WithIgnition(grid.Nx / 2, grid.Ny / 2)
                .WithWind(3.0, new Vector(1, 0, 0))
                .WithVariation(v => v.WindSpeed(1, 5))
                .OverGrid(grid)
                .OverTime(0, 120, 60)
                .RunMonteCarlo(5, seed: 42);
        }

        // ═══════════════════════════════════════════════════════════════
        //  GeoJSON — point features
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void GeoJson_SpreadSnapshot_Contains_FireProperties()
        {
            var (_, snapshots) = RunSmallFire();
            var lastSnap = snapshots[snapshots.Count - 1];

            string json = GeoJsonExporter.ToGeoJson(lastSnap);

            Assert.IsTrue(json.Contains("\"type\":\"FeatureCollection\""));
            Assert.IsTrue(json.Contains("\"burnState\""));
            Assert.IsTrue(json.Contains("\"flameLength\""));
            Assert.IsTrue(json.Contains("\"rateOfSpread\""));
        }

        [TestMethod]
        public void GeoJson_SpreadSnapshot_AllCells()
        {
            var (grid, snapshots) = RunSmallFire();
            string json = GeoJsonExporter.ToGeoJson(snapshots[0]);

            int featureCount = json.Split(
                new[] { "\"type\":\"Feature\"" }, StringSplitOptions.None).Length - 1;
            Assert.AreEqual(grid.CellCount, featureCount);
        }

        [TestMethod]
        public void GeoJson_SpreadTimeSeries_AllTimeSteps()
        {
            var (grid, snapshots) = RunSmallFire();
            string json = GeoJsonExporter.ToGeoJson(snapshots);

            int featureCount = json.Split(
                new[] { "\"type\":\"Feature\"" }, StringSplitOptions.None).Length - 1;
            Assert.AreEqual(grid.CellCount * snapshots.Count, featureCount);
            Assert.IsTrue(json.Contains("\"timeIndex\""));
        }

        // ═══════════════════════════════════════════════════════════════
        //  GeoJSON — fire perimeters
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void GeoJson_FirePerimeters_ValidPolygons()
        {
            var (_, snapshots) = RunSmallFire();
            var result = new WildfireScenarioResult(snapshots, snapshots[0].Grid);
            var perimeters = result.FirePerimeters;

            string json = GeoJsonExporter.FirePerimetersToGeoJson(perimeters, snapshots);

            Assert.IsTrue(json.Contains("\"type\":\"FeatureCollection\""));
            Assert.IsTrue(json.Contains("\"type\":\"Polygon\""));
            Assert.IsTrue(json.Contains("\"firePerimeters\""));
            Assert.IsTrue(json.Contains("\"timeIndex\""));
            Assert.IsTrue(json.Contains("\"areaSquareMetres\""));
        }

        // ═══════════════════════════════════════════════════════════════
        //  GeoJSON — burn probability heatmap
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void GeoJson_BurnProbability_HeatMap()
        {
            var mcResult = RunSmallMC();

            string json = GeoJsonExporter.BurnProbabilityToGeoJson(mcResult);

            Assert.IsTrue(json.Contains("\"type\":\"FeatureCollection\""));
            Assert.IsTrue(json.Contains("\"burnProbability\""));
            Assert.IsTrue(json.Contains("\"iterations\""));
        }

        // ═══════════════════════════════════════════════════════════════
        //  CZML — fire spread
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void Cesium_FireCzml_ContainsDocumentAndEntities()
        {
            var (_, snapshots) = RunSmallFire();
            var tf = new TimeFrame(0, 300, 60);

            string czml = CesiumExporter.ToFireCzml(snapshots, tf);

            Assert.IsTrue(czml.StartsWith("["));
            Assert.IsTrue(czml.EndsWith("]"));
            Assert.IsTrue(czml.Contains("\"id\":\"document\""));
            Assert.IsTrue(czml.Contains("\"id\":\"fire_"));
            Assert.IsTrue(czml.Contains("\"rgba\""));
            Assert.IsTrue(czml.Contains("\"Wildfire Spread\""));
        }

        [TestMethod]
        public void Cesium_FireCzml_TimeIntervals_MatchSimulation()
        {
            var (_, snapshots) = RunSmallFire();
            var tf = new TimeFrame(0, 300, 60);

            string czml = CesiumExporter.ToFireCzml(snapshots, tf);

            // The document clock should reference start/end
            Assert.IsTrue(czml.Contains("\"clock\""));
            Assert.IsTrue(czml.Contains("\"interval\""));
        }

        [TestMethod]
        public void Cesium_BurnStateToColor_Burning_IsRed()
        {
            var (r, g, b, a) = CesiumExporter.BurnStateToColor((int)CellBurnState.Burning);
            Assert.AreEqual(255, r, "Burning should be red.");
            Assert.IsTrue(a > 100, "Burning should be visible.");
        }

        [TestMethod]
        public void Cesium_BurnStateToColor_Burned_IsGrey()
        {
            var (r, g, b, a) = CesiumExporter.BurnStateToColor((int)CellBurnState.Burned);
            Assert.AreEqual(128, r);
            Assert.AreEqual(128, g);
            Assert.AreEqual(128, b);
        }

        [TestMethod]
        public void Cesium_BurnStateToColor_Unburned_IsTransparent()
        {
            var (_, _, _, a) = CesiumExporter.BurnStateToColor((int)CellBurnState.Unburned);
            Assert.AreEqual(0, a);
        }

        [TestMethod]
        public void Cesium_BurnStateToColor_Firebreak_IsBlue()
        {
            var (r, g, b, a) = CesiumExporter.BurnStateToColor((int)CellBurnState.Firebreak);
            Assert.IsTrue(b > r, "Firebreak should be predominantly blue.");
            Assert.IsTrue(a > 0, "Firebreak should be visible.");
        }

        // ═══════════════════════════════════════════════════════════════
        //  Unity binary — fire round-trip
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void UnityBinary_Fire_RoundTrip()
        {
            var (grid, snapshots) = RunSmallFire();
            var tf = new TimeFrame(0, 300, 60);
            string path = Path.Combine(Path.GetTempPath(), $"test_fire_{Guid.NewGuid()}.bin");

            try
            {
                UnityBinaryExporter.SaveFire(snapshots, grid, tf, path);
                Assert.IsTrue(File.Exists(path));

                var data = UnityBinaryExporter.ReadFire(path);

                Assert.AreEqual(grid.Nx, data.Header.Nx);
                Assert.AreEqual(grid.Ny, data.Header.Ny);
                Assert.AreEqual(snapshots.Count, data.Header.TimeStepCount);

                // Last time step should have some burned cells (burnState > 0)
                var lastBurn = data.Concentration[data.Header.TimeStepCount - 1]; // burnState slot
                bool hasBurned = false;
                for (int i = 0; i < lastBurn.Length; i++)
                    if (lastBurn[i] > 0) hasBurned = true;
                Assert.IsTrue(hasBurned, "Last snapshot should have burned cells.");

                // Flame length should be non-negative
                var lastFlame = data.Probability[data.Header.TimeStepCount - 1]; // flameLength slot
                for (int i = 0; i < lastFlame.Length; i++)
                    Assert.IsTrue(lastFlame[i] >= 0, $"Flame length at {i} should be >= 0.");
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [TestMethod]
        [ExpectedException(typeof(InvalidDataException))]
        public void UnityBinary_ReadFire_WrongMagic_Throws()
        {
            string path = Path.Combine(Path.GetTempPath(), $"test_fire_bad_{Guid.NewGuid()}.bin");
            try
            {
                // Write a plume file (GPLM) and try to read as fire (GFIR)
                File.WriteAllBytes(path, new byte[] { (byte)'G', (byte)'P', (byte)'L', (byte)'M' });
                UnityBinaryExporter.ReadFire(path);
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        // ═══════════════════════════════════════════════════════════════
        //  File I/O — GeoJSON
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void GeoJson_Save_SpreadSnapshot_Creates_File()
        {
            var (_, snapshots) = RunSmallFire();
            string path = Path.Combine(Path.GetTempPath(), $"test_fire_geojson_{Guid.NewGuid()}.geojson");

            try
            {
                GeoJsonExporter.Save(snapshots[0], path);
                Assert.IsTrue(File.Exists(path));
                string content = File.ReadAllText(path);
                Assert.IsTrue(content.Contains("\"burnState\""));
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [TestMethod]
        public void Cesium_SaveFireCzml_Creates_File()
        {
            var (_, snapshots) = RunSmallFire();
            var tf = new TimeFrame(0, 300, 60);
            string path = Path.Combine(Path.GetTempPath(), $"test_fire_czml_{Guid.NewGuid()}.czml");

            try
            {
                CesiumExporter.SaveFireCzml(snapshots, tf, path);
                Assert.IsTrue(File.Exists(path));
                string content = File.ReadAllText(path);
                Assert.IsTrue(content.Contains("\"document\""));
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        // ═══════════════════════════════════════════════════════════════
        //  Firebreak mask and isFirebreak export
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void FirebreakMask_IdentifiesNoFuelCells()
        {
            var grid = new GeoGrid(0, 200, 0, 200, 0, 0, 25);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => 0);
            var fuelMap = new FuelMap(grid);
            fuelMap.SetUniformFuel(FuelModelType.ShortGrass);
            fuelMap.SetUniformMoisture(0.05);

            // Set a row of NoFuel cells (water)
            for (int ix = 0; ix < grid.Nx; ix++)
                fuelMap.SetFuel(ix, 0, FuelModelType.NoFuel);

            var result = RiskScenario
                .ForWildfire()
                .WithTerrain(terrain)
                .WithFuel(fuelMap)
                .WithIgnition(grid.Nx / 2, grid.Ny / 2)
                .WithWind(3.0, new Vector(1, 0, 0))
                .OverGrid(grid)
                .OverTime(0, 120, 60)
                .RunSingle();

            var mask = result.FirebreakMask;
            Assert.AreEqual(grid.Nx * grid.Ny, mask.Length);

            // Row 0 should be firebreak
            for (int ix = 0; ix < grid.Nx; ix++)
                Assert.IsTrue(mask[0 * grid.Nx + ix],
                    $"Cell ({ix}, 0) should be firebreak.");

            // Row 1+ should NOT be firebreak (they have fuel)
            Assert.IsFalse(mask[1 * grid.Nx + 0],
                "Cell (0, 1) should not be firebreak.");
        }

        [TestMethod]
        public void MonteCarloFirebreakMask_IdentifiesNoFuelCells()
        {
            var grid = new GeoGrid(0, 200, 0, 200, 0, 0, 25);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => 0);
            var fuelMap = new FuelMap(grid);
            fuelMap.SetUniformFuel(FuelModelType.ShortGrass);
            fuelMap.SetUniformMoisture(0.05);

            // Set a column of NoFuel cells
            for (int iy = 0; iy < grid.Ny; iy++)
                fuelMap.SetFuel(0, iy, FuelModelType.NoFuel);

            var mcResult = RiskScenario
                .ForWildfire()
                .WithTerrain(terrain)
                .WithFuel(fuelMap)
                .WithIgnition(grid.Nx / 2, grid.Ny / 2)
                .WithWind(3.0, new Vector(1, 0, 0))
                .WithVariation(v => v.WindSpeed(1, 5))
                .OverGrid(grid)
                .OverTime(0, 120, 60)
                .RunMonteCarlo(3, seed: 42);

            var mask = mcResult.FirebreakMask;

            // Column 0 should be firebreak
            for (int iy = 0; iy < grid.Ny; iy++)
                Assert.IsTrue(mask[iy * grid.Nx + 0],
                    $"Cell (0, {iy}) should be firebreak.");

            // BurnProbability should be 0 for firebreak cells
            for (int iy = 0; iy < grid.Ny; iy++)
                Assert.AreEqual(0, mcResult.BurnProbability[iy * grid.Nx + 0], 1e-10,
                    $"BurnProbability at firebreak ({iy}) should be 0.");
        }

        [TestMethod]
        public void GeoJson_SpreadSnapshot_Contains_IsFirebreak()
        {
            var grid = new GeoGrid(0, 200, 0, 200, 0, 0, 25);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => 0);
            var fuelMap = new FuelMap(grid);
            fuelMap.SetUniformFuel(FuelModelType.ShortGrass);
            fuelMap.SetUniformMoisture(0.05);
            fuelMap.SetFuel(0, 0, FuelModelType.NoFuel);

            var result = RiskScenario
                .ForWildfire()
                .WithTerrain(terrain)
                .WithFuel(fuelMap)
                .WithIgnition(grid.Nx / 2, grid.Ny / 2)
                .WithWind(3.0, new Vector(1, 0, 0))
                .OverGrid(grid)
                .OverTime(0, 120, 60)
                .RunSingle();

            string json = GeoJsonExporter.ToGeoJson(result.Snapshots[0]);

            Assert.IsTrue(json.Contains("\"isFirebreak\":true"),
                "GeoJSON should contain isFirebreak:true for NoFuel cells.");
            Assert.IsTrue(json.Contains("\"isFirebreak\":false"),
                "GeoJSON should contain isFirebreak:false for normal cells.");
        }

        [TestMethod]
        public void GeoJson_BurnProbability_Contains_IsFirebreak()
        {
            var grid = new GeoGrid(0, 200, 0, 200, 0, 0, 25);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => 0);
            var fuelMap = new FuelMap(grid);
            fuelMap.SetUniformFuel(FuelModelType.ShortGrass);
            fuelMap.SetUniformMoisture(0.05);
            fuelMap.SetFuel(0, 0, FuelModelType.NoFuel);

            var mcResult = RiskScenario
                .ForWildfire()
                .WithTerrain(terrain)
                .WithFuel(fuelMap)
                .WithIgnition(grid.Nx / 2, grid.Ny / 2)
                .WithWind(3.0, new Vector(1, 0, 0))
                .WithVariation(v => v.WindSpeed(1, 5))
                .OverGrid(grid)
                .OverTime(0, 120, 60)
                .RunMonteCarlo(3, seed: 42);

            string json = GeoJsonExporter.BurnProbabilityToGeoJson(mcResult);

            Assert.IsTrue(json.Contains("\"isFirebreak\":true"),
                "Burn probability GeoJSON should flag firebreak cells.");
            Assert.IsTrue(json.Contains("\"isFirebreak\":false"),
                "Burn probability GeoJSON should flag non-firebreak cells.");
        }
    }
}
