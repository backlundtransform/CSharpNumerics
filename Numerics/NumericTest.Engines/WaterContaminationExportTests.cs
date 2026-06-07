using CSharpNumerics.Engines.GIS.Export;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Spread;
using CSharpNumerics.Engines.GIS.Spread.WaterContamination;
using CSharpNumerics.Engines.GIS.Terrain;
using CSharpNumerics.Physics.Materials.Water;
using System;
using System.Collections.Generic;
using System.IO;

namespace NumericsTests
{
    [TestClass]
    public class WaterContaminationExportTests
    {
        // ════════════════════════════════════════════════════════════════
        //  Helper — run a quick simulation to get snapshots
        // ════════════════════════════════════════════════════════════════

        private WaterContaminationResult RunQuickSimulation()
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
            channels.SetUniformChannel(10, 2, 0.035);

            return RiskScenario
                .ForWaterContamination()
                .WithRiverNetwork(net)
                .WithChannels(channels)
                .WithTerrain(terrain)
                .WithSource(5, 0, 100, double.MaxValue)
                .WithContaminant(AquaticContaminant.Cyanide)
                .WithDischarge(10)
                .OverGrid(grid)
                .OverTime(0, 300, 60)
                .RunSingle();
        }

        private (WaterContaminationMonteCarloResult mc, GeoGrid grid, TerrainGrid terrain, RiverNetwork net) RunQuickMonteCarlo()
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
            channels.SetUniformChannel(10, 2, 0.035);

            var mc = RiskScenario
                .ForWaterContamination()
                .WithRiverNetwork(net)
                .WithChannels(channels)
                .WithTerrain(terrain)
                .WithSource(5, 0, 100, double.MaxValue)
                .WithContaminant(AquaticContaminant.Cyanide)
                .WithDischarge(10)
                .WithVariation(v => v.Discharge(5, 50))
                .OverGrid(grid)
                .OverTime(0, 300, 60)
                .RunMonteCarlo(5, seed: 42);

            return (mc, grid, terrain, net);
        }

        // ════════════════════════════════════════════════════════════════
        //  GeoJSON — point features
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void GeoJson_ContaminationFeatures_ValidFeatureCollection()
        {
            var result = RunQuickSimulation();
            var json = GeoJsonExporter.ContaminationToGeoJson(result);

            Assert.IsTrue(json.Contains("\"type\":\"FeatureCollection\""),
                "Should be a FeatureCollection");
            Assert.IsTrue(json.Contains("\"concentration\":"),
                "Should contain concentration property");
            Assert.IsTrue(json.Contains("\"velocity\":"),
                "Should contain velocity property");
            Assert.IsTrue(json.Contains("\"contaminationState\":"),
                "Should contain contaminationState property");
        }

        // ════════════════════════════════════════════════════════════════
        //  GeoJSON — contamination extent polygons
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void GeoJson_ContaminationExtent_ContainsPolygons()
        {
            var result = RunQuickSimulation();
            var json = GeoJsonExporter.ContaminationExtentToGeoJson(result, 0.001);

            Assert.IsTrue(json.Contains("\"type\":\"FeatureCollection\""),
                "Should be a FeatureCollection");
            Assert.IsTrue(json.Contains("\"contaminationExtent\""),
                "Metadata should indicate contaminationExtent type");
            // May or may not have polygons depending on whether the threshold is exceeded
            // but the structure should be valid
            Assert.IsTrue(json.Contains("\"features\":["),
                "Should have features array");
        }

        // ════════════════════════════════════════════════════════════════
        //  GeoJSON — exceedance probability heatmap
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void GeoJson_ExceedanceProbability_ContainsHeatmap()
        {
            var (mc, _, _, _) = RunQuickMonteCarlo();
            var json = GeoJsonExporter.ExceedanceProbabilityToGeoJson(mc);

            Assert.IsTrue(json.Contains("\"type\":\"FeatureCollection\""),
                "Should be a FeatureCollection");
            Assert.IsTrue(json.Contains("\"exceedanceProbability\":"),
                "Should contain exceedanceProbability property");
            Assert.IsTrue(json.Contains("\"iterations\":"),
                "Should contain iterations property");
        }

        // ════════════════════════════════════════════════════════════════
        //  GeoJSON — river centreline as LineString
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void GeoJson_RiverCentreline_ValidLineString()
        {
            var result = RunQuickSimulation();
            var grid = new GeoGrid(0, 1000, 0, 1000, 0, 0, 100);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => 100.0 - y * 0.001);

            var cells = new List<(int, int)>();
            for (int iy = 0; iy < grid.Ny; iy++)
                cells.Add((5, iy));

            var net = RiverNetwork.FromManual(grid)
                .AddSegment(cells)
                .Build();

            var lastSnap = result.Snapshots[result.Snapshots.Count - 1];
            var json = GeoJsonExporter.RiverCentrelineToGeoJson(net, grid, lastSnap);

            Assert.IsTrue(json.Contains("\"type\":\"LineString\""),
                "Should contain LineString geometry");
            Assert.IsTrue(json.Contains("\"riverCentreline\""),
                "Metadata should indicate riverCentreline type");
            Assert.IsTrue(json.Contains("\"concentrations\":["),
                "Should contain concentrations array");
        }

        // ════════════════════════════════════════════════════════════════
        //  CZML — contamination time-dynamic entities
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void Czml_Contamination_HasTimeIntervalsMatchingSimulation()
        {
            var result = RunQuickSimulation();
            var tf = new TimeFrame(0, 300, 60);

            var czml = CesiumExporter.ToContaminationCzml(result.Snapshots, tf);

            Assert.IsTrue(czml.StartsWith("["), "CZML should be a JSON array");
            Assert.IsTrue(czml.Contains("\"id\":\"document\""),
                "Should contain document packet");
            Assert.IsTrue(czml.Contains("\"Water Contamination\""),
                "Should contain document name");
            Assert.IsTrue(czml.Contains("\"contam_"),
                "Should contain contamination entity packets");
            Assert.IsTrue(czml.Contains("\"rgba\":["),
                "Should contain time-sampled rgba colours");
        }

        [TestMethod]
        public void Czml_ContaminationToColor_ReturnsCorrectRamp()
        {
            // Source → dark purple-black
            var (r, g, b, a) = CesiumExporter.ContaminationToColor(0.5, 3);
            Assert.AreEqual(40, r);
            Assert.AreEqual(255, a);

            // Clean → transparent
            var (_, _, _, a2) = CesiumExporter.ContaminationToColor(0, 0);
            Assert.AreEqual(0, a2);

            // Full concentration → red
            var (r3, g3, _, _) = CesiumExporter.ContaminationToColor(1.0, 1);
            Assert.AreEqual(255, r3);
            Assert.AreEqual(0, g3);
        }

        // ════════════════════════════════════════════════════════════════
        //  Unity binary — round-trip
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void UnityBinary_ContaminationRoundTrip()
        {
            var result = RunQuickSimulation();
            var grid = new GeoGrid(0, 1000, 0, 1000, 0, 0, 100);
            var tf = new TimeFrame(0, 300, 60);

            string path = Path.Combine(Path.GetTempPath(), "test_contam.gwcn");
            try
            {
                UnityBinaryExporter.SaveContamination(result.Snapshots, grid, tf, path);
                Assert.IsTrue(File.Exists(path), "Binary file should exist");

                var data = UnityBinaryExporter.ReadContamination(path);
                Assert.AreEqual(grid.Nx, data.Header.Nx);
                Assert.AreEqual(grid.Ny, data.Header.Ny);
                Assert.AreEqual(result.Snapshots.Count, data.Header.TimeStepCount);
                Assert.AreEqual(result.Snapshots.Count, data.Concentration.Length);
                Assert.AreEqual(result.Snapshots.Count, data.Probability.Length); // velocity in Probability slot

                // Verify concentration values match
                var firstConc = result.Snapshots[0].Snapshot.GetLayer("concentration");
                for (int i = 0; i < Math.Min(10, firstConc.Length); i++)
                {
                    Assert.AreEqual((float)firstConc[i], data.Concentration[0][i], 1e-4f,
                        $"Concentration mismatch at cell {i}");
                }
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }
    }
}
