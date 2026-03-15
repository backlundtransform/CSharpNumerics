using CSharpNumerics.Engines.GIS.Analysis;
using CSharpNumerics.Engines.GIS.Export;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Simulation;
using CSharpNumerics.ML.Clustering.Algorithms;
using CSharpNumerics.ML.Clustering.Evaluators;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Enums;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace NumericTest
{
    [TestClass]
    public class ExportTests
    {
        // ═══════════════════════════════════════════════════════════════
        //  Helper: build synthetic data
        // ═══════════════════════════════════════════════════════════════

        private static (GeoGrid grid, TimeFrame tf, MonteCarloScenarioResult mc)
            CreateSyntheticMC(int iterations = 10)
        {
            var grid = new GeoGrid(0, 20, 0, 20, 0, 0, 10); // 3×3×1 = 9 cells
            var tf = new TimeFrame(0, 60, 60);                // 2 time steps

            int cells = grid.CellCount;
            int times = tf.Count;
            int cols = cells * times;

            var data = new double[iterations, cols];
            var snapshots = new List<List<GridSnapshot>>(iterations);

            for (int i = 0; i < iterations; i++)
            {
                for (int t = 0; t < times; t++)
                {
                    int offset = t * cells;
                    for (int c = 0; c < cells; c++)
                        data[i, offset + c] = (c == 0) ? 10.0 + i : 0.1;
                }

                var iterSnaps = new List<GridSnapshot>();
                for (int t = 0; t < times; t++)
                {
                    var vals = new double[cells];
                    for (int c = 0; c < cells; c++)
                        vals[c] = data[i, t * cells + c];
                    iterSnaps.Add(new GridSnapshot(grid, vals, tf.TimeAt(t), t));
                }
                snapshots.Add(iterSnaps);
            }

            var mc = new MonteCarloScenarioResult(
                new Matrix(data), snapshots, grid, tf);

            return (grid, tf, mc);
        }

        // ═══════════════════════════════════════════════════════════════
        //  GeoJsonExporter — snapshot
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void GeoJson_Snapshot_Contains_FeatureCollection()
        {
            var grid = new GeoGrid(0, 10, 0, 10, 0, 0, 10);
            var vals = new double[] { 1.5, 2.5, 3.5, 4.5 };
            var snap = new GridSnapshot(grid, vals, 30.0, 0);

            string json = GeoJsonExporter.ToGeoJson(snap);

            Assert.IsTrue(json.Contains("\"type\":\"FeatureCollection\""));
            Assert.IsTrue(json.Contains("\"type\":\"Feature\""));
            Assert.IsTrue(json.Contains("\"type\":\"Point\""));
            Assert.IsTrue(json.Contains("\"concentration\""));
        }

        [TestMethod]
        public void GeoJson_Snapshot_Contains_All_Cells()
        {
            var grid = new GeoGrid(0, 20, 0, 20, 0, 0, 10); // 9 cells
            var vals = new double[9];
            var snap = new GridSnapshot(grid, vals, 0, 0);

            string json = GeoJsonExporter.ToGeoJson(snap);

            int featureCount = json.Split(new[] { "\"type\":\"Feature\"" }, StringSplitOptions.None).Length - 1;
            Assert.AreEqual(9, featureCount);
        }

        [TestMethod]
        public void GeoJson_Snapshot_With_Metadata()
        {
            var grid = new GeoGrid(0, 10, 0, 10, 0, 0, 10);
            var snap = new GridSnapshot(grid, new double[4], 0, 0);

            var meta = new ExportMetadata
            {
                Simulation = "GaussianPlume",
                EmissionRate = 5.0,
                Unit = "kg/m³"
            };

            string json = GeoJsonExporter.ToGeoJson(snap, meta);

            Assert.IsTrue(json.Contains("\"simulation\":\"GaussianPlume\""));
            Assert.IsTrue(json.Contains("\"emissionRate\":5"));
            Assert.IsTrue(json.Contains("\"unit\":\"kg/m\\u00B3\"") ||
                          json.Contains("\"unit\":\"kg/m³\""));
        }

        // ═══════════════════════════════════════════════════════════════
        //  GeoJsonExporter — ProbabilityMap
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void GeoJson_ProbabilityMap_Contains_Probability()
        {
            var (_, _, mc) = CreateSyntheticMC();
            var pmap = ProbabilityMap.Build(mc, timeIndex: 0, threshold: 5.0);

            string json = GeoJsonExporter.ToGeoJson(pmap);

            Assert.IsTrue(json.Contains("\"probability\""));
            Assert.IsTrue(json.Contains("\"threshold\""));
        }

        // ═══════════════════════════════════════════════════════════════
        //  GeoJsonExporter — combined snapshot + probability
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void GeoJson_Combined_Contains_Both_Layers()
        {
            var (grid, _, mc) = CreateSyntheticMC();
            var snap = mc.Snapshots[0][0];
            var pmap = ProbabilityMap.Build(mc, timeIndex: 0, threshold: 5.0);

            string json = GeoJsonExporter.ToGeoJson(snap, pmap, clusterLabel: 1);

            Assert.IsTrue(json.Contains("\"concentration\""));
            Assert.IsTrue(json.Contains("\"probability\""));
            Assert.IsTrue(json.Contains("\"cluster\":1"));
        }

        // ═══════════════════════════════════════════════════════════════
        //  GeoJsonExporter — TimeAnimator
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void GeoJson_TimeAnimator_Contains_AllTimeSteps()
        {
            var (_, _, mc) = CreateSyntheticMC();
            var anim = TimeAnimator.Build(mc, threshold: 5.0);

            string json = GeoJsonExporter.ToGeoJson(anim);

            // 9 cells × 2 time steps = 18 features
            int featureCount = json.Split(new[] { "\"type\":\"Feature\"" }, StringSplitOptions.None).Length - 1;
            Assert.AreEqual(18, featureCount);
        }

        // ═══════════════════════════════════════════════════════════════
        //  GeoJsonExporter — file I/O
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void GeoJson_Save_Creates_File()
        {
            var grid = new GeoGrid(0, 10, 0, 10, 0, 0, 10);
            var snap = new GridSnapshot(grid, new double[4], 0, 0);
            string path = Path.Combine(Path.GetTempPath(), $"test_geojson_{Guid.NewGuid()}.geojson");

            try
            {
                GeoJsonExporter.Save(snap, path);
                Assert.IsTrue(File.Exists(path));
                string content = File.ReadAllText(path);
                Assert.IsTrue(content.Contains("\"FeatureCollection\""));
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [TestMethod]
        public void GeoJson_SavePerTimeStep_Creates_Multiple_Files()
        {
            var (_, _, mc) = CreateSyntheticMC();
            var anim = TimeAnimator.Build(mc, threshold: 5.0);
            string basePath = Path.Combine(Path.GetTempPath(), $"test_pertime_{Guid.NewGuid()}");

            try
            {
                var paths = GeoJsonExporter.SavePerTimeStep(anim, basePath);
                Assert.AreEqual(2, paths.Count);
                foreach (var p in paths)
                {
                    Assert.IsTrue(File.Exists(p));
                    Assert.IsTrue(p.EndsWith(".geojson"));
                }
            }
            finally
            {
                foreach (var p in Directory.GetFiles(Path.GetTempPath(), Path.GetFileName(basePath) + "*"))
                    File.Delete(p);
            }
        }

        // ═══════════════════════════════════════════════════════════════
        //  UnityBinaryExporter — round-trip
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void UnityBinary_RoundTrip_Header()
        {
            var grid = new GeoGrid(0, 20, 0, 20, 0, 0, 10);
            var tf = new TimeFrame(0, 60, 60);
            var snap = new GridSnapshot(grid, new double[grid.CellCount], 0, 0);
            string path = Path.Combine(Path.GetTempPath(), $"test_unity_{Guid.NewGuid()}.bin");

            try
            {
                UnityBinaryExporter.Save(new[] { snap }, grid, tf, path);
                var header = UnityBinaryExporter.ReadHeader(path);

                Assert.AreEqual(grid.Nx, header.Nx);
                Assert.AreEqual(grid.Ny, header.Ny);
                Assert.AreEqual(grid.Nz, header.Nz);
                Assert.AreEqual(1, header.TimeStepCount);
                Assert.AreEqual(grid.XMin, header.XMin, 1e-10);
                Assert.AreEqual(grid.XMax, header.XMax, 1e-10);
                Assert.AreEqual(tf.Start, header.TStart, 1e-10);
                Assert.AreEqual(tf.StepSeconds, header.TStep, 1e-10);
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [TestMethod]
        public void UnityBinary_RoundTrip_Concentration()
        {
            var grid = new GeoGrid(0, 20, 0, 20, 0, 0, 10); // 9 cells
            var tf = new TimeFrame(0, 60, 60);
            var vals = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
            var snap = new GridSnapshot(grid, vals, 0, 0);
            string path = Path.Combine(Path.GetTempPath(), $"test_unity_conc_{Guid.NewGuid()}.bin");

            try
            {
                UnityBinaryExporter.Save(new[] { snap }, grid, tf, path);
                var data = UnityBinaryExporter.Read(path);

                Assert.AreEqual(1, data.Header.TimeStepCount);
                Assert.AreEqual(9, data.Concentration[0].Length);

                for (int i = 0; i < 9; i++)
                    Assert.AreEqual((float)vals[i], data.Concentration[0][i], 1e-5f);

                // Probability layer should be zeros
                for (int i = 0; i < 9; i++)
                    Assert.AreEqual(0f, data.Probability[0][i]);
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [TestMethod]
        public void UnityBinary_RoundTrip_ConcentrationAndProbability()
        {
            var (grid, tf, mc) = CreateSyntheticMC();
            var anim = TimeAnimator.Build(mc, threshold: 5.0);
            var meanSnaps = new[] { mc.Snapshots[0][0], mc.Snapshots[0][1] };
            string path = Path.Combine(Path.GetTempPath(), $"test_unity_both_{Guid.NewGuid()}.bin");

            try
            {
                UnityBinaryExporter.Save(meanSnaps, anim, path);
                var data = UnityBinaryExporter.Read(path);

                Assert.AreEqual(2, data.Header.TimeStepCount);
                Assert.AreEqual(grid.CellCount, data.Concentration[0].Length);
                Assert.AreEqual(grid.CellCount, data.Probability[0].Length);

                // Concentration at cell 0 should be 10
                Assert.AreEqual(10f, data.Concentration[0][0], 1e-3f);
                // Probability at cell 0 should be 1.0 (all scenarios exceed threshold)
                Assert.AreEqual(1.0f, data.Probability[0][0], 1e-3f);
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [TestMethod]
        public void UnityBinary_ProbabilityOnly()
        {
            var (grid, _, mc) = CreateSyntheticMC();
            var anim = TimeAnimator.Build(mc, threshold: 5.0);
            string path = Path.Combine(Path.GetTempPath(), $"test_unity_prob_{Guid.NewGuid()}.bin");

            try
            {
                UnityBinaryExporter.Save(anim, path);
                var data = UnityBinaryExporter.Read(path);

                Assert.AreEqual(2, data.Header.TimeStepCount);
                // Concentration should be zeros
                Assert.AreEqual(0f, data.Concentration[0][0]);
                // Probability at cell 0 should be 1.0
                Assert.AreEqual(1.0f, data.Probability[0][0], 1e-3f);
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [TestMethod]
        [ExpectedException(typeof(InvalidDataException))]
        public void UnityBinary_InvalidMagic_Throws()
        {
            string path = Path.Combine(Path.GetTempPath(), $"test_unity_bad_{Guid.NewGuid()}.bin");
            try
            {
                File.WriteAllBytes(path, new byte[] { 0, 0, 0, 0 });
                UnityBinaryExporter.ReadHeader(path);
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        // ═══════════════════════════════════════════════════════════════
        //  CesiumExporter — CZML
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void Cesium_ToCzml_ContainsDocumentAndEntities()
        {
            var (_, _, mc) = CreateSyntheticMC();
            var anim = TimeAnimator.Build(mc, threshold: 5.0);

            string czml = CesiumExporter.ToCzml(anim);

            Assert.IsTrue(czml.StartsWith("["));
            Assert.IsTrue(czml.EndsWith("]"));
            Assert.IsTrue(czml.Contains("\"id\":\"document\""));
            Assert.IsTrue(czml.Contains("\"id\":\"cell_0\""));
            Assert.IsTrue(czml.Contains("\"rgba\""));
        }

        [TestMethod]
        public void Cesium_ToCzml_MinProbability_Filters_Cells()
        {
            var (_, _, mc) = CreateSyntheticMC();
            var anim = TimeAnimator.Build(mc, threshold: 5.0);

            string czmlAll = CesiumExporter.ToCzml(anim, minProbability: 0);
            string czmlFiltered = CesiumExporter.ToCzml(anim, minProbability: 0.5);

            // Filtered should have fewer cell entities
            int countAll = czmlAll.Split(new[] { "\"id\":\"cell_" }, StringSplitOptions.None).Length - 1;
            int countFiltered = czmlFiltered.Split(new[] { "\"id\":\"cell_" }, StringSplitOptions.None).Length - 1;

            Assert.IsTrue(countFiltered <= countAll,
                $"Filtered ({countFiltered}) should be <= all ({countAll})");
        }

        [TestMethod]
        public void Cesium_GeoJsonCesium_Contains_CesiumMetadata()
        {
            var (_, _, mc) = CreateSyntheticMC();
            var pmap = ProbabilityMap.Build(mc, timeIndex: 0, threshold: 5.0);

            string json = CesiumExporter.ToGeoJsonCesium(pmap);

            Assert.IsTrue(json.Contains("\"cesium\""));
            Assert.IsTrue(json.Contains("\"heightReference\""));
            Assert.IsTrue(json.Contains("\"cesium:color\""));
        }

        [TestMethod]
        public void Cesium_ProbabilityToColor_Boundaries()
        {
            var (r0, g0, b0, a0) = CesiumExporter.ProbabilityToColor(0);
            Assert.AreEqual(0, a0, "Zero probability → transparent");

            var (r1, g1, b1, a1) = CesiumExporter.ProbabilityToColor(1.0);
            Assert.AreEqual(255, r1, "High probability → red");
            Assert.AreEqual(0, g1, "High probability → no green");
            Assert.AreEqual(255, a1, "High probability → opaque");

            var (rMid, gMid, bMid, aMid) = CesiumExporter.ProbabilityToColor(0.5);
            Assert.IsTrue(rMid > 0 && gMid > 0, "Mid probability → yellow-ish");
        }

        [TestMethod]
        public void Cesium_Save_Creates_File()
        {
            var (_, _, mc) = CreateSyntheticMC();
            var anim = TimeAnimator.Build(mc, threshold: 5.0);
            string path = Path.Combine(Path.GetTempPath(), $"test_cesium_{Guid.NewGuid()}.czml");

            try
            {
                CesiumExporter.Save(anim, path);
                Assert.IsTrue(File.Exists(path));
                var content = File.ReadAllText(path);
                Assert.IsTrue(content.Contains("\"document\""));
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        // ═══════════════════════════════════════════════════════════════
        //  ScenarioResult export convenience methods
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void ScenarioResult_ExportGeoJson_String()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .OverGrid(new GeoGrid(-50, 50, -50, 50, 0, 0, 50))
                .OverTime(0, 60, 60)
                .RunMonteCarlo(10, seed: 42)
                .Build(threshold: 1e-8);

            string json = result.ExportGeoJson();
            Assert.IsTrue(json.Contains("\"FeatureCollection\""));
        }

        [TestMethod]
        public void ScenarioResult_ExportGeoJson_File()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .OverGrid(new GeoGrid(-50, 50, -50, 50, 0, 0, 50))
                .OverTime(0, 60, 60)
                .RunMonteCarlo(10, seed: 42)
                .Build(threshold: 1e-8);

            string path = Path.Combine(Path.GetTempPath(), $"test_result_{Guid.NewGuid()}.geojson");
            try
            {
                result.ExportGeoJson(path);
                Assert.IsTrue(File.Exists(path));
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [TestMethod]
        public void ScenarioResult_ExportUnity_File()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .OverGrid(new GeoGrid(-50, 50, -50, 50, 0, 0, 50))
                .OverTime(0, 60, 60)
                .RunMonteCarlo(10, seed: 42)
                .Build(threshold: 1e-8);

            string path = Path.Combine(Path.GetTempPath(), $"test_result_{Guid.NewGuid()}.bin");
            try
            {
                result.ExportUnity(path);
                Assert.IsTrue(File.Exists(path));

                var data = UnityBinaryExporter.Read(path);
                Assert.AreEqual(2, data.Header.TimeStepCount);
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [TestMethod]
        public void ScenarioResult_ExportCesium_File()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .OverGrid(new GeoGrid(-50, 50, -50, 50, 0, 0, 50))
                .OverTime(0, 60, 60)
                .RunMonteCarlo(10, seed: 42)
                .Build(threshold: 1e-8);

            string path = Path.Combine(Path.GetTempPath(), $"test_result_{Guid.NewGuid()}.czml");
            try
            {
                result.ExportCesium(path);
                Assert.IsTrue(File.Exists(path));
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [TestMethod]
        public void ScenarioResult_Deterministic_ExportGeoJson()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .OverGrid(new GeoGrid(-50, 50, -50, 50, 0, 0, 50))
                .OverTime(0, 60, 60)
                .RunSingle();

            string json = result.ExportGeoJson();
            Assert.IsTrue(json.Contains("\"FeatureCollection\""));
        }

        [TestMethod]
        public void ScenarioResult_Deterministic_ExportUnity()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .OverGrid(new GeoGrid(-50, 50, -50, 50, 0, 0, 50))
                .OverTime(0, 60, 60)
                .RunSingle();

            string path = Path.Combine(Path.GetTempPath(), $"test_det_{Guid.NewGuid()}.bin");
            try
            {
                result.ExportUnity(path);
                Assert.IsTrue(File.Exists(path));
                var data = UnityBinaryExporter.Read(path);
                Assert.AreEqual(2, data.Header.TimeStepCount);
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }
    }
}
