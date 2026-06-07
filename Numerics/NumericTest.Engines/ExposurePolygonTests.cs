using CSharpNumerics.Engines.GIS.Analysis;
using CSharpNumerics.Engines.GIS.Export;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Simulation;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Environmental.Enums;
using CSharpNumerics.Physics.Materials;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Collections.Generic;
using System.Linq;

namespace NumericTest
{
    [TestClass]
    public class ExposurePolygonTests
    {
        // ═══════════════════════════════════════════════════════════════
        //  Helper factories
        // ═══════════════════════════════════════════════════════════════

        private static GeoGrid TestGrid() =>
            new GeoGrid(-500, 500, -500, 500, 0, 50, 50);

        private static TimeFrame TestTime() =>
            new TimeFrame(0, 120, 60);

        private static List<GridSnapshot> CreateSnapshotsWithHotspot(
            GeoGrid grid, TimeFrame tf, double peakConc = 0.01)
        {
            var snaps = new List<GridSnapshot>();

            for (int t = 0; t < tf.Count; t++)
            {
                double[] vals = new double[grid.CellCount];

                // Create a hotspot in the positive-x, y~0 region (ground level)
                for (int iy = 0; iy < grid.Ny; iy++)
                {
                    for (int ix = 0; ix < grid.Nx; ix++)
                    {
                        var pos = grid.CellCentre(ix, iy, 0);
                        double r = System.Math.Sqrt(pos.x * pos.x + pos.y * pos.y);
                        double concentration = peakConc * System.Math.Exp(-r * r / (200 * 200));

                        // Peak grows with time
                        double timeFactor = (t + 1.0) / tf.Count;
                        for (int iz = 0; iz < grid.Nz; iz++)
                        {
                            int idx = grid.Index(ix, iy, iz);
                            vals[idx] = concentration * timeFactor * (iz == 0 ? 1.0 : 0.1);
                        }
                    }
                }

                snaps.Add(new GridSnapshot(grid, vals, tf.TimeAt(t), t));
            }

            return snaps;
        }

        // ═══════════════════════════════════════════════════════════════
        //  PeakExposure — concentration
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void PeakExposure_ReturnsNonEmptyPolygon_WhenCellsExceedThreshold()
        {
            var grid = TestGrid();
            var tf = TestTime();
            var snaps = CreateSnapshotsWithHotspot(grid, tf);

            var polygon = ExposurePolygonGenerator.PeakExposure(snaps, 1e-4);

            Assert.IsNotNull(polygon);
            Assert.IsTrue(polygon.Boundary.Count >= 3,
                $"Expected at least 3 boundary points, got {polygon.Boundary.Count}");
            Assert.IsTrue(polygon.ExceedanceCellCount > 0);
            Assert.IsTrue(polygon.AreaSquareMetres > 0);
            Assert.AreEqual(ExposureType.Peak, polygon.Type);
            Assert.AreEqual("concentration", polygon.LayerName);
        }

        [TestMethod]
        public void PeakExposure_ReturnsEmptyPolygon_WhenNoCellsExceed()
        {
            var grid = TestGrid();
            var tf = TestTime();
            var snaps = CreateSnapshotsWithHotspot(grid, tf, peakConc: 1e-10);

            var polygon = ExposurePolygonGenerator.PeakExposure(snaps, 1.0);

            Assert.IsNotNull(polygon);
            Assert.AreEqual(0, polygon.Boundary.Count);
            Assert.AreEqual(0, polygon.ExceedanceCellCount);
            Assert.AreEqual(0, polygon.AreaSquareMetres);
        }

        [TestMethod]
        public void PeakExposure_LowerThreshold_ProducesLargerArea()
        {
            var grid = TestGrid();
            var tf = TestTime();
            var snaps = CreateSnapshotsWithHotspot(grid, tf);

            var tight = ExposurePolygonGenerator.PeakExposure(snaps, 1e-3);
            var wide = ExposurePolygonGenerator.PeakExposure(snaps, 1e-5);

            Assert.IsTrue(wide.ExceedanceCellCount >= tight.ExceedanceCellCount,
                "Lower threshold should include more cells");
            Assert.IsTrue(wide.AreaSquareMetres >= tight.AreaSquareMetres,
                "Lower threshold should cover more area");
        }

        [TestMethod]
        public void PeakExposure_UsesMaxAcrossTimeSteps()
        {
            var grid = TestGrid();
            var tf = TestTime();
            var snaps = CreateSnapshotsWithHotspot(grid, tf);

            // Peak polygon using all time steps: last step has highest values
            var polygon = ExposurePolygonGenerator.PeakExposure(snaps, 5e-3);
            Assert.IsTrue(polygon.ExceedanceCellCount > 0,
                "Peak should use the maximum across time steps");

            // Using only the first step (lower values) should produce smaller area
            var firstOnly = ExposurePolygonGenerator.PeakExposure(
                new List<GridSnapshot> { snaps[0] }, 5e-3);

            Assert.IsTrue(polygon.ExceedanceCellCount >= firstOnly.ExceedanceCellCount);
        }

        // ═══════════════════════════════════════════════════════════════
        //  IntegratedExposure — concentration
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void IntegratedExposure_ReturnsNonEmptyPolygon()
        {
            var grid = TestGrid();
            var tf = TestTime();
            var snaps = CreateSnapshotsWithHotspot(grid, tf);

            var polygon = ExposurePolygonGenerator.IntegratedExposure(
                snaps, tf.StepSeconds, threshold: 0.001);

            Assert.IsNotNull(polygon);
            Assert.IsTrue(polygon.Boundary.Count >= 3);
            Assert.IsTrue(polygon.ExceedanceCellCount > 0);
            Assert.AreEqual(ExposureType.Integrated, polygon.Type);
        }

        [TestMethod]
        public void IntegratedExposure_AccumulatesOverTime()
        {
            var grid = TestGrid();
            var tf = TestTime();
            var snaps = CreateSnapshotsWithHotspot(grid, tf);

            // Integrated with all steps should give larger area than single step
            var multiStep = ExposurePolygonGenerator.IntegratedExposure(
                snaps, tf.StepSeconds, threshold: 0.001);
            var singleStep = ExposurePolygonGenerator.IntegratedExposure(
                new List<GridSnapshot> { snaps[0] }, tf.StepSeconds, threshold: 0.001);

            Assert.IsTrue(multiStep.ExceedanceCellCount >= singleStep.ExceedanceCellCount,
                "More time steps should accumulate more exposure");
        }

        // ═══════════════════════════════════════════════════════════════
        //  Named layer (activity / dose)
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void PeakExposure_WithNamedLayer_UsesLayerValues()
        {
            var grid = TestGrid();
            var tf = TestTime();
            var snaps = CreateSnapshotsWithHotspot(grid, tf);

            // Add activity layer with higher values than concentration
            foreach (var snap in snaps)
            {
                var vals = snap.GetValues();
                var activity = new double[vals.Length];
                for (int i = 0; i < vals.Length; i++)
                    activity[i] = vals[i] * 1e6; // activity >> concentration
                snap.SetLayer("activity", activity);
            }

            var concPolygon = ExposurePolygonGenerator.PeakExposure(snaps, 1e-3);
            var actPolygon = ExposurePolygonGenerator.PeakExposure(snaps, 1e-3, "activity");

            Assert.AreEqual("concentration", concPolygon.LayerName);
            Assert.AreEqual("activity", actPolygon.LayerName);
            Assert.IsTrue(actPolygon.ExceedanceCellCount >= concPolygon.ExceedanceCellCount,
                "Activity layer (higher values) should exceed threshold in more cells");
        }

        // ═══════════════════════════════════════════════════════════════
        //  ScenarioResult integration (fluent API)
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void ScenarioResult_GeneratePeakExposurePolygon_Deterministic()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .WithStability(StabilityClass.D)
                .OverGrid(new GeoGrid(-500, 500, -500, 500, 0, 100, 50))
                .OverTime(0, 120, 60)
                .RunSingle();

            var polygon = result.GeneratePeakExposurePolygon(1e-6);

            Assert.IsNotNull(polygon);
            Assert.IsTrue(polygon.ExceedanceCellCount > 0,
                "Plume should create cells above 1e-6 kg/m³");
            Assert.IsTrue(polygon.Boundary.Count >= 3);
            Assert.AreEqual(ExposureType.Peak, polygon.Type);
        }

        [TestMethod]
        public void ScenarioResult_GenerateIntegratedExposurePolygon_Deterministic()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .WithStability(StabilityClass.D)
                .OverGrid(new GeoGrid(-500, 500, -500, 500, 0, 100, 50))
                .OverTime(0, 120, 60)
                .RunSingle();

            var polygon = result.GenerateIntegratedExposurePolygon(1e-4);

            Assert.IsNotNull(polygon);
            Assert.IsTrue(polygon.ExceedanceCellCount > 0);
            Assert.AreEqual(ExposureType.Integrated, polygon.Type);
        }

        [TestMethod]
        public void ScenarioResult_GeneratePeakExposurePolygon_WithMaterial()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .WithStability(StabilityClass.D)
                .WithMaterial(Materials.Radioisotope("Cs137"))
                .OverGrid(new GeoGrid(-500, 500, -500, 500, 0, 100, 50))
                .OverTime(0, 120, 60)
                .RunSingle();

            // Peak dose polygon
            var dosePoly = result.GeneratePeakExposurePolygon(0, "dose");
            Assert.IsNotNull(dosePoly);
            Assert.AreEqual("dose", dosePoly.LayerName);

            // Peak activity polygon
            var actPoly = result.GeneratePeakExposurePolygon(0, "activity");
            Assert.IsNotNull(actPoly);
            Assert.AreEqual("activity", actPoly.LayerName);
        }

        [TestMethod]
        public void ScenarioResult_GenerateIntegratedExposurePolygon_WithMaterial()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .WithStability(StabilityClass.D)
                .WithMaterial(Materials.Radioisotope("Cs137"))
                .OverGrid(new GeoGrid(-500, 500, -500, 500, 0, 100, 50))
                .OverTime(0, 120, 60)
                .RunSingle();

            var dosePoly = result.GenerateIntegratedExposurePolygon(0, "dose");
            Assert.IsNotNull(dosePoly);
            Assert.AreEqual("dose", dosePoly.LayerName);
            Assert.AreEqual(ExposureType.Integrated, dosePoly.Type);
        }

        // ═══════════════════════════════════════════════════════════════
        //  GeoJSON polygon export
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void GeoJsonExporter_ToGeoJson_ExposurePolygon()
        {
            var grid = TestGrid();
            var tf = TestTime();
            var snaps = CreateSnapshotsWithHotspot(grid, tf);

            var polygon = ExposurePolygonGenerator.PeakExposure(snaps, 1e-4);
            string json = GeoJsonExporter.ToGeoJson(polygon);

            Assert.IsTrue(json.Contains("\"type\":\"Feature\""));
            Assert.IsTrue(json.Contains("\"type\":\"Polygon\""));
            Assert.IsTrue(json.Contains("\"coordinates\""));
            Assert.IsTrue(json.Contains("\"threshold\""));
            Assert.IsTrue(json.Contains("\"exposureType\":\"peak\""));
            Assert.IsTrue(json.Contains("\"areaSquareMetres\""));
        }

        [TestMethod]
        public void GeoJsonExporter_ToGeoJson_EmptyPolygon_HasNullGeometry()
        {
            var grid = TestGrid();
            var tf = TestTime();
            var snaps = CreateSnapshotsWithHotspot(grid, tf, peakConc: 1e-20);

            var polygon = ExposurePolygonGenerator.PeakExposure(snaps, 1.0);
            string json = GeoJsonExporter.ToGeoJson(polygon);

            Assert.IsTrue(json.Contains("\"geometry\":null"));
        }

        [TestMethod]
        public void GeoJsonExporter_ToGeoJson_MultiplePolygons_AsFeatureCollection()
        {
            var grid = TestGrid();
            var tf = TestTime();
            var snaps = CreateSnapshotsWithHotspot(grid, tf);

            var peak = ExposurePolygonGenerator.PeakExposure(snaps, 1e-4);
            var integrated = ExposurePolygonGenerator.IntegratedExposure(
                snaps, tf.StepSeconds, 0.001);

            string json = GeoJsonExporter.ToGeoJson(
                new List<ExposurePolygon> { peak, integrated });

            Assert.IsTrue(json.Contains("\"type\":\"FeatureCollection\""));
            Assert.IsTrue(json.Contains("\"exposureType\":\"peak\""));
            Assert.IsTrue(json.Contains("\"exposureType\":\"integrated\""));
        }

        // ═══════════════════════════════════════════════════════════════
        //  Boundary geometry
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void PeakExposure_BoundaryIsClosed()
        {
            var grid = TestGrid();
            var tf = TestTime();
            var snaps = CreateSnapshotsWithHotspot(grid, tf);

            var polygon = ExposurePolygonGenerator.PeakExposure(snaps, 1e-4);

            if (polygon.Boundary.Count >= 3)
            {
                var first = polygon.Boundary[0];
                var last = polygon.Boundary[polygon.Boundary.Count - 1];
                double dist = System.Math.Sqrt(
                    (first.x - last.x) * (first.x - last.x) +
                    (first.y - last.y) * (first.y - last.y));
                Assert.IsTrue(dist < 1.0,
                    $"Boundary should be closed (first-last distance = {dist:F3})");
            }
        }

        [TestMethod]
        public void PeakExposure_BoundaryPointsAreAtGroundLevel()
        {
            var grid = TestGrid();
            var tf = TestTime();
            var snaps = CreateSnapshotsWithHotspot(grid, tf);

            var polygon = ExposurePolygonGenerator.PeakExposure(snaps, 1e-4);

            foreach (var pt in polygon.Boundary)
            {
                Assert.AreEqual(grid.ZMin, pt.z, 0.01,
                    "All boundary points should be at ground level (z = zMin)");
            }
        }

        [TestMethod]
        public void PeakExposure_AreaMatchesCellCountTimesStepSquared()
        {
            var grid = TestGrid();
            var tf = TestTime();
            var snaps = CreateSnapshotsWithHotspot(grid, tf);

            var polygon = ExposurePolygonGenerator.PeakExposure(snaps, 1e-4);

            double expectedArea = polygon.ExceedanceCellCount * grid.Step * grid.Step;
            Assert.AreEqual(expectedArea, polygon.AreaSquareMetres, 0.01);
        }
    }
}
