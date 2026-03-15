using CSharpNumerics.Engines.GIS.Export;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Simulation;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Enums;
using CSharpNumerics.Physics.Materials;
using CSharpNumerics.Physics.Materials.Nuclear.Isotopes;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Linq;

namespace NumericTest
{
    [TestClass]
    public class RadioactivePlumeTests
    {
        private static GeoGrid SmallGrid() =>
            new GeoGrid(-200, 200, -200, 200, 0, 50, 50);

        private static TimeFrame ShortTime() =>
            new TimeFrame(0, 120, 60);

        // ═══════════════════════════════════════════════════════════════
        //  Materials factory
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void Materials_Radioisotope_LooksUpCs137()
        {
            var mat = Materials.Radioisotope("Cs137");
            Assert.AreEqual("Cs137", mat.Isotope.Name);
            Assert.IsNotNull(mat.Chain, "Cs137 should have a known decay chain");
            Assert.AreEqual(2, mat.Chain.StepCount);
        }

        [TestMethod]
        public void Materials_Radioisotope_I131()
        {
            var mat = Materials.Radioisotope("I131");
            Assert.AreEqual("I131", mat.Isotope.Name);
            Assert.IsNotNull(mat.Chain);
            Assert.AreEqual(1, mat.Chain.StepCount);
        }

        [TestMethod]
        public void Materials_Radioisotope_FromIsotopeInstance()
        {
            var mat = Materials.Radioisotope(Isotope.Sr90);
            Assert.AreEqual("Sr90", mat.Isotope.Name);
            Assert.IsNull(mat.Chain); // no known chain for Sr90
        }

        // ═══════════════════════════════════════════════════════════════
        //  GridSnapshot layers
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void GridSnapshot_SetLayer_And_GetLayer()
        {
            var grid = SmallGrid();
            var values = new double[grid.CellCount];
            var snapshot = new GridSnapshot(grid, values, 0, 0);

            var activity = new double[grid.CellCount];
            activity[0] = 42.0;
            snapshot.SetLayer("activity", activity);

            Assert.IsTrue(snapshot.HasLayer("activity"));
            Assert.IsFalse(snapshot.HasLayer("dose"));
            Assert.AreEqual(42.0, snapshot.GetLayer("activity")[0]);
        }

        [TestMethod]
        public void GridSnapshot_LayerNames()
        {
            var grid = SmallGrid();
            var snapshot = new GridSnapshot(grid, new double[grid.CellCount], 0, 0);

            snapshot.SetLayer("activity", new double[grid.CellCount]);
            snapshot.SetLayer("dose", new double[grid.CellCount]);

            var names = snapshot.LayerNames.ToList();
            Assert.AreEqual(2, names.Count);
            Assert.IsTrue(names.Contains("activity"));
            Assert.IsTrue(names.Contains("dose"));
        }

        [TestMethod]
        [ExpectedException(typeof(System.Collections.Generic.KeyNotFoundException))]
        public void GridSnapshot_GetLayer_Missing_Throws()
        {
            var grid = SmallGrid();
            var snapshot = new GridSnapshot(grid, new double[grid.CellCount], 0, 0);
            snapshot.GetLayer("nonexistent");
        }

        [TestMethod]
        public void GridSnapshot_LayerNames_EmptyByDefault()
        {
            var grid = SmallGrid();
            var snapshot = new GridSnapshot(grid, new double[grid.CellCount], 0, 0);
            Assert.IsFalse(snapshot.LayerNames.Any());
        }

        // ═══════════════════════════════════════════════════════════════
        //  PlumeSimulator with Material
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void PlumeSimulator_WithMaterial_ProducesActivityLayer()
        {
            var sim = new PlumeSimulator(
                emissionRate: 1.0, windSpeed: 5,
                windDirection: new Vector(1, 0, 0),
                stackHeight: 20,
                sourcePosition: new Vector(0, 0, 20),
                stability: StabilityClass.D);
            sim.Material = Materials.Radioisotope("Cs137");

            var grid = SmallGrid();
            var tf = ShortTime();
            var snapshots = sim.Run(grid, tf);

            Assert.IsTrue(snapshots[0].HasLayer("activity"));
            Assert.IsTrue(snapshots[0].HasLayer("dose"));

            // Where there's concentration, there should be activity
            var actLayer = snapshots[0].GetLayer("activity");
            bool anyActivity = actLayer.Any(a => a > 0);
            Assert.IsTrue(anyActivity, "Expected some cells with activity > 0");
        }

        [TestMethod]
        public void PlumeSimulator_WithoutMaterial_NoLayers()
        {
            var sim = new PlumeSimulator(
                emissionRate: 1.0, windSpeed: 5,
                windDirection: new Vector(1, 0, 0),
                stackHeight: 20,
                sourcePosition: new Vector(0, 0, 20));

            var snapshots = sim.Run(SmallGrid(), ShortTime());
            Assert.IsFalse(snapshots[0].HasLayer("activity"));
        }

        [TestMethod]
        public void PlumeSimulator_ActivityProportionalToConcentration()
        {
            var sim = new PlumeSimulator(
                emissionRate: 1.0, windSpeed: 5,
                windDirection: new Vector(1, 0, 0),
                stackHeight: 20,
                sourcePosition: new Vector(0, 0, 20));
            sim.Material = Materials.Radioisotope("Cs137");

            var grid = SmallGrid();
            var tf = new TimeFrame(0, 0, 60); // single time step at t=0
            var snapshots = sim.Run(grid, tf);

            // At t=0, activity = concentration × specificActivity (no decay correction)
            var snap = snapshots[0];
            var act = snap.GetLayer("activity");
            for (int i = 0; i < snap.Count; i++)
            {
                if (snap[i] > 0)
                {
                    double expected = snap[i] * Isotope.Cs137.SpecificActivity;
                    Assert.AreEqual(expected, act[i], expected * 0.01,
                        $"Cell {i}: conc={snap[i]:E3}, activity={act[i]:E3}");
                    break; // check at least one
                }
            }
        }

        // ═══════════════════════════════════════════════════════════════
        //  Fluent API with material
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void FluentAPI_WithMaterial_RunSingle()
        {
            var result = RiskScenario
                .ForGaussianPlume(1.0)
                .FromSource(new Vector(0, 0, 20))
                .WithWind(5, new Vector(1, 0, 0))
                .WithStability(StabilityClass.D)
                .WithMaterial(Materials.Radioisotope("Cs137"))
                .OverGrid(SmallGrid())
                .OverTime(0, 120, 60)
                .RunSingle();

            Assert.IsNotNull(result.Snapshots);
            Assert.IsTrue(result.Snapshots.Count > 0);
            Assert.IsTrue(result.Snapshots[0].HasLayer("activity"));
            Assert.IsTrue(result.Snapshots[0].HasLayer("dose"));
        }

        [TestMethod]
        public void FluentAPI_WithMaterial_RunMonteCarlo()
        {
            var result = RiskScenario
                .ForGaussianPlume(1.0)
                .FromSource(new Vector(0, 0, 20))
                .WithWind(5, new Vector(1, 0, 0))
                .WithMaterial(Materials.Radioisotope("I131"))
                .WithVariation(v => v.WindSpeed(3, 7))
                .OverGrid(SmallGrid())
                .OverTime(0, 60, 60)
                .RunMonteCarlo(5, seed: 42)
                .Build();

            // MC result should have snapshots with layers
            Assert.IsNotNull(result.MonteCarloResult);
            var firstScenarioSnaps = result.MonteCarloResult.Snapshots[0];
            Assert.IsTrue(firstScenarioSnaps[0].HasLayer("activity"));
        }

        // ═══════════════════════════════════════════════════════════════
        //  GeoJSON export with activity/dose
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void GeoJson_IncludesActivityAndDose()
        {
            var sim = new PlumeSimulator(
                emissionRate: 1.0, windSpeed: 5,
                windDirection: new Vector(1, 0, 0),
                stackHeight: 20,
                sourcePosition: new Vector(0, 0, 20));
            sim.Material = Materials.Radioisotope("Cs137");

            var grid = SmallGrid();
            var tf = new TimeFrame(0, 0, 60);
            var snapshots = sim.Run(grid, tf);

            var json = GeoJsonExporter.ToGeoJson(snapshots[0]);

            Assert.IsTrue(json.Contains("\"activity\":"), "GeoJSON should contain activity");
            Assert.IsTrue(json.Contains("\"dose\":"), "GeoJSON should contain dose");
            Assert.IsTrue(json.Contains("\"concentration\":"), "GeoJSON should contain concentration");
        }

        [TestMethod]
        public void GeoJson_NoMaterial_NoExtraProperties()
        {
            var grid = SmallGrid();
            var snapshot = new GridSnapshot(grid, new double[grid.CellCount], 0, 0);
            var json = GeoJsonExporter.ToGeoJson(snapshot);

            Assert.IsFalse(json.Contains("\"activity\":"));
            Assert.IsFalse(json.Contains("\"dose\":"));
        }
    }
}
