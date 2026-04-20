using CSharpNumerics.Engines.GIS.Export;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Simulation;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Environmental.Enums;
using CSharpNumerics.Physics.Materials;
using CSharpNumerics.Physics.Materials.Biological;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Linq;

namespace NumericTest
{
    [TestClass]
    public class BiologicalMaterialTests
    {
        [TestMethod]
        public void GenericVirus_HasExpectedProperties()
        {
            var virus = BiologicalAgent.GenericVirus;

            Assert.AreEqual("Virus", virus.Code);
            Assert.AreEqual(BiologicalAgentClass.Virus, virus.Classification);
            Assert.AreEqual(0.12, virus.TypicalDiameterMicrons, 1e-12);
            Assert.AreEqual(1e-18, virus.UnitMassKg, 1e-24);
            Assert.AreEqual(6 * 3600, virus.ViabilityHalfLifeSeconds, 1e-12);
        }

        [TestMethod]
        public void BiologicalLibrary_Get_CaseInsensitiveAlias()
        {
            var bacteria = BiologicalLibrary.Get("bacteria");
            var spores = BiologicalLibrary.Get("Spores");

            Assert.AreEqual(BiologicalAgentClass.Bacteria, bacteria.Classification);
            Assert.AreEqual(BiologicalAgentClass.Spore, spores.Classification);
            Assert.AreEqual(3, BiologicalLibrary.All().Count);
        }

        [TestMethod]
        public void BiologicalAgent_UnitConversion_RoundTrip()
        {
            var bacteria = BiologicalAgent.GenericBacteria;
            double kgPerM3 = 2.5e-12;

            double units = bacteria.KgM3ToUnitsPerM3(kgPerM3);
            double roundTrip = bacteria.UnitsPerM3ToKgM3(units);

            Assert.IsTrue(units > 0);
            Assert.AreEqual(kgPerM3, roundTrip, 1e-24);
        }

        [TestMethod]
        public void Materials_Biological_CreatesDescriptor()
        {
            var mat = Materials.Biological("virus");
            var agent = mat.BiologicalAgent ?? throw new AssertFailedException("Expected a biological agent descriptor.");

            Assert.IsTrue(mat.IsBiological);
            Assert.IsFalse(mat.IsChemical);
            Assert.IsFalse(mat.IsNuclear);
            Assert.AreEqual(BiologicalAgentClass.Virus, agent.Classification);
        }

        [TestMethod]
        public void PlumeSimulator_WithBiological_ProducesLayers()
        {
            var grid = new GeoGrid(-200, 200, -200, 200, 0, 50, 50);
            var tf = new TimeFrame(0, 60, 60);

            var sim = new PlumeSimulator(
                2.0, 8, new Vector(1, 0, 0), 50,
                new Vector(0, 0, 50), StabilityClass.D);
            sim.Material = Materials.Biological("Bacteria");

            var snaps = sim.Run(grid, tf);
            Assert.IsTrue(snaps[0].HasLayer("bioUnits"));
            Assert.IsTrue(snaps[0].HasLayer("viableBioUnits"));
            Assert.IsTrue(snaps[0].HasLayer("infectiousDose"));

            var viable = snaps[0].GetLayer("viableBioUnits");
            Assert.IsTrue(viable.Any(v => v > 0), "Some cells should have non-zero viable biological units.");
        }

        [TestMethod]
        public void PlumeSimulator_WithBiological_ViabilityDecaysOverTime()
        {
            var grid = new GeoGrid(-200, 200, -200, 200, 0, 50, 50);
            var tf = new TimeFrame(0, 3600, 3600);

            var sim = new PlumeSimulator(
                2.0, 8, new Vector(1, 0, 0), 50,
                new Vector(0, 0, 50), StabilityClass.D);
            sim.Material = Materials.Biological("Virus");

            var snaps = sim.Run(grid, tf);
            var raw0 = snaps[0].GetLayer("bioUnits");
            var raw1 = snaps[1].GetLayer("bioUnits");
            var viable0 = snaps[0].GetLayer("viableBioUnits");
            var viable1 = snaps[1].GetLayer("viableBioUnits");

            int firstPositive = Array.FindIndex(raw0, value => value > 0);
            Assert.IsTrue(firstPositive >= 0, "Expected at least one populated grid cell.");
            Assert.AreEqual(raw0[firstPositive], raw1[firstPositive], 1e-12, "Steady-state plume should keep the raw layer constant.");
            Assert.IsTrue(viable1[firstPositive] < viable0[firstPositive], "Viable load should decay with elapsed time.");
        }

        [TestMethod]
        public void PlumeSimulator_WithBiological_InfectiousDoseEqualsViableUnitsTimesTime()
        {
            var grid = new GeoGrid(-200, 200, -200, 200, 0, 50, 50);
            var tf = new TimeFrame(0, 120, 60);

            var sim = new PlumeSimulator(
                2.0, 8, new Vector(1, 0, 0), 50,
                new Vector(0, 0, 50), StabilityClass.D);
            sim.Material = Materials.Biological("Spore");

            var snaps = sim.Run(grid, tf);
            var viable = snaps[0].GetLayer("viableBioUnits");
            var dose = snaps[0].GetLayer("infectiousDose");

            for (int i = 0; i < viable.Length; i++)
            {
                Assert.AreEqual(viable[i] * tf.StepSeconds, dose[i], 1e-10);
            }
        }

        [TestMethod]
        public void FluentAPI_WithBiological_RunSingle()
        {
            var result = RiskScenario
                .ForGaussianPlume(2.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(8, new Vector(1, 0, 0))
                .WithStability(StabilityClass.D)
                .WithMaterial(Materials.Biological("Virus"))
                .OverGrid(new GeoGrid(-500, 500, -500, 500, 0, 100, 50))
                .OverTime(0, 120, 60)
                .RunSingle();

            Assert.IsTrue(result.Snapshots[0].HasLayer("bioUnits"));
            Assert.IsTrue(result.Snapshots[0].HasLayer("infectiousDose"));

            var polygon = result.GeneratePeakExposurePolygon(0, "infectiousDose");
            Assert.AreEqual("infectiousDose", polygon.LayerName);
        }

        [TestMethod]
        public void GeoJson_IncludesBiologicalLayers()
        {
            var grid = new GeoGrid(-100, 100, -100, 100, 0, 50, 50);
            var tf = new TimeFrame(0, 60, 60);

            var sim = new PlumeSimulator(
                2.0, 8, new Vector(1, 0, 0), 50,
                new Vector(0, 0, 50), StabilityClass.D);
            sim.Material = Materials.Biological("bacteria");

            var snaps = sim.Run(grid, tf);
            string json = GeoJsonExporter.ToGeoJson(snaps[0]);

            Assert.IsTrue(json.Contains("\"bioUnits\""));
            Assert.IsTrue(json.Contains("\"viableBioUnits\""));
            Assert.IsTrue(json.Contains("\"infectiousDose\""));
        }
    }
}
