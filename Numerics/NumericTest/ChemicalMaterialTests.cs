using CSharpNumerics.Engines.GIS.Export;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Simulation;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Enums;
using CSharpNumerics.Physics.Materials;
using CSharpNumerics.Physics.Materials.Chemical;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Linq;

namespace NumericTest
{
    [TestClass]
    public class ChemicalMaterialTests
    {
        // ═══════════════════════════════════════════════════════════════
        //  ChemicalSubstance — static instances
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void Chlorine_HasCorrectProperties()
        {
            var cl = ChemicalSubstance.Chlorine;
            Assert.AreEqual("Cl2", cl.Formula);
            Assert.AreEqual("Chlorine", cl.Name);
            Assert.AreEqual(70.906, cl.MolarMass, 0.01);
            Assert.AreEqual(2.49, cl.VapourDensity, 0.01);
            Assert.AreEqual(10, cl.IDLH);
            Assert.AreEqual(3, cl.ERPG2);
            Assert.AreEqual(20, cl.ERPG3);
            Assert.AreEqual(PhaseAtSTP.Gas, cl.Phase);
        }

        [TestMethod]
        public void Ammonia_HasCorrectProperties()
        {
            var nh3 = ChemicalSubstance.Ammonia;
            Assert.AreEqual("NH3", nh3.Formula);
            Assert.AreEqual(17.031, nh3.MolarMass, 0.01);
            Assert.AreEqual(0.59, nh3.VapourDensity, 0.01);
            Assert.AreEqual(300, nh3.IDLH);
            Assert.AreEqual(PhaseAtSTP.Gas, nh3.Phase);
        }

        [TestMethod]
        public void HydrogenSulfide_HasCorrectProperties()
        {
            var h2s = ChemicalSubstance.HydrogenSulfide;
            Assert.AreEqual("H2S", h2s.Formula);
            Assert.AreEqual(34.081, h2s.MolarMass, 0.01);
            Assert.AreEqual(50, h2s.IDLH);
            Assert.AreEqual(1, h2s.TLV_TWA);
            Assert.AreEqual(5, h2s.TLV_STEL);
        }

        [TestMethod]
        public void Methane_HasCorrectProperties()
        {
            var ch4 = ChemicalSubstance.Methane;
            Assert.AreEqual("CH4", ch4.Formula);
            Assert.AreEqual("Methane", ch4.Name);
            Assert.AreEqual(16.043, ch4.MolarMass, 0.01);
            Assert.IsTrue(double.IsPositiveInfinity(ch4.IDLH),
                "Methane is a simple asphyxiant — no toxicity IDLH");
            Assert.AreEqual(PhaseAtSTP.Gas, ch4.Phase);
        }

        [TestMethod]
        public void Propane_HasCorrectProperties()
        {
            var c3h8 = ChemicalSubstance.Propane;
            Assert.AreEqual("C3H8", c3h8.Formula);
            Assert.AreEqual(44.096, c3h8.MolarMass, 0.01);
            Assert.AreEqual(2100, c3h8.IDLH);
            Assert.AreEqual(PhaseAtSTP.LiquefiedGas, c3h8.Phase);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Equality
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void ChemicalSubstance_Equality_ByFormula()
        {
            Assert.AreEqual(ChemicalSubstance.Chlorine, ChemicalSubstance.Chlorine);
            Assert.AreNotEqual(ChemicalSubstance.Chlorine, ChemicalSubstance.Ammonia);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Unit conversions
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void KgM3ToPpm_RoundTrip()
        {
            var cl = ChemicalSubstance.Chlorine;
            double kgm3 = 0.001; // 1 mg/m³ chlorine
            double ppm = cl.KgM3ToPpm(kgm3);
            double back = cl.PpmToKgM3(ppm);

            Assert.IsTrue(ppm > 0);
            Assert.AreEqual(kgm3, back, 1e-10, "ppm → kg/m³ round-trip should match");
        }

        [TestMethod]
        public void KgM3ToPpm_Chlorine_OrderOfMagnitude()
        {
            var cl = ChemicalSubstance.Chlorine;
            // 1 ppm Cl2 ≈ 2.95 mg/m³ = 2.95e-6 kg/m³ at 20°C
            double kgm3 = 2.95e-6;
            double ppm = cl.KgM3ToPpm(kgm3);
            Assert.AreEqual(1.0, ppm, 0.05,
                $"2.95 mg/m³ of Cl2 should be ~1 ppm, got {ppm:F3}");
        }

        // ═══════════════════════════════════════════════════════════════
        //  ChemicalLibrary
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void ChemicalLibrary_Get_CaseInsensitive()
        {
            var cl = ChemicalLibrary.Get("cl2");
            Assert.AreEqual("Cl2", cl.Formula);
        }

        [TestMethod]
        public void ChemicalLibrary_Get_AllFiveRegistered()
        {
            Assert.IsTrue(ChemicalLibrary.TryGet("Cl2", out _));
            Assert.IsTrue(ChemicalLibrary.TryGet("NH3", out _));
            Assert.IsTrue(ChemicalLibrary.TryGet("H2S", out _));
            Assert.IsTrue(ChemicalLibrary.TryGet("CH4", out _));
            Assert.IsTrue(ChemicalLibrary.TryGet("C3H8", out _));
            Assert.AreEqual(5, ChemicalLibrary.All().Count);
        }

        [TestMethod]
        [ExpectedException(typeof(System.Collections.Generic.KeyNotFoundException))]
        public void ChemicalLibrary_Get_Unknown_Throws()
        {
            ChemicalLibrary.Get("XenonTrioxide");
        }

        [TestMethod]
        public void ChemicalLibrary_Register_Custom()
        {
            var phosgene = new ChemicalSubstance(
                "COCl2", "Phosgene", "75-44-5",
                98.92, 3.4, 8.3, PhaseAtSTP.Gas,
                2, 0.5, 1.5, 5, 0.1, 0.3);
            ChemicalLibrary.Register(phosgene);

            Assert.IsTrue(ChemicalLibrary.TryGet("COCl2", out var found));
            Assert.AreEqual("Phosgene", found.Name);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Materials.Chemical() factory
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void Materials_Chemical_CreatesDescriptor()
        {
            var mat = Materials.Chemical("Cl2");
            Assert.IsTrue(mat.IsChemical);
            Assert.IsFalse(mat.IsNuclear);
            Assert.AreEqual("Cl2", mat.Substance.Value.Formula);
        }

        [TestMethod]
        public void Materials_Chemical_FromInstance()
        {
            var mat = Materials.Chemical(ChemicalSubstance.Ammonia);
            Assert.AreEqual("NH3", mat.Substance.Value.Formula);
        }

        // ═══════════════════════════════════════════════════════════════
        //  PlumeSimulator — chemical layers (ppm, toxicDose)
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void PlumeSimulator_WithChemical_ProducesPpmLayer()
        {
            var grid = new GeoGrid(-200, 200, -200, 200, 0, 50, 50);
            var tf = new TimeFrame(0, 60, 60);

            var sim = new PlumeSimulator(
                5.0, 10, new Vector(1, 0, 0), 50,
                new Vector(0, 0, 50), StabilityClass.D);
            sim.Material = Materials.Chemical("Cl2");

            var snaps = sim.Run(grid, tf);
            Assert.IsTrue(snaps[0].HasLayer("ppm"), "Should have ppm layer");
            Assert.IsTrue(snaps[0].HasLayer("toxicDose"), "Should have toxicDose layer");

            var ppmVals = snaps[0].GetLayer("ppm");
            Assert.IsTrue(ppmVals.Any(v => v > 0),
                "Some cells should have non-zero ppm");
        }

        [TestMethod]
        public void PlumeSimulator_WithChemical_PpmMatchesConversion()
        {
            var grid = new GeoGrid(-200, 200, -200, 200, 0, 50, 50);
            var tf = new TimeFrame(0, 60, 60);

            var sim = new PlumeSimulator(
                5.0, 10, new Vector(1, 0, 0), 50,
                new Vector(0, 0, 50), StabilityClass.D);
            sim.Material = Materials.Chemical("NH3");

            var snaps = sim.Run(grid, tf);
            var conc = snaps[0].GetValues();
            var ppm = snaps[0].GetLayer("ppm");

            var nh3 = ChemicalSubstance.Ammonia;
            for (int i = 0; i < conc.Length; i++)
            {
                double expected = nh3.KgM3ToPpm(conc[i]);
                Assert.AreEqual(expected, ppm[i], 1e-10,
                    $"ppm[{i}] should match conversion");
            }
        }

        [TestMethod]
        public void PlumeSimulator_WithChemical_ToxicDoseIsPpmTimesTime()
        {
            var grid = new GeoGrid(-200, 200, -200, 200, 0, 50, 50);
            var tf = new TimeFrame(0, 120, 60);

            var sim = new PlumeSimulator(
                5.0, 10, new Vector(1, 0, 0), 50,
                new Vector(0, 0, 50), StabilityClass.D);
            sim.Material = Materials.Chemical("H2S");

            var snaps = sim.Run(grid, tf);
            var ppm = snaps[0].GetLayer("ppm");
            var toxDose = snaps[0].GetLayer("toxicDose");

            for (int i = 0; i < ppm.Length; i++)
            {
                Assert.AreEqual(ppm[i] * tf.StepSeconds, toxDose[i], 1e-10,
                    "toxicDose should equal ppm × stepSeconds");
            }
        }

        // ═══════════════════════════════════════════════════════════════
        //  Fluent API — .WithMaterial(Materials.Chemical("Cl2"))
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void FluentAPI_WithChemical_RunSingle()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .WithStability(StabilityClass.D)
                .WithMaterial(Materials.Chemical("Cl2"))
                .OverGrid(new GeoGrid(-500, 500, -500, 500, 0, 100, 50))
                .OverTime(0, 120, 60)
                .RunSingle();

            Assert.IsNotNull(result.Snapshots);
            Assert.IsTrue(result.Snapshots[0].HasLayer("ppm"));
            Assert.IsTrue(result.Snapshots[0].HasLayer("toxicDose"));
        }

        [TestMethod]
        public void FluentAPI_ExposurePolygon_WithChemicalPpm()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .WithStability(StabilityClass.D)
                .WithMaterial(Materials.Chemical("Cl2"))
                .OverGrid(new GeoGrid(-500, 500, -500, 500, 0, 100, 50))
                .OverTime(0, 120, 60)
                .RunSingle();

            // IDLH polygon for Cl2 (10 ppm)
            var idlhPoly = result.GeneratePeakExposurePolygon(10, "ppm");
            Assert.IsNotNull(idlhPoly);
            Assert.AreEqual("ppm", idlhPoly.LayerName);
        }

        // ═══════════════════════════════════════════════════════════════
        //  GeoJSON export with chemical layers
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void GeoJson_IncludesPpmAndToxicDose()
        {
            var grid = new GeoGrid(-100, 100, -100, 100, 0, 50, 50);
            var tf = new TimeFrame(0, 60, 60);

            var sim = new PlumeSimulator(
                5.0, 10, new Vector(1, 0, 0), 50,
                new Vector(0, 0, 50), StabilityClass.D);
            sim.Material = Materials.Chemical("NH3");

            var snaps = sim.Run(grid, tf);
            string json = GeoJsonExporter.ToGeoJson(snaps[0]);

            Assert.IsTrue(json.Contains("\"ppm\""), "GeoJSON should contain ppm property");
            Assert.IsTrue(json.Contains("\"toxicDose\""), "GeoJSON should contain toxicDose property");
        }

        // ═══════════════════════════════════════════════════════════════
        //  ToString
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void ChemicalSubstance_ToString_ContainsName()
        {
            string s = ChemicalSubstance.Chlorine.ToString();
            Assert.IsTrue(s.Contains("Chlorine"));
            Assert.IsTrue(s.Contains("Cl2"));
            Assert.IsTrue(s.Contains("IDLH=10"));
        }

        [TestMethod]
        public void ChemicalSubstance_ToString_AsphyxiantShowsNA()
        {
            string s = ChemicalSubstance.Methane.ToString();
            Assert.IsTrue(s.Contains("IDLH=N/A"));
        }
    }
}
