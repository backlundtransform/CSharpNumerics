using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Spread.Wildfire;
using CSharpNumerics.Engines.GIS.Spread.Wildfire.Enums;
using CSharpNumerics.Engines.GIS.Terrain;
using CSharpNumerics.ML.Clustering.Algorithms;
using CSharpNumerics.ML.Clustering.Evaluators;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Materials.Fire.Enums;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Linq;

namespace NumericTest
{
    [TestClass]
    public class WildfireScenarioTests
    {
        // ═══════════════════════════════════════════════════════════════
        //  Helpers
        // ═══════════════════════════════════════════════════════════════

        private static GeoGrid MakeGrid(double size, double step) =>
            new GeoGrid(0, size, 0, size, 0, 0, step);

        private static TerrainGrid FlatTerrain(GeoGrid grid) =>
            TerrainGrid.FromFunction(grid, (x, y) => 0);

        private static FuelMap UniformGrass(GeoGrid grid, double moisture = 0.05)
        {
            var map = new FuelMap(grid);
            map.SetUniformFuel(FuelModelType.ShortGrass);
            map.SetUniformMoisture(moisture);
            return map;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Fluent API — deterministic round-trip
        // ═══════════════════════════════════════════════════════════════

        #region Fluent API deterministic

        [TestMethod]
        public void FluentAPI_RunSingle_Deterministic()
        {
            var grid = MakeGrid(500, 25);
            var terrain = FlatTerrain(grid);
            var fuelMap = UniformGrass(grid);

            var result = RiskScenario
                .ForWildfire()
                .WithTerrain(terrain)
                .WithFuel(fuelMap)
                .WithIgnition(grid.Nx / 2, grid.Ny / 2)
                .WithWind(3.0, new Vector(1, 0, 0))
                .WithMoisture(0.05)
                .OverGrid(grid)
                .OverTime(0, 600, 60)
                .RunSingle();

            Assert.IsNotNull(result);
            Assert.IsTrue(result.Snapshots.Count > 0, "Should produce snapshots.");
            Assert.IsTrue(result.FinalBurnedArea > 0, "Should have burned area.");
            Assert.IsTrue(result.MaxFlameLength > 0, "Should have positive flame length.");
        }

        [TestMethod]
        public void FluentAPI_RunSingle_NoMoistureOverride()
        {
            var grid = MakeGrid(300, 25);
            var terrain = FlatTerrain(grid);
            var fuelMap = UniformGrass(grid);

            var result = RiskScenario
                .ForWildfire()
                .WithTerrain(terrain)
                .WithFuel(fuelMap)
                .WithIgnition(grid.Nx / 2, grid.Ny / 2)
                .OverGrid(grid)
                .OverTime(0, 300, 60)
                .RunSingle();

            Assert.IsNotNull(result);
            Assert.IsTrue(result.Snapshots.Count > 0);
        }

        [TestMethod]
        public void FluentAPI_FirePerimeter_ValidPolygon()
        {
            var grid = MakeGrid(500, 25);
            var terrain = FlatTerrain(grid);
            var fuelMap = UniformGrass(grid);

            var result = RiskScenario
                .ForWildfire()
                .WithTerrain(terrain)
                .WithFuel(fuelMap)
                .WithIgnition(grid.Nx / 2, grid.Ny / 2)
                .WithWind(2.0, new Vector(1, 0, 0))
                .OverGrid(grid)
                .OverTime(0, 600, 60)
                .RunSingle();

            // Last time step should have a meaningful perimeter
            var lastIdx = result.Snapshots.Count - 1;
            var perimeter = result.GenerateFirePerimeter(lastIdx);

            Assert.IsNotNull(perimeter, "Perimeter should not be null.");
            Assert.IsTrue(perimeter.Boundary.Count > 0,
                "Perimeter boundary should have vertices.");
            Assert.IsTrue(perimeter.AreaSquareMetres > 0,
                "Perimeter area should be positive.");
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Monte Carlo
        // ═══════════════════════════════════════════════════════════════

        #region Monte Carlo

        [TestMethod]
        public void MonteCarlo_BurnProbability_InRange()
        {
            var grid = MakeGrid(300, 25);
            var terrain = FlatTerrain(grid);
            var fuelMap = UniformGrass(grid);

            var mcResult = RiskScenario
                .ForWildfire()
                .WithTerrain(terrain)
                .WithFuel(fuelMap)
                .WithIgnition(grid.Nx / 2, grid.Ny / 2)
                .WithWind(3.0, new Vector(1, 0, 0))
                .WithVariation(v => v
                    .WindSpeed(2, 6)
                    .Moisture(0.04, 0.10))
                .OverGrid(grid)
                .OverTime(0, 300, 60)
                .RunMonteCarlo(20, seed: 42);

            Assert.AreEqual(20, mcResult.Iterations);
            Assert.AreEqual(grid.Nx * grid.Ny, mcResult.BurnProbability.Length);

            // All probabilities should be in [0, 1]
            for (int i = 0; i < mcResult.BurnProbability.Length; i++)
            {
                Assert.IsTrue(mcResult.BurnProbability[i] >= 0 && mcResult.BurnProbability[i] <= 1,
                    $"Burn probability at cell {i} = {mcResult.BurnProbability[i]} out of range.");
            }

            // Ignition cell should have probability 1.0 (always burns)
            int cx = grid.Nx / 2;
            int cy = grid.Ny / 2;
            Assert.AreEqual(1.0, mcResult.BurnProbability[cy * grid.Nx + cx], 1e-10,
                "Ignition cell should always burn.");
        }

        [TestMethod]
        public void MonteCarlo_BurnedAreas_Reasonable()
        {
            var grid = MakeGrid(300, 25);
            var terrain = FlatTerrain(grid);
            var fuelMap = UniformGrass(grid);

            var mcResult = RiskScenario
                .ForWildfire()
                .WithTerrain(terrain)
                .WithFuel(fuelMap)
                .WithIgnition(grid.Nx / 2, grid.Ny / 2)
                .WithWind(3.0, new Vector(1, 0, 0))
                .WithVariation(v => v.WindSpeed(1, 5))
                .OverGrid(grid)
                .OverTime(0, 300, 60)
                .RunMonteCarlo(10, seed: 99);

            Assert.AreEqual(10, mcResult.BurnedAreas.Length);
            Assert.IsTrue(mcResult.MeanBurnedArea > 0,
                "Mean burned area should be positive.");
            Assert.IsTrue(mcResult.MaxBurnedArea >= mcResult.MeanBurnedArea,
                "Max should be >= mean.");
        }

        [TestMethod]
        public void MonteCarlo_WindDirectionJitter()
        {
            var grid = MakeGrid(300, 25);
            var terrain = FlatTerrain(grid);
            var fuelMap = UniformGrass(grid);

            var mcResult = RiskScenario
                .ForWildfire()
                .WithTerrain(terrain)
                .WithFuel(fuelMap)
                .WithIgnition(grid.Nx / 2, grid.Ny / 2)
                .WithWind(4.0, new Vector(1, 0, 0))
                .WithVariation(v => v.WindDirectionJitter(30))
                .OverGrid(grid)
                .OverTime(0, 300, 60)
                .RunMonteCarlo(10, seed: 7);

            Assert.AreEqual(10, mcResult.Iterations);
            // With direction jitter, different iterations should produce some variation
            Assert.IsTrue(mcResult.AllSnapshots.Count == 10);
        }

        [TestMethod]
        public void MonteCarlo_IgnitionOffset()
        {
            var grid = MakeGrid(300, 25);
            var terrain = FlatTerrain(grid);
            var fuelMap = UniformGrass(grid);

            var mcResult = RiskScenario
                .ForWildfire()
                .WithTerrain(terrain)
                .WithFuel(fuelMap)
                .WithIgnition(grid.Nx / 2, grid.Ny / 2)
                .WithWind(2.0, new Vector(1, 0, 0))
                .WithVariation(v => v.IgnitionOffset(2))
                .OverGrid(grid)
                .OverTime(0, 300, 60)
                .RunMonteCarlo(10, seed: 13);

            Assert.AreEqual(10, mcResult.Iterations);
            // Different ignition locations should produce some spread in burn prob
            int nonZeroCount = 0;
            for (int i = 0; i < mcResult.BurnProbability.Length; i++)
                if (mcResult.BurnProbability[i] > 0) nonZeroCount++;

            Assert.IsTrue(nonZeroCount > 1, "Multiple cells should have non-zero burn probability.");
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Validation errors
        // ═══════════════════════════════════════════════════════════════

        #region Validation

        [TestMethod]
        [ExpectedException(typeof(InvalidOperationException))]
        public void FluentAPI_MissingGrid_Throws()
        {
            var grid = MakeGrid(100, 25);
            RiskScenario
                .ForWildfire()
                .WithTerrain(FlatTerrain(grid))
                .WithFuel(UniformGrass(grid))
                .WithIgnition(1, 1)
                .OverTime(0, 60, 60)
                .RunSingle();
        }

        [TestMethod]
        [ExpectedException(typeof(InvalidOperationException))]
        public void FluentAPI_MissingIgnition_Throws()
        {
            var grid = MakeGrid(100, 25);
            RiskScenario
                .ForWildfire()
                .WithTerrain(FlatTerrain(grid))
                .WithFuel(UniformGrass(grid))
                .OverGrid(grid)
                .OverTime(0, 60, 60)
                .RunSingle();
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Clustering integration — AnalyzeWith / Build
        // ═══════════════════════════════════════════════════════════════

        #region Clustering

        [TestMethod]
        public void AnalyzeWith_Clustering_IdentifiesRegimes()
        {
            // Create scenarios with wind variation — high wind (8-10)
            // vs low wind (1-2) should produce distinguishable clusters.
            var grid = MakeGrid(300, 25);
            var terrain = FlatTerrain(grid);
            var fuelMap = UniformGrass(grid);

            var mcResult = RiskScenario
                .ForWildfire()
                .WithTerrain(terrain)
                .WithFuel(fuelMap)
                .WithIgnition(grid.Nx / 2, grid.Ny / 2)
                .WithWind(5.0, new Vector(1, 0, 0))
                .WithVariation(v => v.WindSpeed(1, 10).Moisture(0.04, 0.12))
                .OverGrid(grid)
                .OverTime(0, 300, 60)
                .RunMonteCarlo(30, seed: 42);

            var analysis = mcResult.AnalyzeWith(
                new KMeans { Seed = 42 },
                new SilhouetteEvaluator(),
                minK: 2,
                maxK: 4);

            Assert.IsNotNull(analysis.ClusterAnalysis);
            Assert.IsTrue(analysis.ClusterAnalysis.BestClusterCount >= 2,
                "Should find at least 2 clusters.");
            Assert.IsTrue(analysis.ClusterAnalysis.Labels.Length == 30,
                "One label per iteration.");
        }

        [TestMethod]
        public void AnalyzeWith_DominantClusterBurnProbability()
        {
            var grid = MakeGrid(300, 25);
            var terrain = FlatTerrain(grid);
            var fuelMap = UniformGrass(grid);

            var mcResult = RiskScenario
                .ForWildfire()
                .WithTerrain(terrain)
                .WithFuel(fuelMap)
                .WithIgnition(grid.Nx / 2, grid.Ny / 2)
                .WithWind(5.0, new Vector(1, 0, 0))
                .WithVariation(v => v.WindSpeed(1, 10))
                .OverGrid(grid)
                .OverTime(0, 300, 60)
                .RunMonteCarlo(20, seed: 99);

            var analysis = mcResult.AnalyzeWith(
                new KMeans { Seed = 99 },
                new SilhouetteEvaluator(),
                minK: 2, maxK: 3);

            var prob = analysis.GetClusterBurnProbability(analysis.ClusterAnalysis.DominantCluster);
            Assert.AreEqual(grid.Nx * grid.Ny, prob.Length);

            for (int i = 0; i < prob.Length; i++)
                Assert.IsTrue(prob[i] >= 0 && prob[i] <= 1,
                    $"Cluster burn prob at {i} = {prob[i]} out of [0,1].");
        }

        [TestMethod]
        public void AnalyzeWith_Build_ProducesProbabilityMap()
        {
            var grid = MakeGrid(300, 25);
            var terrain = FlatTerrain(grid);
            var fuelMap = UniformGrass(grid);

            var mcResult = RiskScenario
                .ForWildfire()
                .WithTerrain(terrain)
                .WithFuel(fuelMap)
                .WithIgnition(grid.Nx / 2, grid.Ny / 2)
                .WithWind(5.0, new Vector(1, 0, 0))
                .WithVariation(v => v.WindSpeed(2, 8).Moisture(0.04, 0.10))
                .OverGrid(grid)
                .OverTime(0, 300, 60)
                .RunMonteCarlo(20, seed: 77);

            var result = mcResult
                .AnalyzeWith(
                    new KMeans { Seed = 77 },
                    new SilhouetteEvaluator(),
                    minK: 2, maxK: 3)
                .Build();

            Assert.IsNotNull(result);
            Assert.IsTrue(result.Snapshots.Count > 0, "Should have snapshots.");

            // burnState layer should now contain probability values (0-1)
            var lastSnap = result.Snapshots[result.Snapshots.Count - 1];
            var burnProb = lastSnap.Snapshot.GetLayer("burnState");

            bool hasIntermediate = false;
            for (int i = 0; i < burnProb.Length; i++)
            {
                Assert.IsTrue(burnProb[i] >= 0 && burnProb[i] <= 1,
                    $"Burn probability at {i} = {burnProb[i]} out of [0,1].");
                if (burnProb[i] > 0 && burnProb[i] < 1)
                    hasIntermediate = true;
            }

            Assert.IsTrue(hasIntermediate,
                "Probability map should have intermediate values (not just 0/1) " +
                "when averaging across multiple iterations.");
        }

        [TestMethod]
        public void ToMonteCarloScenarioResult_ValidMatrix()
        {
            var grid = MakeGrid(200, 25);
            var terrain = FlatTerrain(grid);
            var fuelMap = UniformGrass(grid);

            var mcResult = RiskScenario
                .ForWildfire()
                .WithTerrain(terrain)
                .WithFuel(fuelMap)
                .WithIgnition(grid.Nx / 2, grid.Ny / 2)
                .WithWind(3.0, new Vector(1, 0, 0))
                .WithVariation(v => v.WindSpeed(1, 5))
                .OverGrid(grid)
                .OverTime(0, 120, 60)
                .RunMonteCarlo(5, seed: 123);

            var scenarioResult = mcResult.ToMonteCarloScenarioResult();

            Assert.AreEqual(5, scenarioResult.Iterations);
            Assert.AreEqual(grid.CellCount * mcResult.AllSnapshots[0].Count,
                scenarioResult.FeatureCount);
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Firebreak / water exclusion
        // ═══════════════════════════════════════════════════════════════

        #region Firebreak exclusion

        [TestMethod]
        public void FirePerimeter_ExcludesFirebreakCells()
        {
            // Grid with a strip of NoFuel (water) cells down the middle column.
            // Fire ignites on the left side. The perimeter should not include
            // the firebreak cells even though their burnState numeric value (3)
            // exceeds the 0.5 threshold.
            var grid = MakeGrid(250, 25);
            var terrain = FlatTerrain(grid);
            var fuelMap = UniformGrass(grid);

            // Set the middle column to NoFuel (→ Firebreak)
            int midX = grid.Nx / 2;
            for (int iy = 0; iy < grid.Ny; iy++)
                fuelMap.SetFuel(midX, iy, FuelModelType.NoFuel);

            // Ignite on the left side
            var result = RiskScenario
                .ForWildfire()
                .WithTerrain(terrain)
                .WithFuel(fuelMap)
                .WithIgnition(1, grid.Ny / 2)
                .WithWind(3.0, new Vector(1, 0, 0))
                .OverGrid(grid)
                .OverTime(0, 600, 60)
                .RunSingle();

            var lastIdx = result.Snapshots.Count - 1;
            var perimeter = result.GenerateFirePerimeter(lastIdx);

            Assert.IsNotNull(perimeter);

            // Verify that firebreak cells are NOT counted in the perimeter.
            // The exceedance count should only include Burning/Burned cells,
            // not Firebreak cells.
            var lastSnap = result.Snapshots[lastIdx];
            int firebreakCount = 0;
            int burnedCount = 0;
            var bs = lastSnap.Snapshot.GetLayer("burnState");
            for (int i = 0; i < bs.Length; i++)
            {
                int s = (int)bs[i];
                if (s == (int)CellBurnState.Firebreak) firebreakCount++;
                if (s == (int)CellBurnState.Burning || s == (int)CellBurnState.Burned) burnedCount++;
            }

            Assert.IsTrue(firebreakCount > 0, "Should have firebreak cells.");
            Assert.IsTrue(burnedCount > 0, "Should have burned cells.");

            // ExceedanceCellCount should match only burned/burning cells
            Assert.IsTrue(perimeter.ExceedanceCellCount <= burnedCount,
                $"Perimeter exceedance count ({perimeter.ExceedanceCellCount}) should not exceed " +
                $"burned cell count ({burnedCount}). Firebreak cells must be excluded.");
        }

        #endregion
    }
}
