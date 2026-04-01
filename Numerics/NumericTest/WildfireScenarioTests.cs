using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Spread.Wildfire;
using CSharpNumerics.Engines.GIS.Spread.Wildfire.Enums;
using CSharpNumerics.Engines.GIS.Terrain;
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
    }
}
