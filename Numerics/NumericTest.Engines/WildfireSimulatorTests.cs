using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Spread;
using CSharpNumerics.Engines.GIS.Spread.Wildfire;
using CSharpNumerics.Engines.GIS.Spread.Wildfire.Enums;
using CSharpNumerics.Engines.GIS.Terrain;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Materials.Fire;
using CSharpNumerics.Physics.Materials.Fire.Enums;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;

namespace NumericTest
{
    [TestClass]
    public class WildfireSimulatorTests
    {
        // ═══════════════════════════════════════════════════════════════
        //  Helpers
        // ═══════════════════════════════════════════════════════════════

        /// <summary>Create a flat 2-D grid.</summary>
        private static GeoGrid MakeGrid(double size, double step) =>
            new GeoGrid(0, size, 0, size, 0, 0, step);

        /// <summary>Flat terrain (elevation = 0 everywhere).</summary>
        private static TerrainGrid FlatTerrain(GeoGrid grid) =>
            TerrainGrid.FromFunction(grid, (x, y) => 0);

        /// <summary>Uniform short-grass fuel map with given moisture.</summary>
        private static FuelMap UniformGrass(GeoGrid grid, double moisture = 0.05)
        {
            var map = new FuelMap(grid);
            map.SetUniformFuel(FuelModelType.ShortGrass);
            map.SetUniformMoisture(moisture);
            return map;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Circular spread (no wind, flat)
        // ═══════════════════════════════════════════════════════════════

        #region Circular spread

        [TestMethod]
        public void NoWind_Flat_NearCircularSpread()
        {
            // 500m grid, 25m cells, ignition at centre
            var grid = MakeGrid(500, 25);
            var terrain = FlatTerrain(grid);
            var fuelMap = UniformGrass(grid);

            int cx = grid.Nx / 2;
            int cy = grid.Ny / 2;

            var parameters = new WildfireParameters(
                ignitionPoints: new[] { (cx, cy) },
                midflameWindSpeed: 0,
                burnDuration: 600);

            var simulator = new WildfireSimulator(parameters);
            var timeFrame = new TimeFrame(0, 600, 60); // 10 min, 1-min steps
            var snapshots = simulator.Run(grid, terrain, fuelMap, timeFrame);

            Assert.IsTrue(snapshots.Count > 0, "Should produce snapshots.");

            // After several steps, fire should have spread
            var last = snapshots[snapshots.Count - 1];
            int burningOrBurned = last.BurningCellCount + last.BurnedCellCount;
            Assert.IsTrue(burningOrBurned > 1,
                $"Fire should spread from ignition: {burningOrBurned} cells affected.");

            // Check roughly circular: measure max x/y extents of burned cells
            // and verify they're similar in all directions
            double minX = double.MaxValue, maxX = double.MinValue;
            double minY = double.MaxValue, maxY = double.MinValue;
            for (int iy = 0; iy < grid.Ny; iy++)
            {
                for (int ix = 0; ix < grid.Nx; ix++)
                {
                    var s = last.GetBurnState(ix, iy);
                    if (s == CellBurnState.Burning || s == CellBurnState.Burned)
                    {
                        if (ix < minX) minX = ix;
                        if (ix > maxX) maxX = ix;
                        if (iy < minY) minY = iy;
                        if (iy > maxY) maxY = iy;
                    }
                }
            }

            double xExtent = maxX - minX;
            double yExtent = maxY - minY;

            // For circular spread, x and y extents should be similar (within 50%)
            if (xExtent > 0 && yExtent > 0)
            {
                double ratio = Math.Min(xExtent, yExtent) / Math.Max(xExtent, yExtent);
                Assert.IsTrue(ratio > 0.5,
                    $"Near-circular spread expected: xExtent={xExtent}, yExtent={yExtent}, ratio={ratio}");
            }
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Elliptical spread (with wind)
        // ═══════════════════════════════════════════════════════════════

        #region Wind-driven spread

        [TestMethod]
        public void Wind_Flat_EllipticalSpread()
        {
            // Strong wind in +x direction should elongate spread eastward
            var grid = MakeGrid(1000, 25);
            var terrain = FlatTerrain(grid);
            var fuelMap = UniformGrass(grid);

            int cx = grid.Nx / 2;
            int cy = grid.Ny / 2;

            var parameters = new WildfireParameters(
                ignitionPoints: new[] { (cx, cy) },
                midflameWindSpeed: 5.0,                     // 5 m/s
                windDirection: new Vector(1, 0, 0),          // blowing east
                burnDuration: 1200);

            var simulator = new WildfireSimulator(parameters);
            var timeFrame = new TimeFrame(0, 600, 60);
            var snapshots = simulator.Run(grid, terrain, fuelMap, timeFrame);

            var last = snapshots[snapshots.Count - 1];

            // Find max east and west extents from centre
            int maxEast = 0, maxWest = 0;
            for (int ix = 0; ix < grid.Nx; ix++)
            {
                var s = last.GetBurnState(ix, cy);
                if (s == CellBurnState.Burning || s == CellBurnState.Burned)
                {
                    int dist = ix - cx;
                    if (dist > maxEast) maxEast = dist;
                    if (-dist > maxWest) maxWest = -dist;
                }
            }

            // Downwind (east) should extend further than upwind (west)
            Assert.IsTrue(maxEast > maxWest,
                $"Downwind extent ({maxEast}) should exceed upwind ({maxWest}).");
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Slope effects
        // ═══════════════════════════════════════════════════════════════

        #region Slope effects

        [TestMethod]
        public void UphillSlope_AcceleratesSpread()
        {
            // z = 0.3*y → uphill to the north
            var grid = MakeGrid(500, 25);
            var slopeTerrain = TerrainGrid.FromFunction(grid, (x, y) => 0.3 * y);
            var flatTerrain = FlatTerrain(grid);
            var fuelMap = UniformGrass(grid);

            int cx = grid.Nx / 2;
            int cy = grid.Ny / 2;

            var parameters = new WildfireParameters(
                ignitionPoints: new[] { (cx, cy) },
                midflameWindSpeed: 0,
                burnDuration: 1200);

            // Run on slope
            var simSlope = new WildfireSimulator(parameters);
            var tf = new TimeFrame(0, 600, 60);
            var snapsSlope = simSlope.Run(grid, slopeTerrain, fuelMap, tf);

            // Run on flat
            var simFlat = new WildfireSimulator(parameters);
            var snapsFlat = simFlat.Run(grid, flatTerrain, fuelMap, tf);

            var lastSlope = snapsSlope[snapsSlope.Count - 1];
            var lastFlat = snapsFlat[snapsFlat.Count - 1];

            // On slope, the northward (uphill) extent from centre should be
            // greater than on flat terrain
            int maxNorthSlope = 0, maxNorthFlat = 0;
            for (int iy = cy; iy < grid.Ny; iy++)
            {
                var sSlope = lastSlope.GetBurnState(cx, iy);
                if (sSlope == CellBurnState.Burning || sSlope == CellBurnState.Burned)
                    maxNorthSlope = iy - cy;

                var sFlat = lastFlat.GetBurnState(cx, iy);
                if (sFlat == CellBurnState.Burning || sFlat == CellBurnState.Burned)
                    maxNorthFlat = iy - cy;
            }

            Assert.IsTrue(maxNorthSlope >= maxNorthFlat,
                $"Uphill spread ({maxNorthSlope}) should be >= flat spread ({maxNorthFlat}).");
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Firebreaks and NoFuel
        // ═══════════════════════════════════════════════════════════════

        #region Firebreaks

        [TestMethod]
        public void NoFuel_CellsBlockFire()
        {
            // Small grid, barrier of NoFuel across the middle
            var grid = MakeGrid(200, 25);
            var terrain = FlatTerrain(grid);
            var fuelMap = UniformGrass(grid);

            // Create a NoFuel barrier at ix = grid.Nx/2 for all iy
            int barrierX = grid.Nx / 2;
            for (int iy = 0; iy < grid.Ny; iy++)
                fuelMap.SetFuel(barrierX, iy, FuelModelType.NoFuel);

            // Ignite west of barrier
            var parameters = new WildfireParameters(
                ignitionPoints: new[] { (barrierX - 2, grid.Ny / 2) },
                midflameWindSpeed: 3.0,
                windDirection: new Vector(1, 0, 0), // wind blowing east
                burnDuration: 1200);

            var simulator = new WildfireSimulator(parameters);
            var tf = new TimeFrame(0, 600, 60);
            var snapshots = simulator.Run(grid, terrain, fuelMap, tf);

            var last = snapshots[snapshots.Count - 1];

            // No cell east of the barrier should be burning or burned
            for (int iy = 0; iy < grid.Ny; iy++)
            {
                for (int ix = barrierX + 1; ix < grid.Nx; ix++)
                {
                    var s = last.GetBurnState(ix, iy);
                    Assert.AreEqual(CellBurnState.Unburned, s,
                        $"Cell ({ix},{iy}) east of NoFuel barrier should be unburned.");
                }
            }
        }

        [TestMethod]
        public void Firebreak_CellState_BlocksSpread()
        {
            var grid = MakeGrid(200, 25);
            var terrain = FlatTerrain(grid);
            var fuelMap = UniformGrass(grid);

            // NoFuel cells become Firebreak state
            int barrierX = grid.Nx / 2;
            for (int iy = 0; iy < grid.Ny; iy++)
                fuelMap.SetFuel(barrierX, iy, FuelModelType.NoFuel);

            var parameters = new WildfireParameters(
                ignitionPoints: new[] { (1, grid.Ny / 2) },
                midflameWindSpeed: 0,
                burnDuration: 1200);

            var simulator = new WildfireSimulator(parameters);
            var tf = new TimeFrame(0, 300, 60);
            var snapshots = simulator.Run(grid, terrain, fuelMap, tf);

            // Barrier cells should be in Firebreak state
            var last = snapshots[snapshots.Count - 1];
            Assert.AreEqual(CellBurnState.Firebreak, last.GetBurnState(barrierX, grid.Ny / 2));
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Burned area monotonicity
        // ═══════════════════════════════════════════════════════════════

        #region Monotonicity

        [TestMethod]
        public void BurnedAreaHectares_IncreasesMonotonically()
        {
            var grid = MakeGrid(500, 25);
            var terrain = FlatTerrain(grid);
            var fuelMap = UniformGrass(grid);

            var parameters = new WildfireParameters(
                ignitionPoints: new[] { (grid.Nx / 2, grid.Ny / 2) },
                midflameWindSpeed: 2.0,
                windDirection: new Vector(1, 1, 0),
                burnDuration: 1200);

            var simulator = new WildfireSimulator(parameters);
            var tf = new TimeFrame(0, 600, 60);
            var snapshots = simulator.Run(grid, terrain, fuelMap, tf);

            double prevArea = 0;
            foreach (var snap in snapshots)
            {
                double area = snap.BurnedAreaHectares;
                Assert.IsTrue(area >= prevArea,
                    $"Burned area should be monotonically non-decreasing: " +
                    $"t={snap.Time}s area={area} < prev={prevArea}");
                prevArea = area;
            }
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  High moisture — fire does not spread
        // ═══════════════════════════════════════════════════════════════

        #region Saturated fuel

        [TestMethod]
        public void HighMoisture_FireDoesNotSpread()
        {
            var grid = MakeGrid(200, 25);
            var terrain = FlatTerrain(grid);
            var fuelMap = new FuelMap(grid);
            fuelMap.SetUniformFuel(FuelModelType.ShortGrass);
            fuelMap.SetUniformMoisture(0.99); // way above extinction (0.12)

            int cx = grid.Nx / 2;
            int cy = grid.Ny / 2;

            var parameters = new WildfireParameters(
                ignitionPoints: new[] { (cx, cy) },
                midflameWindSpeed: 5.0,
                windDirection: new Vector(1, 0, 0),
                burnDuration: 1200);

            var simulator = new WildfireSimulator(parameters);
            var tf = new TimeFrame(0, 600, 60);
            var snapshots = simulator.Run(grid, terrain, fuelMap, tf);

            // Only the ignition cell should have burned/burning,
            // no neighbours should ignite
            var last = snapshots[snapshots.Count - 1];
            int affected = 0;
            for (int iy = 0; iy < grid.Ny; iy++)
                for (int ix = 0; ix < grid.Nx; ix++)
                {
                    var s = last.GetBurnState(ix, iy);
                    if (s == CellBurnState.Burning || s == CellBurnState.Burned)
                        affected++;
                }

            Assert.AreEqual(1, affected,
                $"Only ignition cell should be affected with saturated fuel: {affected}");
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  SpreadSnapshot layer integrity
        // ═══════════════════════════════════════════════════════════════

        #region Snapshot layers

        [TestMethod]
        public void Snapshot_HasAllLayers()
        {
            var grid = MakeGrid(200, 25);
            var terrain = FlatTerrain(grid);
            var fuelMap = UniformGrass(grid);

            var parameters = new WildfireParameters(
                ignitionPoints: new[] { (grid.Nx / 2, grid.Ny / 2) },
                burnDuration: 600);

            var simulator = new WildfireSimulator(parameters);
            var tf = new TimeFrame(0, 120, 60);
            var snapshots = simulator.Run(grid, terrain, fuelMap, tf);

            Assert.IsTrue(snapshots.Count > 0);
            var snap = snapshots[0];

            Assert.IsTrue(snap.Snapshot.HasLayer("burnState"));
            Assert.IsTrue(snap.Snapshot.HasLayer("flameLength"));
            Assert.IsTrue(snap.Snapshot.HasLayer("rateOfSpread"));
            Assert.IsTrue(snap.Snapshot.HasLayer("burnTime"));
        }

        [TestMethod]
        public void Snapshot_BurnTime_PositiveForBurningCells()
        {
            var grid = MakeGrid(500, 25);
            var terrain = FlatTerrain(grid);
            var fuelMap = UniformGrass(grid);

            var parameters = new WildfireParameters(
                ignitionPoints: new[] { (grid.Nx / 2, grid.Ny / 2) },
                midflameWindSpeed: 2.0,
                windDirection: new Vector(1, 0, 0),
                burnDuration: 1200);

            var simulator = new WildfireSimulator(parameters);
            var tf = new TimeFrame(0, 300, 60);
            var snapshots = simulator.Run(grid, terrain, fuelMap, tf);

            // In later snapshots, burning cells should have positive burn time
            var last = snapshots[snapshots.Count - 1];
            var burnTime = last.Snapshot.GetLayer("burnTime");
            var burnState = last.Snapshot.GetLayer("burnState");

            for (int i = 0; i < burnState.Length; i++)
            {
                if ((int)burnState[i] == (int)CellBurnState.Burning)
                    Assert.IsTrue(burnTime[i] >= 0,
                        $"Burning cell {i} should have non-negative burn time: {burnTime[i]}");
            }
        }

        #endregion
    }
}
