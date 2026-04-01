using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Terrain;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Materials.Fire;
using CSharpNumerics.Physics.Materials.Fire.Enums;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;

namespace NumericTest
{
    [TestClass]
    public class TerrainTests
    {
        // ═══════════════════════════════════════════════════════════════
        //  Helper: create a simple 2-D grid (Nz=1)
        // ═══════════════════════════════════════════════════════════════

        private static GeoGrid MakeGrid(double size, double step) =>
            new GeoGrid(0, size, 0, size, 0, 0, step);

        // ═══════════════════════════════════════════════════════════════
        //  TerrainGrid — flat surface
        // ═══════════════════════════════════════════════════════════════

        #region Flat surface

        [TestMethod]
        public void TerrainGrid_FlatSurface_SlopeIsZero()
        {
            var grid = MakeGrid(100, 10);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => 50.0);

            for (int iy = 0; iy < grid.Ny; iy++)
                for (int ix = 0; ix < grid.Nx; ix++)
                    Assert.AreEqual(0.0, terrain.Slope(ix, iy), 1e-10,
                        $"Flat surface slope at ({ix},{iy}) should be 0.");
        }

        [TestMethod]
        public void TerrainGrid_FlatSurface_AspectIsZero()
        {
            var grid = MakeGrid(100, 10);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => 50.0);

            // Aspect is defined as 0 for flat terrain (no gradient)
            for (int iy = 0; iy < grid.Ny; iy++)
                for (int ix = 0; ix < grid.Nx; ix++)
                    Assert.AreEqual(0.0, terrain.Aspect(ix, iy), 1e-10);
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  TerrainGrid — tilted plane
        // ═══════════════════════════════════════════════════════════════

        #region Tilted plane

        [TestMethod]
        public void TerrainGrid_TiltedPlane_EastSlope_CorrectAngle()
        {
            // z = x → slope 45° everywhere, aspect = East (π/2) since
            // the surface rises towards east, downslope is west (3π/2).
            // Actually: dz/dx = 1, dz/dy = 0 → downslope is -x direction
            // aspect = atan2(-(-1), 0) ... let's compute:
            // downslope: negate gradient → (-1, 0)
            // aspect = atan2(-dzdx, -dzdy) = atan2(-1, 0) = -π/2 → +2π = 3π/2 (West)
            var grid = MakeGrid(100, 10);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => x);

            // Check an interior cell for slope = atan(1) = π/4 = 45°
            int midX = grid.Nx / 2;
            int midY = grid.Ny / 2;
            Assert.AreEqual(Math.PI / 4.0, terrain.Slope(midX, midY), 0.01,
                "z=x should give 45° slope.");

            // Aspect = 3π/2 (West, downslope direction)
            Assert.AreEqual(3 * Math.PI / 2, terrain.Aspect(midX, midY), 0.01,
                "z=x surface slopes down to west.");
        }

        [TestMethod]
        public void TerrainGrid_TiltedPlane_NorthSlope_CorrectAngle()
        {
            // z = y → dz/dx=0, dz/dy=1
            // downslope = (0, -1), aspect = atan2(0, -1) = π (South)
            var grid = MakeGrid(100, 10);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => y);

            int midX = grid.Nx / 2;
            int midY = grid.Ny / 2;

            Assert.AreEqual(Math.PI / 4.0, terrain.Slope(midX, midY), 0.01,
                "z=y should give 45° slope.");
            Assert.AreEqual(Math.PI, terrain.Aspect(midX, midY), 0.01,
                "z=y surface slopes down to south.");
        }

        [TestMethod]
        public void TerrainGrid_GentleSlope_CorrectAngle()
        {
            // z = 0.1 * x → dz/dx = 0.1, slope = atan(0.1) ≈ 5.71°
            var grid = MakeGrid(100, 10);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => 0.1 * x);

            int midX = grid.Nx / 2;
            int midY = grid.Ny / 2;

            double expectedSlope = Math.Atan(0.1);
            Assert.AreEqual(expectedSlope, terrain.Slope(midX, midY), 0.001);
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  TerrainGrid — SlopeInDirection
        // ═══════════════════════════════════════════════════════════════

        #region SlopeInDirection

        [TestMethod]
        public void SlopeInDirection_AlongAspect_MatchesFullSlope()
        {
            // z = x → slope in +x direction = atan(1) = 45°
            // slope in +x means direction = (1, 0, 0)
            var grid = MakeGrid(100, 10);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => x);

            int midX = grid.Nx / 2;
            int midY = grid.Ny / 2;

            double slopeUphill = terrain.SlopeInDirection(midX, midY, new Vector(1, 0, 0));
            Assert.AreEqual(Math.PI / 4.0, slopeUphill, 0.01,
                "Slope in uphill direction should match full slope.");
        }

        [TestMethod]
        public void SlopeInDirection_Perpendicular_IsZero()
        {
            // z = x → gradient along y is 0
            var grid = MakeGrid(100, 10);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => x);

            int midX = grid.Nx / 2;
            int midY = grid.Ny / 2;

            double slopePerp = terrain.SlopeInDirection(midX, midY, new Vector(0, 1, 0));
            Assert.AreEqual(0.0, slopePerp, 0.01,
                "Slope perpendicular to gradient should be 0.");
        }

        [TestMethod]
        public void SlopeInDirection_Downhill_IsNegative()
        {
            // z = x → going in -x direction is downhill
            var grid = MakeGrid(100, 10);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => x);

            int midX = grid.Nx / 2;
            int midY = grid.Ny / 2;

            double slopeDownhill = terrain.SlopeInDirection(midX, midY, new Vector(-1, 0, 0));
            Assert.IsTrue(slopeDownhill < 0,
                $"Downhill slope should be negative: {slopeDownhill}");
        }

        [TestMethod]
        public void SlopeInDirection_Diagonal_IntermediateValue()
        {
            // z = x + y → gradient = (1, 1), slope along (1,0) = atan(1) = 45°
            // slope along (1,1) = atan(√2) ≈ 54.7°
            var grid = MakeGrid(100, 10);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => x + y);

            int midX = grid.Nx / 2;
            int midY = grid.Ny / 2;

            double slopeDiag = terrain.SlopeInDirection(midX, midY, new Vector(1, 1, 0));
            // directional derivative = (1)(1/√2) + (1)(1/√2) = √2
            double expected = Math.Atan(Math.Sqrt(2));
            Assert.AreEqual(expected, slopeDiag, 0.01);
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  TerrainGrid — FromArray
        // ═══════════════════════════════════════════════════════════════

        #region FromArray

        [TestMethod]
        public void TerrainGrid_FromArray_RoundTrip()
        {
            var grid = MakeGrid(20, 10); // 3×3 grid (0,10,20 → 3 cells)
            var elev = new double[grid.Ny, grid.Nx];
            for (int iy = 0; iy < grid.Ny; iy++)
                for (int ix = 0; ix < grid.Nx; ix++)
                    elev[iy, ix] = iy * 100 + ix;

            var terrain = TerrainGrid.FromArray(grid, elev);

            for (int iy = 0; iy < grid.Ny; iy++)
                for (int ix = 0; ix < grid.Nx; ix++)
                    Assert.AreEqual(iy * 100 + ix, terrain[ix, iy], 1e-10);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void TerrainGrid_FromArray_WrongDimensions_Throws()
        {
            var grid = MakeGrid(20, 10);
            var badElev = new double[1, 1];
            TerrainGrid.FromArray(grid, badElev);
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  FuelMap — set/get round-trip
        // ═══════════════════════════════════════════════════════════════

        #region FuelMap

        [TestMethod]
        public void FuelMap_DefaultsToNoFuel()
        {
            var grid = MakeGrid(50, 10);
            var map = new FuelMap(grid);

            Assert.AreEqual(FuelModelType.NoFuel, map.GetFuelType(0, 0));
            Assert.AreEqual(FuelModelType.NoFuel, map.GetFuelType(grid.Nx - 1, grid.Ny - 1));
        }

        [TestMethod]
        public void FuelMap_DefaultMoistureIs008()
        {
            var grid = MakeGrid(50, 10);
            var map = new FuelMap(grid);

            Assert.AreEqual(FuelMap.DefaultMoisture, map.GetMoisture(0, 0), 1e-10);
        }

        [TestMethod]
        public void FuelMap_SetGet_SingleCell()
        {
            var grid = MakeGrid(50, 10);
            var map = new FuelMap(grid);

            map.SetFuel(2, 3, FuelModelType.ShortGrass);
            Assert.AreEqual(FuelModelType.ShortGrass, map.GetFuelType(2, 3));

            var fuel = map.GetFuel(2, 3);
            Assert.AreEqual("Short Grass", fuel.Name);
        }

        [TestMethod]
        public void FuelMap_SetUniformFuel()
        {
            var grid = MakeGrid(50, 10);
            var map = new FuelMap(grid);

            map.SetUniformFuel(FuelModelType.Chaparral);

            for (int iy = 0; iy < grid.Ny; iy++)
                for (int ix = 0; ix < grid.Nx; ix++)
                    Assert.AreEqual(FuelModelType.Chaparral, map.GetFuelType(ix, iy));
        }

        [TestMethod]
        public void FuelMap_MoistureSetGet()
        {
            var grid = MakeGrid(50, 10);
            var map = new FuelMap(grid);

            map.SetMoisture(1, 2, 0.15);
            Assert.AreEqual(0.15, map.GetMoisture(1, 2), 1e-10);

            // Other cells unchanged
            Assert.AreEqual(FuelMap.DefaultMoisture, map.GetMoisture(0, 0), 1e-10);
        }

        [TestMethod]
        public void FuelMap_SetUniformMoisture()
        {
            var grid = MakeGrid(50, 10);
            var map = new FuelMap(grid);

            map.SetUniformMoisture(0.03);

            for (int iy = 0; iy < grid.Ny; iy++)
                for (int ix = 0; ix < grid.Nx; ix++)
                    Assert.AreEqual(0.03, map.GetMoisture(ix, iy), 1e-10);
        }

        [TestMethod]
        public void FuelMap_SetFuelByElevation()
        {
            var grid = MakeGrid(100, 10);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => y); // elevation = y
            var map = new FuelMap(grid);

            var ranges = new List<(double, double, FuelModelType)>
            {
                (0, 50, FuelModelType.ShortGrass),
                (50, 200, FuelModelType.Chaparral)
            };

            map.SetFuelByElevation(terrain, ranges);

            // y=0..40 → ShortGrass, y=50..100 → Chaparral
            // Cell centres at y = 0, 10, 20, 30, 40, 50, 60, ...
            Assert.AreEqual(FuelModelType.ShortGrass, map.GetFuelType(0, 0));   // y=0
            Assert.AreEqual(FuelModelType.ShortGrass, map.GetFuelType(0, 4));   // y=40
            Assert.AreEqual(FuelModelType.Chaparral, map.GetFuelType(0, 5));    // y=50
            Assert.AreEqual(FuelModelType.Chaparral, map.GetFuelType(0, grid.Ny - 1)); // y=100
        }

        #endregion
    }
}
