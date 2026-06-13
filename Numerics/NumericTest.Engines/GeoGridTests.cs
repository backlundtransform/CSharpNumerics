using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Linq;

namespace NumericsTests
{
    [TestClass]
    public class GeoGridTests
    {
        // ════════════════════════════════════════════
        //  GeoGrid — construction & dimensions
        // ════════════════════════════════════════════

        [TestMethod]
        public void GeoGrid_Constructor_ComputesCellCounts()
        {
            var grid = new GeoGrid(-500, 500, -500, 500, 0, 100, 10);

            Assert.AreEqual(101, grid.Nx);  // (-500..500)/10 + 1
            Assert.AreEqual(101, grid.Ny);
            Assert.AreEqual(11, grid.Nz);   // (0..100)/10 + 1
            Assert.AreEqual(101 * 101 * 11, grid.CellCount);
        }

        [TestMethod]
        public void GeoGrid_SingleLayerZ_Works()
        {
            // zMin == zMax → single z-layer
            var grid = new GeoGrid(-100, 100, -100, 100, 0, 0, 50);

            Assert.AreEqual(5, grid.Nx);   // (-100..100)/50 + 1
            Assert.AreEqual(5, grid.Ny);
            Assert.AreEqual(1, grid.Nz);   // zMax == zMin → 1
            Assert.AreEqual(25, grid.CellCount);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void GeoGrid_ZeroStep_Throws()
        {
            new GeoGrid(0, 100, 0, 100, 0, 10, 0);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void GeoGrid_InvertedXRange_Throws()
        {
            new GeoGrid(100, -100, 0, 100, 0, 10, 10);
        }

        // ════════════════════════════════════════════
        //  Index round-trip
        // ════════════════════════════════════════════

        [TestMethod]
        public void GeoGrid_Index_RoundTrips()
        {
            var grid = new GeoGrid(0, 40, 0, 30, 0, 20, 10);
            // Nx=5, Ny=4, Nz=3

            for (int iz = 0; iz < grid.Nz; iz++)
                for (int iy = 0; iy < grid.Ny; iy++)
                    for (int ix = 0; ix < grid.Nx; ix++)
                    {
                        int flat = grid.Index(ix, iy, iz);
                        var (rx, ry, rz) = grid.Index3D(flat);
                        Assert.AreEqual(ix, rx, $"ix mismatch at ({ix},{iy},{iz})");
                        Assert.AreEqual(iy, ry, $"iy mismatch at ({ix},{iy},{iz})");
                        Assert.AreEqual(iz, rz, $"iz mismatch at ({ix},{iy},{iz})");
                    }
        }

        [TestMethod]
        public void GeoGrid_FlatIndex_IsContiguous()
        {
            var grid = new GeoGrid(0, 20, 0, 10, 0, 0, 10);
            // Nx=3, Ny=2, Nz=1
            Assert.AreEqual(0, grid.Index(0, 0, 0));
            Assert.AreEqual(1, grid.Index(1, 0, 0));
            Assert.AreEqual(2, grid.Index(2, 0, 0));
            Assert.AreEqual(3, grid.Index(0, 1, 0));
            Assert.AreEqual(5, grid.Index(2, 1, 0));
        }

        // ════════════════════════════════════════════
        //  Cell centre positions
        // ════════════════════════════════════════════

        [TestMethod]
        public void GeoGrid_CellCentre_ReturnsCorrectPositions()
        {
            var grid = new GeoGrid(-100, 100, -100, 100, 0, 50, 50);

            var origin = grid.CellCentre(0, 0, 0);
            Assert.AreEqual(-100, origin.x);
            Assert.AreEqual(-100, origin.y);
            Assert.AreEqual(0, origin.z);

            var mid = grid.CellCentre(2, 2, 1);
            Assert.AreEqual(0, mid.x);
            Assert.AreEqual(0, mid.y);
            Assert.AreEqual(50, mid.z);
        }

        [TestMethod]
        public void GeoGrid_CellCentre_ByFlatIndex_MatchesByTriple()
        {
            var grid = new GeoGrid(0, 30, 0, 20, 0, 10, 10);

            for (int i = 0; i < grid.CellCount; i++)
            {
                var (ix, iy, iz) = grid.Index3D(i);
                var a = grid.CellCentre(ix, iy, iz);
                var b = grid.CellCentre(i);
                Assert.AreEqual(a.x, b.x);
                Assert.AreEqual(a.y, b.y);
                Assert.AreEqual(a.z, b.z);
            }
        }

        // ════════════════════════════════════════════
        //  Nearest index lookup
        // ════════════════════════════════════════════

        [TestMethod]
        public void GeoGrid_NearestIndex_ExactHit()
        {
            var grid = new GeoGrid(0, 100, 0, 100, 0, 50, 10);

            var (ix, iy, iz) = grid.NearestIndex(new Vector(30, 40, 20));
            Assert.AreEqual(3, ix);
            Assert.AreEqual(4, iy);
            Assert.AreEqual(2, iz);
        }

        [TestMethod]
        public void GeoGrid_NearestIndex_ClampsOutOfBounds()
        {
            var grid = new GeoGrid(0, 100, 0, 100, 0, 50, 10);

            var (ix, iy, iz) = grid.NearestIndex(new Vector(-999, 999, -5));
            Assert.AreEqual(0, ix);
            Assert.AreEqual(grid.Ny - 1, iy);
            Assert.AreEqual(0, iz);
        }

        // ════════════════════════════════════════════
        //  Enumeration
        // ════════════════════════════════════════════

        [TestMethod]
        public void GeoGrid_Enumeration_ReturnsAllCells()
        {
            var grid = new GeoGrid(0, 20, 0, 10, 0, 0, 10);
            var cells = grid.ToList();

            Assert.AreEqual(grid.CellCount, cells.Count);

            // First cell should be origin
            Assert.AreEqual(0, cells[0].x);
            Assert.AreEqual(0, cells[0].y);
            Assert.AreEqual(0, cells[0].z);
        }

        // ════════════════════════════════════════════
        //  GeoCell
        // ════════════════════════════════════════════

        [TestMethod]
        public void GeoCell_StoresValues()
        {
            var pos = new Vector(10, 20, 30);
            var cell = new GeoCell(pos, 1.5e-4, 3, 42);

            Assert.AreEqual(10, cell.Position.x);
            Assert.AreEqual(20, cell.Position.y);
            Assert.AreEqual(30, cell.Position.z);
            Assert.AreEqual(1.5e-4, cell.Value, 1e-10);
            Assert.AreEqual(3, cell.TimeIndex);
            Assert.AreEqual(42, cell.GridIndex);
        }

        // ════════════════════════════════════════════
        //  GridSnapshot
        // ════════════════════════════════════════════

        [TestMethod]
        public void GridSnapshot_IndexAccess_Works()
        {
            var grid = new GeoGrid(0, 20, 0, 10, 0, 0, 10);
            // Nx=3, Ny=2, Nz=1 → 6 cells
            var values = new double[] { 1, 2, 3, 4, 5, 6 };
            var snap = new GridSnapshot(grid, values, time: 60.0, timeIndex: 1);

            Assert.AreEqual(1, snap[0]);
            Assert.AreEqual(6, snap[5]);
            Assert.AreEqual(4, snap[0, 1, 0]);  // ix=0, iy=1 → flat=3 → value=4
            Assert.AreEqual(60.0, snap.Time);
            Assert.AreEqual(1, snap.TimeIndex);
        }

        [TestMethod]
        public void GridSnapshot_ValueAt_FindsNearest()
        {
            var grid = new GeoGrid(0, 20, 0, 10, 0, 0, 10);
            var values = new double[] { 10, 20, 30, 40, 50, 60 };
            var snap = new GridSnapshot(grid, values, 0, 0);

            // Position (12, 8, 0) → nearest ix=1, iy=1 → flat=4 → value=50
            double val = snap.ValueAt(new Vector(12, 8, 0));
            Assert.AreEqual(50, val);
        }

        [TestMethod]
        public void GridSnapshot_MinMax()
        {
            var grid = new GeoGrid(0, 20, 0, 10, 0, 0, 10);
            var values = new double[] { -3, 0, 7, 2, 9, 1 };
            var snap = new GridSnapshot(grid, values, 0, 0);

            Assert.AreEqual(-3, snap.Min());
            Assert.AreEqual(9, snap.Max());
        }

        [TestMethod]
        public void GridSnapshot_CellsAbove_FiltersCorrectly()
        {
            var grid = new GeoGrid(0, 20, 0, 10, 0, 0, 10);
            var values = new double[] { 0.1, 0.5, 0.01, 0.8, 0.05, 0.9 };
            var snap = new GridSnapshot(grid, values, 0, 0);

            var above = snap.CellsAbove(0.4);
            Assert.AreEqual(3, above.Count);
        }

        [TestMethod]
        public void GridSnapshot_AllCells_ReturnsAllWithPositions()
        {
            var grid = new GeoGrid(0, 10, 0, 10, 0, 0, 10);
            // Nx=2, Ny=2, Nz=1 → 4 cells
            var values = new double[] { 1, 2, 3, 4 };
            var snap = new GridSnapshot(grid, values, 30.0, 2);

            var cells = snap.AllCells().ToList();
            Assert.AreEqual(4, cells.Count);

            // Every cell should have TimeIndex = 2
            foreach (var c in cells)
                Assert.AreEqual(2, c.TimeIndex);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void GridSnapshot_WrongLength_Throws()
        {
            var grid = new GeoGrid(0, 20, 0, 10, 0, 0, 10);
            new GridSnapshot(grid, new double[] { 1, 2, 3 }, 0, 0);
        }

        [TestMethod]
        public void GridSnapshot_GetValues_ReturnsClone()
        {
            var grid = new GeoGrid(0, 10, 0, 10, 0, 0, 10);
            // Nx=2, Ny=2, Nz=1 → 4 cells
            var values = new double[] { 5, 10, 15, 20 };
            var snap = new GridSnapshot(grid, values, 0, 0);

            var copy = snap.GetValues();
            copy[0] = 999;

            Assert.AreEqual(5, snap[0]); // original unchanged
        }

        // ════════════════════════════════════════════
        //  TimeFrame
        // ════════════════════════════════════════════

        [TestMethod]
        public void TimeFrame_CountAndValues()
        {
            var tf = new TimeFrame(0, 3600, 60);

            Assert.AreEqual(61, tf.Count);  // 0, 60, 120, ..., 3600
            Assert.AreEqual(0, tf.TimeAt(0));
            Assert.AreEqual(60, tf.TimeAt(1));
            Assert.AreEqual(3600, tf.TimeAt(60));
        }

        [TestMethod]
        public void TimeFrame_ToArray_MatchesTimeAt()
        {
            var tf = new TimeFrame(100, 400, 100);
            var arr = tf.ToArray();

            Assert.AreEqual(tf.Count, arr.Length);
            for (int i = 0; i < tf.Count; i++)
                Assert.AreEqual(tf.TimeAt(i), arr[i]);
        }

        [TestMethod]
        public void TimeFrame_NearestIndex_Rounds()
        {
            var tf = new TimeFrame(0, 300, 60);
            // Times: 0, 60, 120, 180, 240, 300

            Assert.AreEqual(0, tf.NearestIndex(25));   // closer to 0
            Assert.AreEqual(1, tf.NearestIndex(50));   // closer to 60
            Assert.AreEqual(2, tf.NearestIndex(120));  // exact
            Assert.AreEqual(5, tf.NearestIndex(9999)); // clamped
            Assert.AreEqual(0, tf.NearestIndex(-10));  // clamped
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void TimeFrame_ZeroStep_Throws()
        {
            new TimeFrame(0, 100, 0);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void TimeFrame_EndBeforeStart_Throws()
        {
            new TimeFrame(200, 100, 10);
        }

        [TestMethod]
        public void TimeFrame_SinglePoint()
        {
            var tf = new TimeFrame(0, 0, 1);

            Assert.AreEqual(1, tf.Count);
            Assert.AreEqual(0, tf.TimeAt(0));
        }
    }
}
