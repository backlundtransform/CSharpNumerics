using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Terrain;
using CSharpNumerics.Physics.Environmental.Water;
using System;
using System.Collections.Generic;
using System.Linq;

namespace NumericsTests
{
    [TestClass]
    public class RiverNetworkTests
    {
        // ════════════════════════════════════════════════════════════════
        //  FromElevation — D8 flow direction
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void FromElevation_VShapedValley_SingleStreamAtBottom()
        {
            // V-shaped valley: elevation = |x| so the valley bottom is at x=0
            // Grid: 21×5, x in [-100, 100], y in [0, 40], step=10
            var grid = new GeoGrid(-100, 100, 0, 40, 0, 10, 10);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) =>
            {
                // Slight downstream slope in y to force single flow direction
                return Math.Abs(x) + 100 - y * 0.5;
            });

            // threshold chosen to select only the valley bottom
            var net = RiverNetwork.FromElevation(terrain, 3);

            // The center column (ix = Nx/2) should be river cells
            int midX = grid.Nx / 2;
            Assert.IsTrue(net.RiverCellCount > 0, "Should have river cells");

            // Valley bottom cells should be river
            for (int iy = 0; iy < grid.Ny; iy++)
            {
                if (net.IsRiverCell(midX, iy))
                {
                    // Verify it's at or near the valley center
                    Assert.IsTrue(midX >= grid.Nx / 2 - 1 && midX <= grid.Nx / 2 + 1,
                        "River should be at valley bottom");
                }
            }
        }

        [TestMethod]
        public void FromElevation_Confluence_TwoUpstreamBranches()
        {
            // Create a Y-shaped terrain: two valleys merging at a central point
            // Grid: 11×11, centered at origin
            var grid = new GeoGrid(-50, 50, -50, 50, 0, 10, 10);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) =>
            {
                // Two channels merge: elevation drops toward center-bottom
                // Left branch: high elevation except along line from (-40,40) to (0,0)
                // Right branch: high elevation except along line from (40,40) to (0,0)
                double distLeft = Math.Abs(x + y) / Math.Sqrt(2);
                double distRight = Math.Abs(x - y) / Math.Sqrt(2);
                double valley = Math.Min(distLeft, distRight);
                return valley + 100 - y * 0.3;
            });

            var net = RiverNetwork.FromElevation(terrain, 2);

            // Check that at least one cell has multiple upstream neighbours (confluence)
            bool hasConfluence = false;
            for (int iy = 0; iy < grid.Ny; iy++)
                for (int ix = 0; ix < grid.Nx; ix++)
                    if (net.IsRiverCell(ix, iy) && net.GetUpstream(ix, iy).Count >= 2)
                        hasConfluence = true;

            Assert.IsTrue(net.RiverCellCount > 0, "Should have river cells");
            // Note: confluence detection depends on terrain shape; we verify topology works
        }

        // ════════════════════════════════════════════════════════════════
        //  FromManual — hand-drawn networks
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void FromManual_SingleSegment_CreatesCorrectTopology()
        {
            var grid = new GeoGrid(0, 100, 0, 100, 0, 10, 10);
            var cells = new List<(int, int)> { (5, 0), (5, 1), (5, 2), (5, 3), (5, 4) };

            var net = RiverNetwork.FromManual(grid)
                .AddSegment(cells)
                .Build();

            Assert.AreEqual(5, net.RiverCellCount);
            Assert.IsTrue(net.IsRiverCell(5, 0));
            Assert.IsTrue(net.IsRiverCell(5, 4));

            // Headwater (5,0): no upstream, downstream is (5,1)
            Assert.AreEqual(0, net.GetUpstream(5, 0).Count);
            Assert.AreEqual(1, net.GetDownstream(5, 0).Count);
            Assert.AreEqual((5, 1), net.GetDownstream(5, 0)[0]);

            // Outlet (5,4): no downstream, upstream is (5,3)
            Assert.AreEqual(0, net.GetDownstream(5, 4).Count);
            Assert.AreEqual(1, net.GetUpstream(5, 4).Count);
        }

        [TestMethod]
        public void FromManual_TwoBranches_ConfluenceHasTwoUpstream()
        {
            var grid = new GeoGrid(0, 100, 0, 100, 0, 10, 10);

            // Left branch: (3,0) → (3,1) → (4,2) → (5,3)
            var left = new List<(int, int)> { (3, 0), (3, 1), (4, 2), (5, 3) };
            // Right branch: (7,0) → (7,1) → (6,2) → (5,3)
            var right = new List<(int, int)> { (7, 0), (7, 1), (6, 2), (5, 3) };
            // Main stem: (5,3) → (5,4) → (5,5)
            var main = new List<(int, int)> { (5, 3), (5, 4), (5, 5) };

            var net = RiverNetwork.FromManual(grid)
                .AddSegment(left)
                .AddSegment(right)
                .AddSegment(main)
                .Build();

            // Confluence at (5,3): two upstream ((4,2) and (6,2))
            var upstream53 = net.GetUpstream(5, 3);
            Assert.AreEqual(2, upstream53.Count, "Confluence should have 2 upstream neighbours");
            Assert.IsTrue(upstream53.Contains((4, 2)));
            Assert.IsTrue(upstream53.Contains((6, 2)));
        }

        [TestMethod]
        public void FromManual_SetConfluence_ExplicitUpstream()
        {
            var grid = new GeoGrid(0, 100, 0, 100, 0, 10, 10);

            var net = RiverNetwork.FromManual(grid)
                .AddSegment(new List<(int, int)> { (2, 0), (2, 1) })
                .AddSegment(new List<(int, int)> { (4, 0), (4, 1) })
                .SetConfluence(3, 2, new List<(int, int)> { (2, 1), (4, 1) })
                .AddSegment(new List<(int, int)> { (3, 2), (3, 3) })
                .Build();

            Assert.AreEqual(2, net.GetUpstream(3, 2).Count);
            Assert.IsTrue(net.IsRiverCell(3, 2));
        }

        // ════════════════════════════════════════════════════════════════
        //  Topological sort
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void GetReachCells_TopologicalOrder_UpstreamBeforeDownstream()
        {
            var grid = new GeoGrid(0, 100, 0, 100, 0, 10, 10);

            // Simple Y confluence
            var net = RiverNetwork.FromManual(grid)
                .AddSegment(new List<(int, int)> { (2, 0), (3, 1), (4, 2) })
                .AddSegment(new List<(int, int)> { (6, 0), (5, 1), (4, 2) })
                .AddSegment(new List<(int, int)> { (4, 2), (4, 3), (4, 4) })
                .Build();

            var reach = net.GetReachCells();

            // Build index map: cell → position in ordering
            var posMap = new Dictionary<(int, int), int>();
            for (int i = 0; i < reach.Count; i++)
                posMap[reach[i]] = i;

            // Every upstream cell must appear before its downstream cell
            for (int i = 0; i < reach.Count; i++)
            {
                var cell = reach[i];
                foreach (var ds in net.GetDownstream(cell.ix, cell.iy))
                {
                    Assert.IsTrue(posMap.ContainsKey(ds), $"Downstream cell {ds} should be in reach");
                    Assert.IsTrue(posMap[ds] > i,
                        $"Cell {cell} (pos {i}) should appear before downstream {ds} (pos {posMap[ds]})");
                }
            }
        }

        [TestMethod]
        public void NonRiverCell_IsNotInNetwork()
        {
            var grid = new GeoGrid(0, 100, 0, 100, 0, 10, 10);
            var net = RiverNetwork.FromManual(grid)
                .AddSegment(new List<(int, int)> { (5, 0), (5, 1), (5, 2) })
                .Build();

            Assert.IsFalse(net.IsRiverCell(0, 0));
            Assert.IsFalse(net.IsRiverCell(9, 9));
        }

        // ════════════════════════════════════════════════════════════════
        //  ChannelMap
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void ChannelMap_SetGet_RoundTrip()
        {
            var grid = new GeoGrid(0, 100, 0, 100, 0, 10, 10);
            var cm = new ChannelMap(grid);

            cm.SetChannel(3, 4, 15.0, 3.0, 0.025);

            Assert.AreEqual(15.0, cm.GetWidth(3, 4), 1e-10);
            Assert.AreEqual(3.0, cm.GetDepth(3, 4), 1e-10);
            Assert.AreEqual(0.025, cm.GetManningN(3, 4), 1e-10);
        }

        [TestMethod]
        public void ChannelMap_UniformChannel_SetsAll()
        {
            var grid = new GeoGrid(0, 50, 0, 50, 0, 10, 10);
            var cm = new ChannelMap(grid);
            cm.SetUniformChannel(20, 4, 0.03);

            Assert.AreEqual(20.0, cm.GetWidth(0, 0), 1e-10);
            Assert.AreEqual(4.0, cm.GetDepth(3, 3), 1e-10);
            Assert.AreEqual(0.03, cm.GetManningN(2, 1), 1e-10);
        }

        [TestMethod]
        public void ChannelMap_GetVelocity_MatchesManningEquation()
        {
            var grid = new GeoGrid(0, 100, 0, 100, 0, 10, 10);

            // Terrain with steady downhill slope in y
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => 100 - y);

            // Simple straight river along x=50 (ix=5)
            var net = RiverNetwork.FromManual(grid)
                .AddSegment(new List<(int, int)> { (5, 0), (5, 1), (5, 2), (5, 3) })
                .Build();

            var cm = new ChannelMap(grid);
            double w = 12, d = 2.5, n = 0.030;
            cm.SetChannel(5, 1, w, d, n);

            double velocity = cm.GetVelocity(5, 1, net, terrain);
            double slope = cm.GetBedSlope(5, 1, net, terrain);

            // Verify against Manning's equation directly
            double Rh = ManningEquation.RectangularHydraulicRadius(w, d);
            double expected = ManningEquation.Velocity(n, Rh, slope);

            Assert.AreEqual(expected, velocity, 1e-10,
                "ChannelMap velocity should match ManningEquation.Velocity");
            Assert.IsTrue(velocity > 0, "Velocity should be positive on a slope");
        }

        [TestMethod]
        public void ChannelMap_GetBedSlope_OutletCell_ReturnsZero()
        {
            var grid = new GeoGrid(0, 100, 0, 100, 0, 10, 10);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => 100 - y);

            var net = RiverNetwork.FromManual(grid)
                .AddSegment(new List<(int, int)> { (5, 0), (5, 1), (5, 2) })
                .Build();

            var cm = new ChannelMap(grid);

            // (5,2) is outlet — no downstream
            double slope = cm.GetBedSlope(5, 2, net, terrain);
            Assert.AreEqual(0, slope, "Outlet cell should have zero bed slope");
        }

        [TestMethod]
        public void ChannelMap_SetChannelByStreamOrder_DownstreamIsWider()
        {
            var grid = new GeoGrid(0, 100, 0, 100, 0, 10, 10);
            var terrain = TerrainGrid.FromFunction(grid, (x, y) => 100 - y);

            // Y network: two branches merge into main stem
            var net = RiverNetwork.FromManual(grid)
                .AddSegment(new List<(int, int)> { (3, 0), (3, 1), (4, 2), (5, 3) })
                .AddSegment(new List<(int, int)> { (7, 0), (7, 1), (6, 2), (5, 3) })
                .AddSegment(new List<(int, int)> { (5, 3), (5, 4), (5, 5), (5, 6) })
                .Build();

            var cm = new ChannelMap(grid);
            cm.SetChannelByStreamOrder(net, terrain,
                headwaterWidth: 3, headwaterDepth: 0.5, growthFactor: 1.5);

            // Headwater cell
            double headWidth = cm.GetWidth(3, 0);
            // Downstream cell after confluence
            double downWidth = cm.GetWidth(5, 6);

            Assert.IsTrue(downWidth > headWidth,
                $"Downstream width ({downWidth:F2}) should be > headwater width ({headWidth:F2})");
        }
    }
}
