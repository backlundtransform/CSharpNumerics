using CSharpNumerics.Engines.GIS.Grid;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Engines.GIS.Terrain;

/// <summary>
/// A directed flow graph over a <see cref="GeoGrid"/> representing a river/stream
/// network. Each river cell knows its upstream and downstream neighbours, enabling
/// topologically-ordered contaminant transport.
/// <para>
/// Create via <see cref="FromElevation"/> (D8 flow-direction algorithm) or
/// <see cref="FromManual"/> (hand-drawn river paths).
/// </para>
/// </summary>
public class RiverNetwork
{
    private readonly bool[] _isRiver;
    private readonly List<(int ix, int iy)>[] _downstream;
    private readonly List<(int ix, int iy)>[] _upstream;
    private List<(int ix, int iy)> _topologicalOrder;

    /// <summary>The underlying spatial grid.</summary>
    public GeoGrid Grid { get; }

    /// <summary>Number of cells along x.</summary>
    public int Nx => Grid.Nx;

    /// <summary>Number of cells along y.</summary>
    public int Ny => Grid.Ny;

    /// <summary>Total number of river cells.</summary>
    public int RiverCellCount { get; private set; }

    private RiverNetwork(GeoGrid grid)
    {
        Grid = grid ?? throw new ArgumentNullException(nameof(grid));
        int count = grid.Nx * grid.Ny;
        _isRiver = new bool[count];
        _downstream = new List<(int, int)>[count];
        _upstream = new List<(int, int)>[count];
        for (int i = 0; i < count; i++)
        {
            _downstream[i] = new List<(int, int)>();
            _upstream[i] = new List<(int, int)>();
        }
    }

    private int Idx(int ix, int iy) => iy * Nx + ix;

    // ═══════════════════════════════════════════════════════════════
    //  Queries
    // ═══════════════════════════════════════════════════════════════

    /// <summary>Returns true if cell (ix, iy) is part of the river network.</summary>
    public bool IsRiverCell(int ix, int iy) => _isRiver[Idx(ix, iy)];

    /// <summary>
    /// Returns the downstream neighbours of cell (ix, iy).
    /// Usually 1 cell; more than 1 at a bifurcation.
    /// Empty for outlet cells.
    /// </summary>
    public IReadOnlyList<(int ix, int iy)> GetDownstream(int ix, int iy) =>
        _downstream[Idx(ix, iy)];

    /// <summary>
    /// Returns the upstream neighbours of cell (ix, iy).
    /// More than 1 at a confluence. Empty for headwater cells.
    /// </summary>
    public IReadOnlyList<(int ix, int iy)> GetUpstream(int ix, int iy) =>
        _upstream[Idx(ix, iy)];

    /// <summary>
    /// Returns all river cells in topological order (upstream → downstream).
    /// Suitable for forward-marching advection schemes.
    /// </summary>
    public IReadOnlyList<(int ix, int iy)> GetReachCells()
    {
        if (_topologicalOrder == null)
            _topologicalOrder = ComputeTopologicalOrder();
        return _topologicalOrder;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Factory: FromElevation (D8 flow direction)
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Builds a river network from terrain elevation using the D8
    /// flow-direction algorithm. Each cell drains to its steepest
    /// downhill neighbour (8-connected). Flow accumulation is computed
    /// and cells with accumulation ≥ <paramref name="flowThreshold"/>
    /// are marked as river.
    /// </summary>
    /// <param name="terrain">Elevation surface.</param>
    /// <param name="flowThreshold">Minimum flow accumulation (in cell units)
    /// to classify a cell as river. Higher values produce sparser networks.</param>
    public static RiverNetwork FromElevation(TerrainGrid terrain, double flowThreshold)
    {
        if (terrain == null) throw new ArgumentNullException(nameof(terrain));

        var grid = terrain.Grid;
        var net = new RiverNetwork(grid);
        int nx = grid.Nx;
        int ny = grid.Ny;
        int count = nx * ny;

        // Step 1: D8 flow direction — each cell points to its steepest downhill neighbour
        var flowDir = new int[count]; // flat index of target, or -1 if pit/outlet
        for (int i = 0; i < count; i++) flowDir[i] = -1;

        int[] dix = { -1, 0, 1, -1, 1, -1, 0, 1 };
        int[] diy = { -1, -1, -1, 0, 0, 1, 1, 1 };
        double diag = Math.Sqrt(2);
        double[] dist = { diag, 1, diag, 1, 1, diag, 1, diag };

        for (int iy = 0; iy < ny; iy++)
        {
            for (int ix = 0; ix < nx; ix++)
            {
                double elev = terrain[ix, iy];
                double steepest = 0;
                int bestIdx = -1;

                for (int d = 0; d < 8; d++)
                {
                    int nix = ix + dix[d];
                    int niy = iy + diy[d];
                    if (nix < 0 || nix >= nx || niy < 0 || niy >= ny) continue;

                    double drop = (elev - terrain[nix, niy]) / (dist[d] * grid.Step);
                    if (drop > steepest)
                    {
                        steepest = drop;
                        bestIdx = niy * nx + nix;
                    }
                }

                flowDir[iy * nx + ix] = bestIdx;
            }
        }

        // Step 2: Flow accumulation — count upstream contributing cells
        var accumulation = new int[count];
        for (int i = 0; i < count; i++) accumulation[i] = 1; // self

        // Sort cells by elevation descending so we process high cells first
        var indices = Enumerable.Range(0, count)
            .OrderByDescending(i => terrain[i % nx, i / nx])
            .ToArray();

        foreach (int i in indices)
        {
            if (flowDir[i] >= 0)
                accumulation[flowDir[i]] += accumulation[i];
        }

        // Step 3: Mark river cells and build directed graph
        for (int i = 0; i < count; i++)
        {
            if (accumulation[i] >= flowThreshold)
            {
                int ix = i % nx;
                int iy = i / nx;
                net._isRiver[i] = true;
                net.RiverCellCount++;
            }
        }

        // Build edges only between river cells
        for (int i = 0; i < count; i++)
        {
            if (!net._isRiver[i]) continue;
            int target = flowDir[i];
            if (target < 0 || !net._isRiver[target]) continue;

            int srcX = i % nx, srcY = i / nx;
            int tgtX = target % nx, tgtY = target / nx;

            net._downstream[i].Add((tgtX, tgtY));
            net._upstream[target].Add((srcX, srcY));
        }

        return net;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Factory: FromManual (hand-drawn river paths)
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Returns a <see cref="ManualRiverBuilder"/> for constructing a river
    /// network by explicitly adding segments and confluences.
    /// </summary>
    public static ManualRiverBuilder FromManual(GeoGrid grid) =>
        new ManualRiverBuilder(grid);

    // ═══════════════════════════════════════════════════════════════
    //  Topological sort (Kahn's algorithm)
    // ═══════════════════════════════════════════════════════════════

    private List<(int ix, int iy)> ComputeTopologicalOrder()
    {
        var result = new List<(int ix, int iy)>(RiverCellCount);
        var inDegree = new int[Nx * Ny];

        // Count in-degree (number of upstream neighbours) for each river cell
        for (int iy = 0; iy < Ny; iy++)
            for (int ix = 0; ix < Nx; ix++)
                if (_isRiver[Idx(ix, iy)])
                    inDegree[Idx(ix, iy)] = _upstream[Idx(ix, iy)].Count;

        // Start with headwater cells (in-degree = 0)
        var queue = new Queue<(int ix, int iy)>();
        for (int iy = 0; iy < Ny; iy++)
            for (int ix = 0; ix < Nx; ix++)
                if (_isRiver[Idx(ix, iy)] && inDegree[Idx(ix, iy)] == 0)
                    queue.Enqueue((ix, iy));

        while (queue.Count > 0)
        {
            var cell = queue.Dequeue();
            result.Add(cell);

            foreach (var ds in _downstream[Idx(cell.ix, cell.iy)])
            {
                int dsIdx = Idx(ds.ix, ds.iy);
                inDegree[dsIdx]--;
                if (inDegree[dsIdx] == 0)
                    queue.Enqueue(ds);
            }
        }

        return result;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Manual builder
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Fluent builder for constructing a <see cref="RiverNetwork"/> by hand.
    /// Add river segments as ordered lists of cells (upstream → downstream),
    /// then call <see cref="Build"/>.
    /// </summary>
    public class ManualRiverBuilder
    {
        private readonly RiverNetwork _net;

        internal ManualRiverBuilder(GeoGrid grid)
        {
            _net = new RiverNetwork(grid);
        }

        /// <summary>
        /// Adds a river segment as an ordered list of cells from upstream to downstream.
        /// Consecutive cells in the list are linked as upstream → downstream.
        /// </summary>
        public ManualRiverBuilder AddSegment(IReadOnlyList<(int ix, int iy)> cells)
        {
            if (cells == null) throw new ArgumentNullException(nameof(cells));

            for (int i = 0; i < cells.Count; i++)
            {
                int idx = _net.Idx(cells[i].ix, cells[i].iy);
                if (!_net._isRiver[idx])
                {
                    _net._isRiver[idx] = true;
                    _net.RiverCellCount++;
                }

                if (i < cells.Count - 1)
                {
                    var from = cells[i];
                    var to = cells[i + 1];
                    int fromIdx = _net.Idx(from.ix, from.iy);
                    int toIdx = _net.Idx(to.ix, to.iy);

                    if (!_net._downstream[fromIdx].Contains(to))
                        _net._downstream[fromIdx].Add(to);
                    if (!_net._upstream[toIdx].Contains(from))
                        _net._upstream[toIdx].Add(from);
                }
            }

            return this;
        }

        /// <summary>
        /// Explicitly sets a confluence: multiple upstream cells feed into one downstream cell.
        /// All cells involved are marked as river cells.
        /// </summary>
        public ManualRiverBuilder SetConfluence(int ix, int iy, IReadOnlyList<(int ix, int iy)> upstreamCells)
        {
            if (upstreamCells == null) throw new ArgumentNullException(nameof(upstreamCells));

            int dsIdx = _net.Idx(ix, iy);
            if (!_net._isRiver[dsIdx])
            {
                _net._isRiver[dsIdx] = true;
                _net.RiverCellCount++;
            }

            foreach (var us in upstreamCells)
            {
                int usIdx = _net.Idx(us.ix, us.iy);
                if (!_net._isRiver[usIdx])
                {
                    _net._isRiver[usIdx] = true;
                    _net.RiverCellCount++;
                }

                if (!_net._downstream[usIdx].Contains((ix, iy)))
                    _net._downstream[usIdx].Add((ix, iy));
                if (!_net._upstream[dsIdx].Contains(us))
                    _net._upstream[dsIdx].Add(us);
            }

            return this;
        }

        /// <summary>
        /// Finalises and returns the river network.
        /// </summary>
        public RiverNetwork Build() => _net;
    }
}
