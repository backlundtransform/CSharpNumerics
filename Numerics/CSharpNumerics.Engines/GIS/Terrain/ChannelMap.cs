using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Physics.Environmental.Water;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Terrain;

/// <summary>
/// Per-cell channel hydraulic properties on a <see cref="GeoGrid"/>.
/// Assigns width, depth, and Manning's n to each river cell, and derives
/// velocity and bed slope from the terrain and river network.
/// </summary>
public class ChannelMap
{
    private readonly double[] _width;
    private readonly double[] _depth;
    private readonly double[] _manningN;

    /// <summary>The underlying spatial grid.</summary>
    public GeoGrid Grid { get; }

    /// <summary>Number of cells along x.</summary>
    public int Nx => Grid.Nx;

    /// <summary>Number of cells along y.</summary>
    public int Ny => Grid.Ny;

    /// <summary>Default Manning's n for new cells.</summary>
    public const double DefaultManningN = 0.035;

    /// <summary>Default channel width (m).</summary>
    public const double DefaultWidth = 10.0;

    /// <summary>Default channel depth (m).</summary>
    public const double DefaultDepth = 2.0;

    /// <summary>
    /// Creates a channel map for the given grid. All cells default to
    /// <see cref="DefaultWidth"/>, <see cref="DefaultDepth"/>, and <see cref="DefaultManningN"/>.
    /// </summary>
    public ChannelMap(GeoGrid grid)
    {
        Grid = grid ?? throw new ArgumentNullException(nameof(grid));
        int count = grid.Nx * grid.Ny;
        _width = new double[count];
        _depth = new double[count];
        _manningN = new double[count];
        for (int i = 0; i < count; i++)
        {
            _width[i] = DefaultWidth;
            _depth[i] = DefaultDepth;
            _manningN[i] = DefaultManningN;
        }
    }

    private int Idx(int ix, int iy) => iy * Nx + ix;

    // ═══════════════════════════════════════════════════════════════
    //  Channel assignment
    // ═══════════════════════════════════════════════════════════════

    /// <summary>Sets channel properties for cell (ix, iy).</summary>
    public void SetChannel(int ix, int iy, double width, double depth, double manningN)
    {
        int idx = Idx(ix, iy);
        _width[idx] = width;
        _depth[idx] = depth;
        _manningN[idx] = manningN;
    }

    /// <summary>Sets all cells to the same channel properties.</summary>
    public void SetUniformChannel(double width, double depth, double manningN)
    {
        for (int i = 0; i < _width.Length; i++)
        {
            _width[i] = width;
            _depth[i] = depth;
            _manningN[i] = manningN;
        }
    }

    /// <summary>
    /// Assigns channel geometry by Strahler stream order approximation.
    /// Cells further downstream (more upstream contributing cells) get
    /// wider and deeper channels.
    /// </summary>
    /// <param name="net">River network for topology.</param>
    /// <param name="terrain">Terrain grid (unused in this overload but reserved for future slope-based sizing).</param>
    /// <param name="headwaterWidth">Width (m) for headwater cells.</param>
    /// <param name="headwaterDepth">Depth (m) for headwater cells.</param>
    /// <param name="growthFactor">Multiplicative factor per confluence level (typical: 1.3–1.8).</param>
    /// <param name="manningN">Manning's roughness coefficient for all cells.</param>
    public void SetChannelByStreamOrder(
        RiverNetwork net, TerrainGrid terrain,
        double headwaterWidth = 3.0, double headwaterDepth = 0.5,
        double growthFactor = 1.5, double manningN = DefaultManningN)
    {
        if (net == null) throw new ArgumentNullException(nameof(net));

        // Compute upstream count for each river cell (number of contributing cells)
        var reachCells = net.GetReachCells();
        var upstreamCount = new int[Nx * Ny];
        foreach (var cell in reachCells)
            upstreamCount[Idx(cell.ix, cell.iy)] = 1; // self

        // Accumulate upstream counts in topological order
        foreach (var cell in reachCells)
        {
            int idx = Idx(cell.ix, cell.iy);
            foreach (var ds in net.GetDownstream(cell.ix, cell.iy))
                upstreamCount[Idx(ds.ix, ds.iy)] += upstreamCount[idx];
        }

        // Find max for normalization
        int maxCount = 1;
        foreach (var cell in reachCells)
        {
            int c = upstreamCount[Idx(cell.ix, cell.iy)];
            if (c > maxCount) maxCount = c;
        }

        // Assign geometry: scale by log of relative upstream count
        foreach (var cell in reachCells)
        {
            int idx = Idx(cell.ix, cell.iy);
            int count = upstreamCount[idx];

            // Stream order proxy: log-scale growth
            double order = Math.Log(count) / Math.Log(Math.Max(maxCount, 2));
            double scale = Math.Pow(growthFactor, order * 4); // 4 approximate Strahler levels

            _width[idx] = headwaterWidth * scale;
            _depth[idx] = headwaterDepth * scale;
            _manningN[idx] = manningN;
        }
    }

    // ═══════════════════════════════════════════════════════════════
    //  Getters
    // ═══════════════════════════════════════════════════════════════

    /// <summary>Returns the channel width (m) at cell (ix, iy).</summary>
    public double GetWidth(int ix, int iy) => _width[Idx(ix, iy)];

    /// <summary>Returns the channel depth (m) at cell (ix, iy).</summary>
    public double GetDepth(int ix, int iy) => _depth[Idx(ix, iy)];

    /// <summary>Returns Manning's n at cell (ix, iy).</summary>
    public double GetManningN(int ix, int iy) => _manningN[Idx(ix, iy)];

    /// <summary>
    /// Computes the local bed slope S₀ at cell (ix, iy) as the
    /// elevation drop per unit reach length to the first downstream cell.
    /// Returns 0 if the cell has no downstream neighbour.
    /// </summary>
    public double GetBedSlope(int ix, int iy, RiverNetwork net, TerrainGrid terrain)
    {
        if (net == null || terrain == null) return 0;

        var downstream = net.GetDownstream(ix, iy);
        if (downstream.Count == 0) return 0;

        var ds = downstream[0];
        double elevHere = terrain[ix, iy];
        double elevDown = terrain[ds.ix, ds.iy];
        double drop = elevHere - elevDown;
        if (drop <= 0) return 0;

        // Distance between cell centres
        double dx = (ds.ix - ix) * Grid.Step;
        double dy = (ds.iy - iy) * Grid.Step;
        double dist = Math.Sqrt(dx * dx + dy * dy);
        if (dist <= 0) return 0;

        return drop / dist;
    }

    /// <summary>
    /// Computes the cross-sectional mean velocity (m/s) at cell (ix, iy)
    /// using Manning's equation with the local channel properties and bed slope.
    /// </summary>
    public double GetVelocity(int ix, int iy, RiverNetwork net, TerrainGrid terrain)
    {
        double n = GetManningN(ix, iy);
        double w = GetWidth(ix, iy);
        double d = GetDepth(ix, iy);
        double slope = GetBedSlope(ix, iy, net, terrain);

        double Rh = ManningEquation.RectangularHydraulicRadius(w, d);
        return ManningEquation.Velocity(n, Rh, slope);
    }
}
