using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Engines.GIS.Terrain;

/// <summary>
/// A 2-D elevation surface layered on top of a <see cref="GeoGrid"/>.
/// Provides slope, aspect, and directional slope queries needed by
/// fire and terrain-spread models.
/// <para>
/// The grid is treated as a ground plane (Nz=1). One elevation value
/// per (ix, iy) cell.
/// </para>
/// </summary>
public class TerrainGrid
{
    private readonly double[] _elevation;

    /// <summary>The underlying spatial grid.</summary>
    public GeoGrid Grid { get; }

    /// <summary>Number of cells along x.</summary>
    public int Nx => Grid.Nx;

    /// <summary>Number of cells along y.</summary>
    public int Ny => Grid.Ny;

    /// <summary>Cell spacing in metres.</summary>
    public double Step => Grid.Step;

    private TerrainGrid(GeoGrid grid, double[] elevation)
    {
        Grid = grid ?? throw new ArgumentNullException(nameof(grid));
        _elevation = elevation ?? throw new ArgumentNullException(nameof(elevation));
    }

    /// <summary>
    /// Gets or sets the elevation at cell (ix, iy) in metres.
    /// </summary>
    public double this[int ix, int iy]
    {
        get => _elevation[iy * Nx + ix];
        set => _elevation[iy * Nx + ix] = value;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Factory methods
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Creates a terrain grid from a procedural elevation function f(x, y) → z.
    /// </summary>
    public static TerrainGrid FromFunction(GeoGrid grid, Func<double, double, double> elevationFn)
    {
        if (elevationFn == null) throw new ArgumentNullException(nameof(elevationFn));

        var elevation = new double[grid.Nx * grid.Ny];
        for (int iy = 0; iy < grid.Ny; iy++)
        {
            for (int ix = 0; ix < grid.Nx; ix++)
            {
                var centre = grid.CellCentre(ix, iy, 0);
                elevation[iy * grid.Nx + ix] = elevationFn(centre.x, centre.y);
            }
        }

        return new TerrainGrid(grid, elevation);
    }

    /// <summary>
    /// Creates a terrain grid from a 2-D array of elevation values.
    /// Array layout: row = iy (south→north), col = ix (west→east).
    /// </summary>
    public static TerrainGrid FromArray(GeoGrid grid, double[,] elevation)
    {
        if (elevation == null) throw new ArgumentNullException(nameof(elevation));
        if (elevation.GetLength(0) != grid.Ny || elevation.GetLength(1) != grid.Nx)
            throw new ArgumentException(
                $"Elevation array dimensions [{elevation.GetLength(0)},{elevation.GetLength(1)}] " +
                $"do not match grid [Ny={grid.Ny}, Nx={grid.Nx}].");

        var flat = new double[grid.Nx * grid.Ny];
        for (int iy = 0; iy < grid.Ny; iy++)
            for (int ix = 0; ix < grid.Nx; ix++)
                flat[iy * grid.Nx + ix] = elevation[iy, ix];

        return new TerrainGrid(grid, flat);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Slope & Aspect
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Terrain slope at cell (ix, iy) in radians.
    /// Uses central-difference gradient:
    /// tan(θ) = √((∂z/∂x)² + (∂z/∂y)²)
    /// At grid boundaries, uses one-sided differences.
    /// </summary>
    public double Slope(int ix, int iy)
    {
        double dzdx = GradientX(ix, iy);
        double dzdy = GradientY(ix, iy);
        return Math.Atan(Math.Sqrt(dzdx * dzdx + dzdy * dzdy));
    }

    /// <summary>
    /// Downslope aspect direction at cell (ix, iy) in radians.
    /// Convention: 0 = North (+y), π/2 = East (+x), π = South (−y), 3π/2 = West (−x).
    /// Returns 0 for flat terrain (no gradient).
    /// </summary>
    public double Aspect(int ix, int iy)
    {
        double dzdx = GradientX(ix, iy);
        double dzdy = GradientY(ix, iy);

        if (Math.Abs(dzdx) < 1e-15 && Math.Abs(dzdy) < 1e-15)
            return 0;

        // Downslope direction: negate the gradient
        // atan2 gives angle from +x axis; rotate so 0=North (+y direction)
        // aspect = atan2(-dzdx, -dzdy) in standard geographic convention (0=N, CW)
        double angle = Math.Atan2(-dzdx, -dzdy);
        if (angle < 0) angle += 2 * Math.PI;
        return angle;
    }

    /// <summary>
    /// Slope component along a specific direction at cell (ix, iy).
    /// Returns the slope angle (radians) in the given heading direction.
    /// Positive slope means uphill in that direction.
    /// </summary>
    /// <param name="ix">Cell x-index.</param>
    /// <param name="iy">Cell y-index.</param>
    /// <param name="direction">Horizontal direction vector (only x, y components used).</param>
    public double SlopeInDirection(int ix, int iy, Vector direction)
    {
        double dx = direction.x;
        double dy = direction.y;
        double len = Math.Sqrt(dx * dx + dy * dy);
        if (len < 1e-15) return 0;

        // Unit direction
        double ux = dx / len;
        double uy = dy / len;

        // Directional derivative dz/ds = (∂z/∂x)·ux + (∂z/∂y)·uy
        double dzds = GradientX(ix, iy) * ux + GradientY(ix, iy) * uy;

        // Slope angle in this direction: atan(dz/ds)
        return Math.Atan(dzds);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Gradient helpers (central differences, one-sided at edges)
    // ═══════════════════════════════════════════════════════════════

    private double GradientX(int ix, int iy)
    {
        if (Nx < 2) return 0;
        if (ix == 0)
            return (this[ix + 1, iy] - this[ix, iy]) / Step;
        if (ix == Nx - 1)
            return (this[ix, iy] - this[ix - 1, iy]) / Step;
        return (this[ix + 1, iy] - this[ix - 1, iy]) / (2 * Step);
    }

    private double GradientY(int ix, int iy)
    {
        if (Ny < 2) return 0;
        if (iy == 0)
            return (this[ix, iy + 1] - this[ix, iy]) / Step;
        if (iy == Ny - 1)
            return (this[ix, iy] - this[ix, iy - 1]) / Step;
        return (this[ix, iy + 1] - this[ix, iy - 1]) / (2 * Step);
    }
}
