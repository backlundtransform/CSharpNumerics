using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Engines.Game.Terrain;

/// <summary>
/// Heightmap-based collision for game terrain.
/// Stores a 2D grid of elevation values and provides:
///   - Point-on-terrain height queries
///   - Surface normal estimation
///   - Collision resolution (push objects above terrain)
///
/// Separate from the GIS <c>TerrainGrid</c> — this is optimized for
/// game-engine collision and does not depend on <c>GeoGrid</c>.
/// </summary>
public class TerrainCollider
{
    private readonly double[] _heights;
    private readonly int _resX, _resY;
    private readonly double _originX, _originY;
    private readonly double _cellSize;

    /// <summary>Grid resolution X.</summary>
    public int ResX => _resX;

    /// <summary>Grid resolution Y.</summary>
    public int ResY => _resY;

    /// <summary>World X origin.</summary>
    public double OriginX => _originX;

    /// <summary>World Y origin.</summary>
    public double OriginY => _originY;

    /// <summary>Cell size in world units.</summary>
    public double CellSize => _cellSize;

    /// <summary>World extent in X.</summary>
    public double Width => (_resX - 1) * _cellSize;

    /// <summary>World extent in Y.</summary>
    public double Depth => (_resY - 1) * _cellSize;

    /// <summary>
    /// Creates a terrain collider from a heightmap.
    /// </summary>
    /// <param name="heights">Flat array of height values [resX * resY], row-major.</param>
    /// <param name="resX">Number of vertices along X.</param>
    /// <param name="resY">Number of vertices along Y.</param>
    /// <param name="cellSize">Distance between adjacent vertices.</param>
    /// <param name="originX">World X of the (0,0) vertex.</param>
    /// <param name="originY">World Y of the (0,0) vertex.</param>
    public TerrainCollider(double[] heights, int resX, int resY, double cellSize,
        double originX = 0, double originY = 0)
    {
        if (heights.Length != resX * resY)
            throw new ArgumentException("Heights array size must equal resX * resY.");
        _heights = (double[])heights.Clone();
        _resX = resX;
        _resY = resY;
        _cellSize = cellSize;
        _originX = originX;
        _originY = originY;
    }

    /// <summary>
    /// Creates a flat terrain at the given height.
    /// </summary>
    public static TerrainCollider CreateFlat(int resX, int resY, double cellSize,
        double height = 0, double originX = 0, double originY = 0)
    {
        var h = new double[resX * resY];
        for (int i = 0; i < h.Length; i++) h[i] = height;
        return new TerrainCollider(h, resX, resY, cellSize, originX, originY);
    }

    /// <summary>
    /// Creates terrain from a height function f(worldX, worldY) → height.
    /// </summary>
    public static TerrainCollider FromFunction(int resX, int resY, double cellSize,
        Func<double, double, double> heightFn, double originX = 0, double originY = 0)
    {
        var h = new double[resX * resY];
        for (int j = 0; j < resY; j++)
            for (int i = 0; i < resX; i++)
                h[i + j * resX] = heightFn(originX + i * cellSize, originY + j * cellSize);
        return new TerrainCollider(h, resX, resY, cellSize, originX, originY);
    }

    /// <summary>
    /// Get the height value at grid vertex (ix, iy).
    /// </summary>
    public double GetHeight(int ix, int iy) => _heights[ix + iy * _resX];

    /// <summary>
    /// Sample terrain height at world position (wx, wy) using bilinear interpolation.
    /// Returns 0 if outside the terrain bounds.
    /// </summary>
    public double SampleHeight(double wx, double wy)
    {
        double gx = (wx - _originX) / _cellSize;
        double gy = (wy - _originY) / _cellSize;

        if (gx < 0 || gx >= _resX - 1 || gy < 0 || gy >= _resY - 1) return 0;

        int ix = (int)gx;
        int iy = (int)gy;
        double fx = gx - ix;
        double fy = gy - iy;

        double h00 = _heights[ix + iy * _resX];
        double h10 = _heights[ix + 1 + iy * _resX];
        double h01 = _heights[ix + (iy + 1) * _resX];
        double h11 = _heights[ix + 1 + (iy + 1) * _resX];

        return (1 - fx) * (1 - fy) * h00 + fx * (1 - fy) * h10
             + (1 - fx) * fy * h01 + fx * fy * h11;
    }

    /// <summary>
    /// Estimate the surface normal at world position (wx, wy) using central differences.
    /// </summary>
    public Vector SurfaceNormal(double wx, double wy)
    {
        double eps = _cellSize * 0.5;
        double hL = SampleHeight(wx - eps, wy);
        double hR = SampleHeight(wx + eps, wy);
        double hD = SampleHeight(wx, wy - eps);
        double hU = SampleHeight(wx, wy + eps);

        double dhdx = (hR - hL) / (2 * eps);
        double dhdy = (hU - hD) / (2 * eps);

        // Normal = (-dh/dx, -dh/dy, 1), normalized
        var n = new Vector(-dhdx, -dhdy, 1);
        double mag = n.GetMagnitude();
        return (1.0 / mag) * n;
    }

    /// <summary>
    /// Test if a sphere at (position, radius) is colliding with the terrain.
    /// Returns true if the sphere bottom is below the terrain surface.
    /// </summary>
    /// <param name="position">Sphere center world position.</param>
    /// <param name="radius">Sphere radius.</param>
    /// <param name="penetration">Penetration depth (positive = overlapping).</param>
    /// <param name="normal">Contact normal (terrain surface normal).</param>
    public bool TestSphere(Vector position, double radius, out double penetration, out Vector normal)
    {
        double terrainH = SampleHeight(position.x, position.y);
        double bottom = position.z - radius;
        penetration = terrainH - bottom;
        normal = SurfaceNormal(position.x, position.y);

        return penetration > 0;
    }

    /// <summary>
    /// Resolve a sphere collision by pushing it above the terrain and reflecting velocity.
    /// </summary>
    /// <param name="position">Sphere position (modified in place).</param>
    /// <param name="velocity">Sphere velocity (modified in place).</param>
    /// <param name="radius">Sphere radius.</param>
    /// <param name="restitution">Bounce coefficient.</param>
    /// <returns>True if a collision was resolved.</returns>
    public bool ResolveSphere(ref Vector position, ref Vector velocity, double radius, double restitution = 0.3)
    {
        if (!TestSphere(position, radius, out double pen, out Vector normal))
            return false;

        // Push out
        position = position + pen * normal;

        // Reflect velocity component along normal
        double vn = velocity.Dot(normal);
        if (vn < 0)
            velocity = velocity - (1.0 + restitution) * vn * normal;

        return true;
    }
}
