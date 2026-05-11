using System;

namespace CSharpNumerics.Engines.Game.Fluids;

/// <summary>
/// Represents a static or dynamic obstacle in the fluid domain.
/// Obstacles block flow by zeroing velocity inside their region and
/// enforcing no-penetration boundary conditions on adjacent cells.
/// </summary>
public class FluidObstacle
{
    /// <summary>Minimum corner of the axis-aligned obstacle in grid coordinates.</summary>
    public int MinX { get; set; }
    /// <summary>Minimum Y.</summary>
    public int MinY { get; set; }
    /// <summary>Minimum Z (3D only).</summary>
    public int MinZ { get; set; }

    /// <summary>Maximum corner X (inclusive).</summary>
    public int MaxX { get; set; }
    /// <summary>Maximum Y (inclusive).</summary>
    public int MaxY { get; set; }
    /// <summary>Maximum Z (inclusive, 3D only).</summary>
    public int MaxZ { get; set; }

    /// <summary>Whether this obstacle is active.</summary>
    public bool Active { get; set; } = true;

    /// <summary>
    /// Creates a 2D box obstacle in grid coordinates.
    /// </summary>
    public static FluidObstacle Box2D(int minX, int minY, int maxX, int maxY)
    {
        return new FluidObstacle { MinX = minX, MinY = minY, MaxX = maxX, MaxY = maxY };
    }

    /// <summary>
    /// Creates a 3D box obstacle in grid coordinates.
    /// </summary>
    public static FluidObstacle Box3D(int minX, int minY, int minZ, int maxX, int maxY, int maxZ)
    {
        return new FluidObstacle
        {
            MinX = minX, MinY = minY, MinZ = minZ,
            MaxX = maxX, MaxY = maxY, MaxZ = maxZ
        };
    }

    /// <summary>
    /// Applies obstacle boundary conditions to a 2D velocity field.
    /// Zeros velocity inside the obstacle and reflects velocities at boundaries.
    /// </summary>
    public void Apply2D(double[] u, double[] v, int nx, int ny)
    {
        if (!Active) return;

        for (int j = Math.Max(0, MinY); j <= Math.Min(ny - 1, MaxY); j++)
        {
            for (int i = Math.Max(0, MinX); i <= Math.Min(nx - 1, MaxX); i++)
            {
                int idx = i + j * nx;
                u[idx] = 0;
                v[idx] = 0;
            }
        }

        // Reflect normal velocity at obstacle boundaries
        // Left face
        if (MinX > 0)
        {
            for (int j = Math.Max(0, MinY); j <= Math.Min(ny - 1, MaxY); j++)
            {
                int idx = (MinX - 1) + j * nx;
                if (u[idx] > 0) u[idx] = 0;
            }
        }
        // Right face
        if (MaxX < nx - 1)
        {
            for (int j = Math.Max(0, MinY); j <= Math.Min(ny - 1, MaxY); j++)
            {
                int idx = (MaxX + 1) + j * nx;
                if (u[idx] < 0) u[idx] = 0;
            }
        }
        // Bottom face
        if (MinY > 0)
        {
            for (int i = Math.Max(0, MinX); i <= Math.Min(nx - 1, MaxX); i++)
            {
                int idx = i + (MinY - 1) * nx;
                if (v[idx] > 0) v[idx] = 0;
            }
        }
        // Top face
        if (MaxY < ny - 1)
        {
            for (int i = Math.Max(0, MinX); i <= Math.Min(nx - 1, MaxX); i++)
            {
                int idx = i + (MaxY + 1) * nx;
                if (v[idx] < 0) v[idx] = 0;
            }
        }
    }

    /// <summary>
    /// Applies obstacle boundary conditions to a 3D velocity field.
    /// </summary>
    public void Apply3D(double[] u, double[] v, double[] w, int nx, int ny, int nz)
    {
        if (!Active) return;

        for (int k = Math.Max(0, MinZ); k <= Math.Min(nz - 1, MaxZ); k++)
        {
            for (int j = Math.Max(0, MinY); j <= Math.Min(ny - 1, MaxY); j++)
            {
                for (int i = Math.Max(0, MinX); i <= Math.Min(nx - 1, MaxX); i++)
                {
                    int idx = i + j * nx + k * nx * ny;
                    u[idx] = 0;
                    v[idx] = 0;
                    w[idx] = 0;
                }
            }
        }
    }
}
