using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Engines.Game.Fluids;

/// <summary>
/// Represents a source that injects velocity and/or density into the fluid field.
/// Used for exhaust plumes, wind sources, explosions, smoke generators, etc.
/// </summary>
public class FluidEmitter
{
    /// <summary>Position in grid coordinates (fractional).</summary>
    public Vector Position { get; set; }

    /// <summary>Radius of influence in grid cells.</summary>
    public double Radius { get; set; } = 1.5;

    /// <summary>Density injection rate per second.</summary>
    public double DensityRate { get; set; }

    /// <summary>Temperature injection (for buoyancy-driven flows).</summary>
    public double Temperature { get; set; } = 500.0;

    /// <summary>Velocity injection direction and magnitude (grid cells / second).</summary>
    public Vector Velocity { get; set; } = new Vector(0, 0, 0);

    /// <summary>Whether this emitter is active.</summary>
    public bool Active { get; set; } = true;

    /// <summary>
    /// Creates an emitter at the given grid position.
    /// </summary>
    public FluidEmitter(Vector position, double densityRate = 5.0, double radius = 1.5)
    {
        Position = position;
        DensityRate = densityRate;
        Radius = radius;
    }

    /// <summary>
    /// Injects density and velocity into a 2D fluid field for one time step.
    /// </summary>
    public void Apply2D(double[] density, double[] u, double[] v, int nx, int ny, double dt)
    {
        if (!Active) return;

        int cx = (int)Position.x;
        int cy = (int)Position.y;
        int r = (int)Math.Ceiling(Radius);

        for (int j = Math.Max(1, cy - r); j <= Math.Min(ny - 2, cy + r); j++)
        {
            for (int i = Math.Max(1, cx - r); i <= Math.Min(nx - 2, cx + r); i++)
            {
                double dx = i - Position.x;
                double dy = j - Position.y;
                double dist = Math.Sqrt(dx * dx + dy * dy);
                if (dist > Radius) continue;

                double falloff = 1.0 - dist / Radius;
                int idx = i + j * nx;

                density[idx] += DensityRate * falloff * dt;
                u[idx] += Velocity.x * falloff * dt;
                v[idx] += Velocity.y * falloff * dt;
            }
        }
    }

    /// <summary>
    /// Injects density and velocity into a 3D fluid field for one time step.
    /// </summary>
    public void Apply3D(double[] density, double[] u, double[] v, double[] w,
        int nx, int ny, int nz, double dt)
    {
        if (!Active) return;

        int cx = (int)Position.x;
        int cy = (int)Position.y;
        int cz = (int)Position.z;
        int r = (int)Math.Ceiling(Radius);

        for (int k = Math.Max(1, cz - r); k <= Math.Min(nz - 2, cz + r); k++)
        {
            for (int j = Math.Max(1, cy - r); j <= Math.Min(ny - 2, cy + r); j++)
            {
                for (int i = Math.Max(1, cx - r); i <= Math.Min(nx - 2, cx + r); i++)
                {
                    double dx = i - Position.x;
                    double dy = j - Position.y;
                    double dz = k - Position.z;
                    double dist = Math.Sqrt(dx * dx + dy * dy + dz * dz);
                    if (dist > Radius) continue;

                    double falloff = 1.0 - dist / Radius;
                    int idx = i + j * nx + k * nx * ny;

                    density[idx] += DensityRate * falloff * dt;
                    u[idx] += Velocity.x * falloff * dt;
                    v[idx] += Velocity.y * falloff * dt;
                    w[idx] += Velocity.z * falloff * dt;
                }
            }
        }
    }
}
