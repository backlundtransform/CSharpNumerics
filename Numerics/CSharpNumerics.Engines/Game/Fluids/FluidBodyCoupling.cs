using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Engines.Game.Fluids;

/// <summary>
/// Two-way coupling between rigid bodies and the fluid field.
/// 
/// Fluid → Body: body experiences drag/lift from the local fluid velocity.
/// Body → Fluid: body displaces fluid, injecting velocity at its occupied cells.
/// </summary>
public static class FluidBodyCoupling
{
    /// <summary>
    /// Computes the drag force on a spherical body immersed in the fluid field.
    /// F_drag = −½ · ρ · Cd · A · |v_rel| · v_rel
    /// </summary>
    /// <param name="bodyVelocity">Body velocity in world coordinates (m/s).</param>
    /// <param name="fluidVelocity">Fluid velocity at body position (m/s).</param>
    /// <param name="dragCoefficient">Drag coefficient Cd (sphere ≈ 0.47).</param>
    /// <param name="crossSectionArea">Frontal cross-section area (m²).</param>
    /// <param name="fluidDensity">Fluid density (kg/m³).</param>
    /// <returns>Drag force vector in Newtons (opposes relative motion).</returns>
    public static Vector ComputeDragForce(
        Vector bodyVelocity,
        Vector fluidVelocity,
        double dragCoefficient,
        double crossSectionArea,
        double fluidDensity)
    {
        var vRel = bodyVelocity - fluidVelocity;
        double speed = vRel.GetMagnitude();
        if (speed < 1e-10) return new Vector(0, 0, 0);

        double forceMag = 0.5 * fluidDensity * dragCoefficient * crossSectionArea * speed * speed;
        return -(forceMag / speed) * vRel;
    }

    /// <summary>
    /// Injects body-displacement velocity into a 2D fluid field.
    /// Cells overlapping the body are set to the body's velocity,
    /// simulating the body pushing fluid out of the way.
    /// </summary>
    /// <param name="u">X-velocity field.</param>
    /// <param name="v">Y-velocity field.</param>
    /// <param name="nx">Grid X size.</param>
    /// <param name="ny">Grid Y size.</param>
    /// <param name="bodyGridX">Body centre X in grid coordinates.</param>
    /// <param name="bodyGridY">Body centre Y in grid coordinates.</param>
    /// <param name="bodyRadius">Body radius in grid cells.</param>
    /// <param name="bodyVelX">Body X-velocity in grid cells/s.</param>
    /// <param name="bodyVelY">Body Y-velocity in grid cells/s.</param>
    public static void DisplaceFluid2D(
        double[] u, double[] v, int nx, int ny,
        double bodyGridX, double bodyGridY, double bodyRadius,
        double bodyVelX, double bodyVelY)
    {
        int r = (int)Math.Ceiling(bodyRadius);
        int cx = (int)bodyGridX, cy = (int)bodyGridY;

        for (int j = Math.Max(1, cy - r); j <= Math.Min(ny - 2, cy + r); j++)
        {
            for (int i = Math.Max(1, cx - r); i <= Math.Min(nx - 2, cx + r); i++)
            {
                double dx = i - bodyGridX;
                double dy = j - bodyGridY;
                if (dx * dx + dy * dy > bodyRadius * bodyRadius) continue;

                int idx = i + j * nx;
                u[idx] = bodyVelX;
                v[idx] = bodyVelY;
            }
        }
    }

    /// <summary>
    /// Injects body-displacement velocity into a 3D fluid field.
    /// </summary>
    public static void DisplaceFluid3D(
        double[] u, double[] v, double[] w, int nx, int ny, int nz,
        double bodyGridX, double bodyGridY, double bodyGridZ, double bodyRadius,
        double bodyVelX, double bodyVelY, double bodyVelZ)
    {
        int nxy = nx * ny;
        int r = (int)Math.Ceiling(bodyRadius);
        int cx = (int)bodyGridX, cy = (int)bodyGridY, cz = (int)bodyGridZ;
        double r2 = bodyRadius * bodyRadius;

        for (int k = Math.Max(1, cz - r); k <= Math.Min(nz - 2, cz + r); k++)
            for (int j = Math.Max(1, cy - r); j <= Math.Min(ny - 2, cy + r); j++)
                for (int i = Math.Max(1, cx - r); i <= Math.Min(nx - 2, cx + r); i++)
                {
                    double dx = i - bodyGridX;
                    double dy = j - bodyGridY;
                    double dz = k - bodyGridZ;
                    if (dx * dx + dy * dy + dz * dz > r2) continue;

                    int idx = i + j * nx + k * nxy;
                    u[idx] = bodyVelX;
                    v[idx] = bodyVelY;
                    w[idx] = bodyVelZ;
                }
    }
}
