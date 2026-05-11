using System;

namespace CSharpNumerics.Engines.Game.Terrain;

/// <summary>
/// Computes wind deflection over terrain using a simple boundary-layer model.
///
/// Wind is accelerated over ridges (speed-up) and decelerated in valleys (sheltering).
/// The model modifies a 2D wind field based on terrain slope in the wind direction.
///
/// Based on the Jackson-Hunt theory for flow over low hills:
///   speed-up factor ≈ 1 + Δs * (H / L)
/// where H is the hill height and L is the half-length.
///
/// For game purposes, this uses a simplified local slope model.
/// </summary>
public static class WindOverTerrain
{
    /// <summary>
    /// Modify a 2D wind velocity field based on terrain elevation.
    /// Wind accelerates going uphill (speed-up) and decelerates in the lee (sheltering).
    /// </summary>
    /// <param name="windX">X-component of wind field (modified in place).</param>
    /// <param name="windY">Y-component of wind field (modified in place).</param>
    /// <param name="terrain">Terrain collider for height queries.</param>
    /// <param name="gridNx">Wind grid X resolution.</param>
    /// <param name="gridNy">Wind grid Y resolution.</param>
    /// <param name="gridOriginX">World X of wind grid origin.</param>
    /// <param name="gridOriginY">World Y of wind grid origin.</param>
    /// <param name="gridCellSize">Wind grid cell size.</param>
    /// <param name="speedUpFactor">Multiplier for terrain-induced speed-up. Default 1.5.</param>
    public static void ApplyTerrainEffect(
        double[] windX, double[] windY,
        TerrainCollider terrain,
        int gridNx, int gridNy,
        double gridOriginX, double gridOriginY,
        double gridCellSize,
        double speedUpFactor = 1.5)
    {
        double eps = gridCellSize;

        for (int j = 1; j < gridNy - 1; j++)
        {
            for (int i = 1; i < gridNx - 1; i++)
            {
                int idx = i + j * gridNx;
                double wx = gridOriginX + i * gridCellSize;
                double wy = gridOriginY + j * gridCellSize;

                double h = terrain.SampleHeight(wx, wy);
                double hE = terrain.SampleHeight(wx + eps, wy);
                double hW = terrain.SampleHeight(wx - eps, wy);
                double hN = terrain.SampleHeight(wx, wy + eps);
                double hS = terrain.SampleHeight(wx, wy - eps);

                // Local terrain slopes
                double slopeX = (hE - hW) / (2 * eps);
                double slopeY = (hN - hS) / (2 * eps);

                // Wind speed-up: wind in the direction of upslope gets faster
                double ux = windX[idx];
                double uy = windY[idx];

                // Dot product of wind with slope gives uphill component
                // Positive = wind going uphill → speed up
                // Negative = wind going downhill → slow down (sheltering)
                double uphillComponent = ux * slopeX + uy * slopeY;

                // Speed-up factor: 1 + factor * |slope| for uphill
                // Sheltering: 1 / (1 + factor * |slope|) for downhill
                double slopeMag = Math.Sqrt(slopeX * slopeX + slopeY * slopeY);
                if (slopeMag < 1e-6) continue;

                double modifier;
                if (uphillComponent > 0)
                {
                    // Going uphill: accelerate
                    modifier = 1.0 + speedUpFactor * slopeMag;
                }
                else
                {
                    // Going downhill / lee side: decelerate
                    modifier = 1.0 / (1.0 + speedUpFactor * slopeMag);
                }

                windX[idx] *= modifier;
                windY[idx] *= modifier;

                // Add vertical deflection component (wind pushed upward over ridges)
                // This is stored as a modification to the horizontal wind magnitude
                // since we're working with a 2D field
            }
        }
    }

    /// <summary>
    /// Apply terrain blocking: zero out wind in cells below the terrain surface.
    /// This prevents wind from flowing through solid terrain.
    /// </summary>
    /// <param name="windX">X-component of wind field.</param>
    /// <param name="windY">Y-component of wind field.</param>
    /// <param name="terrain">Terrain collider.</param>
    /// <param name="gridNx">Wind grid X resolution.</param>
    /// <param name="gridNy">Wind grid Y resolution.</param>
    /// <param name="gridOriginX">World X of wind grid origin.</param>
    /// <param name="gridOriginY">World Y of wind grid origin.</param>
    /// <param name="gridCellSize">Wind grid cell size.</param>
    /// <param name="windHeight">Height above ground at which wind field is evaluated.</param>
    public static void ApplyTerrainBlocking(
        double[] windX, double[] windY,
        TerrainCollider terrain,
        int gridNx, int gridNy,
        double gridOriginX, double gridOriginY,
        double gridCellSize,
        double windHeight = 10.0)
    {
        for (int j = 0; j < gridNy; j++)
        {
            for (int i = 0; i < gridNx; i++)
            {
                int idx = i + j * gridNx;
                double wx = gridOriginX + i * gridCellSize;
                double wy = gridOriginY + j * gridCellSize;

                double terrainH = terrain.SampleHeight(wx, wy);
                if (windHeight < terrainH)
                {
                    windX[idx] = 0;
                    windY[idx] = 0;
                }
            }
        }
    }
}
