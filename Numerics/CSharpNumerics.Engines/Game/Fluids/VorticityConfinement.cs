using System;

namespace CSharpNumerics.Engines.Game.Fluids;

/// <summary>
/// Vorticity confinement adds back small-scale rotational energy lost to numerical diffusion.
/// This sharpens smoke curls and makes turbulent features persist longer — essential for
/// visually compelling fluid simulations in games.
/// 
/// Algorithm (Steinhoff &amp; Underhill):
///   1. Compute vorticity ω at each cell
///   2. Compute the gradient of |ω|: η = ∇|ω|
///   3. Normalize: N = η / |η|
///   4. Apply confinement force: f = ε · h · (N × ω)
/// </summary>
public static class VorticityConfinement
{
    /// <summary>
    /// Applies vorticity confinement to a 2D velocity field.
    /// In 2D, vorticity is a scalar: ω = dv/dx − du/dy.
    /// Confinement force pushes velocity toward vortex centres.
    /// </summary>
    /// <param name="u">X-velocity field (nx × ny).</param>
    /// <param name="v">Y-velocity field (nx × ny).</param>
    /// <param name="nx">Grid width.</param>
    /// <param name="ny">Grid height.</param>
    /// <param name="epsilon">Confinement strength (typical 0.1–5.0).</param>
    /// <param name="dt">Time step.</param>
    public static void Apply2D(double[] u, double[] v, int nx, int ny, double epsilon, double dt)
    {
        if (epsilon <= 0) return;

        int size = nx * ny;
        var omega = new double[size];     // vorticity scalar
        var absOmega = new double[size];  // |ω|

        // Step 1: compute vorticity ω = dv/dx - du/dy
        for (int j = 1; j < ny - 1; j++)
        {
            for (int i = 1; i < nx - 1; i++)
            {
                int idx = i + j * nx;
                double dvdx = (v[idx + 1] - v[idx - 1]) * 0.5;
                double dudy = (u[idx + nx] - u[idx - nx]) * 0.5;
                omega[idx] = dvdx - dudy;
                absOmega[idx] = Math.Abs(omega[idx]);
            }
        }

        // Step 2-4: compute gradient of |ω|, normalize, apply force
        for (int j = 2; j < ny - 2; j++)
        {
            for (int i = 2; i < nx - 2; i++)
            {
                int idx = i + j * nx;

                // Gradient of |ω|
                double dAbsdx = (absOmega[idx + 1] - absOmega[idx - 1]) * 0.5;
                double dAbsdy = (absOmega[idx + nx] - absOmega[idx - nx]) * 0.5;

                double len = Math.Sqrt(dAbsdx * dAbsdx + dAbsdy * dAbsdy) + 1e-10;
                double nx2 = dAbsdx / len;
                double ny2 = dAbsdy / len;

                // In 2D, N × ω gives force components:
                // f_x = ε · (N_y · ω), f_y = ε · (−N_x · ω)
                u[idx] += epsilon * dt * ny2 * omega[idx];
                v[idx] -= epsilon * dt * nx2 * omega[idx];
            }
        }
    }

    /// <summary>
    /// Applies vorticity confinement to a 3D velocity field.
    /// </summary>
    /// <param name="u">X-velocity (nx × ny × nz).</param>
    /// <param name="v">Y-velocity (nx × ny × nz).</param>
    /// <param name="w">Z-velocity (nx × ny × nz).</param>
    /// <param name="nx">Grid X size.</param>
    /// <param name="ny">Grid Y size.</param>
    /// <param name="nz">Grid Z size.</param>
    /// <param name="epsilon">Confinement strength.</param>
    /// <param name="dt">Time step.</param>
    public static void Apply3D(double[] u, double[] v, double[] w,
        int nx, int ny, int nz, double epsilon, double dt)
    {
        if (epsilon <= 0) return;

        int nxy = nx * ny;
        int size = nxy * nz;
        var omegaX = new double[size];
        var omegaY = new double[size];
        var omegaZ = new double[size];
        var absOmega = new double[size];

        // Step 1: compute vorticity ω = ∇ × v
        for (int k = 1; k < nz - 1; k++)
        {
            for (int j = 1; j < ny - 1; j++)
            {
                for (int i = 1; i < nx - 1; i++)
                {
                    int idx = i + j * nx + k * nxy;

                    double dwdy = (w[idx + nx] - w[idx - nx]) * 0.5;
                    double dvdz = (v[idx + nxy] - v[idx - nxy]) * 0.5;
                    omegaX[idx] = dwdy - dvdz;

                    double dudz = (u[idx + nxy] - u[idx - nxy]) * 0.5;
                    double dwdx = (w[idx + 1] - w[idx - 1]) * 0.5;
                    omegaY[idx] = dudz - dwdx;

                    double dvdx = (v[idx + 1] - v[idx - 1]) * 0.5;
                    double dudy = (u[idx + nx] - u[idx - nx]) * 0.5;
                    omegaZ[idx] = dvdx - dudy;

                    absOmega[idx] = Math.Sqrt(
                        omegaX[idx] * omegaX[idx] +
                        omegaY[idx] * omegaY[idx] +
                        omegaZ[idx] * omegaZ[idx]);
                }
            }
        }

        // Steps 2-4: gradient of |ω|, normalize, apply N × ω force
        for (int k = 2; k < nz - 2; k++)
        {
            for (int j = 2; j < ny - 2; j++)
            {
                for (int i = 2; i < nx - 2; i++)
                {
                    int idx = i + j * nx + k * nxy;

                    double dAdx = (absOmega[idx + 1] - absOmega[idx - 1]) * 0.5;
                    double dAdy = (absOmega[idx + nx] - absOmega[idx - nx]) * 0.5;
                    double dAdz = (absOmega[idx + nxy] - absOmega[idx - nxy]) * 0.5;

                    double len = Math.Sqrt(dAdx * dAdx + dAdy * dAdy + dAdz * dAdz) + 1e-10;
                    double Nx = dAdx / len;
                    double Ny = dAdy / len;
                    double Nz = dAdz / len;

                    // f = ε · (N × ω)
                    double fx = epsilon * (Ny * omegaZ[idx] - Nz * omegaY[idx]);
                    double fy = epsilon * (Nz * omegaX[idx] - Nx * omegaZ[idx]);
                    double fz = epsilon * (Nx * omegaY[idx] - Ny * omegaX[idx]);

                    u[idx] += fx * dt;
                    v[idx] += fy * dt;
                    w[idx] += fz * dt;
                }
            }
        }
    }
}
