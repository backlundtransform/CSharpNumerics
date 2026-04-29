using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Numerics.FiniteDifference;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Multiphysics.Solvers;

/// <summary>
/// Solves 2D plane-stress linear elasticity using an iterative finite-difference approach.
/// Equilibrium equations (body forces neglected, point loads applied):
///   ∂σxx/∂x + ∂τxy/∂y = 0
///   ∂τxy/∂x + ∂σyy/∂y = 0
/// with constitutive law (plane stress):
///   σxx = E/(1−ν²)(εxx + ν·εyy)
///   σyy = E/(1−ν²)(εyy + ν·εxx)
///   τxy  = E/(2(1+ν))·γxy
/// Solved iteratively for displacements (ux, uy) using a relaxation method.
/// </summary>
internal class PlaneStressSolver : IMultiphysicsSolver
{
    public SimulationResult Solve(SimulationBuilder cfg)
    {
        var mat = cfg.Material.Value;
        double E = mat.YoungsModulus;
        double nu = mat.PoissonsRatio;

        double dx = cfg.GeomWidth / cfg.Nx;
        double dy = cfg.GeomHeight / cfg.Ny;
        var grid = new Grid2D(cfg.Nx, cfg.Ny, dx, dy);
        int nx = cfg.Nx, ny = cfg.Ny;

        // Plane-stress stiffness coefficients
        double c11 = E / (1.0 - nu * nu);          // σxx per εxx
        double c12 = nu * c11;                      // coupling
        double c33 = E / (2.0 * (1.0 + nu));       // shear modulus G

        // Displacement fields
        var ux = new double[nx, ny];
        var uy = new double[nx, ny];

        // External force field (from point sources)
        // Sources2D: (ix, iy, value) — value is force magnitude applied in x-direction
        // For y-forces, use negative iy encoding or separate. Here, assume Sources2D is Fx.
        var fx = new double[nx, ny];
        var fy = new double[nx, ny];
        foreach (var (ix, iy, value) in cfg.Sources2D)
        {
            if (ix >= 0 && ix < nx && iy >= 0 && iy < ny)
                fx[ix, iy] = value;
        }

        // Apply body force from UniformLoad (in y-direction, like gravity/distributed)
        if (cfg.UniformLoad != 0)
        {
            for (int iy = 0; iy < ny; iy++)
                for (int ix = 0; ix < nx; ix++)
                    fy[ix, iy] = cfg.UniformLoad;
        }

        // Boundary conditions:
        // Left boundary (ix=0): fixed (ux=uy=0) — Dirichlet
        // Right boundary: traction from RightBC (applied force density)
        // Top/Bottom: free or constrained based on TopBC/BottomBC
        // LeftBC: displacement constraint (0 = fixed)

        // Iterative Gauss-Seidel relaxation for the Navier-Cauchy equations
        double invDx2 = 1.0 / (dx * dx);
        double invDy2 = 1.0 / (dy * dy);
        double inv4DxDy = 1.0 / (4.0 * dx * dy);
        double omega = 1.4; // SOR relaxation factor

        int iterations = 0;
        double residual = double.MaxValue;

        for (int iter = 0; iter < cfg.MaxIterations && residual > cfg.Tolerance; iter++)
        {
            residual = 0;
            iterations = iter + 1;

            for (int iy = 1; iy < ny - 1; iy++)
            {
                for (int ix = 1; ix < nx - 1; ix++)
                {
                    // Skip left boundary (fixed)
                    if (ix == 0) continue;

                    // Navier-Cauchy eq 1 (x-equilibrium):
                    // c11·∂²ux/∂x² + c33·∂²ux/∂y² + (c12+c33)·∂²uy/∂x∂y + fx = 0
                    double uxOld = ux[ix, iy];
                    double d2uxdx2 = (ux[ix + 1, iy] - 2.0 * ux[ix, iy] + ux[ix - 1, iy]) * invDx2;
                    double d2uxdy2 = (ux[ix, iy + 1] - 2.0 * ux[ix, iy] + ux[ix, iy - 1]) * invDy2;
                    double d2uydxdy = (uy[ix + 1, iy + 1] - uy[ix - 1, iy + 1]
                                     - uy[ix + 1, iy - 1] + uy[ix - 1, iy - 1]) * inv4DxDy;

                    double rhs_x = -(c12 + c33) * d2uydxdy - fx[ix, iy];
                    double denom_x = -2.0 * c11 * invDx2 - 2.0 * c33 * invDy2;
                    double numX = rhs_x - c11 * (ux[ix + 1, iy] + ux[ix - 1, iy]) * invDx2
                                        - c33 * (ux[ix, iy + 1] + ux[ix, iy - 1]) * invDy2;
                    double uxNew = numX / denom_x;
                    ux[ix, iy] = uxOld + omega * (uxNew - uxOld);

                    // Navier-Cauchy eq 2 (y-equilibrium):
                    // c33·∂²uy/∂x² + c11·∂²uy/∂y² + (c12+c33)·∂²ux/∂x∂y + fy = 0
                    double uyOld = uy[ix, iy];
                    double d2uydx2 = (uy[ix + 1, iy] - 2.0 * uy[ix, iy] + uy[ix - 1, iy]) * invDx2;
                    double d2uydy2 = (uy[ix, iy + 1] - 2.0 * uy[ix, iy] + uy[ix, iy - 1]) * invDy2;
                    double d2uxdxdy = (ux[ix + 1, iy + 1] - ux[ix - 1, iy + 1]
                                     - ux[ix + 1, iy - 1] + ux[ix - 1, iy - 1]) * inv4DxDy;

                    double rhs_y = -(c12 + c33) * d2uxdxdy - fy[ix, iy];
                    double denom_y = -2.0 * c33 * invDx2 - 2.0 * c11 * invDy2;
                    double numY = rhs_y - c33 * (uy[ix + 1, iy] + uy[ix - 1, iy]) * invDx2
                                        - c11 * (uy[ix, iy + 1] + uy[ix, iy - 1]) * invDy2;
                    double uyNew = numY / denom_y;
                    uy[ix, iy] = uyOld + omega * (uyNew - uyOld);

                    residual += (ux[ix, iy] - uxOld) * (ux[ix, iy] - uxOld)
                              + (uy[ix, iy] - uyOld) * (uy[ix, iy] - uyOld);
                }
            }

            // Apply BCs: left = fixed
            for (int iy = 0; iy < ny; iy++)
            {
                ux[0, iy] = 0;
                uy[0, iy] = 0;
            }

            // Right boundary: Neumann (free or traction)
            for (int iy = 0; iy < ny; iy++)
            {
                ux[nx - 1, iy] = ux[nx - 2, iy];
                uy[nx - 1, iy] = uy[nx - 2, iy];
            }

            // Bottom boundary: fixed in y (roller)
            for (int ix = 0; ix < nx; ix++)
            {
                uy[ix, 0] = 0;
                ux[ix, 0] = ux[ix, 1]; // free in x (roller)
            }

            // Top boundary: free
            for (int ix = 0; ix < nx; ix++)
            {
                ux[ix, ny - 1] = ux[ix, ny - 2];
                uy[ix, ny - 1] = uy[ix, ny - 2];
            }

            residual = Math.Sqrt(residual / (nx * ny));
        }

        // Compute stress fields from final displacement
        var sxx = new double[nx, ny];
        var syy = new double[nx, ny];
        var sxy = new double[nx, ny];

        for (int iy = 1; iy < ny - 1; iy++)
        {
            for (int ix = 1; ix < nx - 1; ix++)
            {
                double exx = (ux[ix + 1, iy] - ux[ix - 1, iy]) / (2.0 * dx);
                double eyy = (uy[ix, iy + 1] - uy[ix, iy - 1]) / (2.0 * dy);
                double gxy = (ux[ix, iy + 1] - ux[ix, iy - 1]) / (2.0 * dy)
                           + (uy[ix + 1, iy] - uy[ix - 1, iy]) / (2.0 * dx);

                sxx[ix, iy] = c11 * exx + c12 * eyy;
                syy[ix, iy] = c12 * exx + c11 * eyy;
                sxy[ix, iy] = c33 * gxy;
            }
        }

        // Primary field: von Mises stress for min/max
        double[,] vonMises = new double[nx, ny];
        double minV = double.MaxValue, maxV = double.MinValue;
        for (int iy = 0; iy < ny; iy++)
            for (int ix = 0; ix < nx; ix++)
            {
                double vm = Math.Sqrt(sxx[ix, iy] * sxx[ix, iy]
                                    - sxx[ix, iy] * syy[ix, iy]
                                    + syy[ix, iy] * syy[ix, iy]
                                    + 3.0 * sxy[ix, iy] * sxy[ix, iy]);
                vonMises[ix, iy] = vm;
                if (vm < minV) minV = vm;
                if (vm > maxV) maxV = vm;
            }

        return new SimulationResult
        {
            Type = MultiphysicsType.PlaneStress,
            Field = vonMises,
            Ux = ux,
            Uy = uy,
            StressXX = sxx,
            StressYY = syy,
            StressXY = sxy,
            MaxValue = maxV,
            MinValue = minV,
            Iterations = iterations,
            FinalTime = residual
        };
    }
}
