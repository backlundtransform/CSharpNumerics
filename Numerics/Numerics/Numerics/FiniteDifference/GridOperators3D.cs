using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Numerics.FiniteDifference;

/// <summary>
/// Discrete finite-difference operators on a <see cref="Grid3D"/>.
/// 3D extension of <see cref="GridOperators"/>: Laplacian, gradient, divergence,
/// advection, and Poisson solver using the 7-point stencil.
/// </summary>
public static class GridOperators3D
{
    // ──────────────────────────────────────────────
    //  Laplacian: ∇²u = ∂²u/∂x² + ∂²u/∂y² + ∂²u/∂z²
    // ──────────────────────────────────────────────

    /// <summary>
    /// 3D Laplacian using the 7-point stencil.
    /// </summary>
    public static VectorN Laplacian3D(VectorN u, Grid3D grid,
        BoundaryCondition bc = BoundaryCondition.Dirichlet)
    {
        int nx = grid.Nx, ny = grid.Ny, nz = grid.Nz;
        var result = new double[grid.Length];
        double invDx2 = 1.0 / (grid.Dx * grid.Dx);
        double invDy2 = 1.0 / (grid.Dy * grid.Dy);
        double invDz2 = 1.0 / (grid.Dz * grid.Dz);

        for (int iz = 0; iz < nz; iz++)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                for (int ix = 0; ix < nx; ix++)
                {
                    int idx = grid.Index(ix, iy, iz);
                    double center = u[idx];

                    double xp = GetValue(u, ix + 1, iy, iz, grid, bc);
                    double xm = GetValue(u, ix - 1, iy, iz, grid, bc);
                    double yp = GetValue(u, ix, iy + 1, iz, grid, bc);
                    double ym = GetValue(u, ix, iy - 1, iz, grid, bc);
                    double zp = GetValue(u, ix, iy, iz + 1, grid, bc);
                    double zm = GetValue(u, ix, iy, iz - 1, grid, bc);

                    result[idx] = (xp - 2.0 * center + xm) * invDx2
                                + (yp - 2.0 * center + ym) * invDy2
                                + (zp - 2.0 * center + zm) * invDz2;
                }
            }
        }

        return new VectorN(result);
    }

    // ──────────────────────────────────────────────
    //  Gradient: (∂u/∂x, ∂u/∂y, ∂u/∂z)
    // ──────────────────────────────────────────────

    /// <summary>
    /// 3D gradient using central differences.
    /// </summary>
    public static (VectorN dux, VectorN duy, VectorN duz) Gradient3D(VectorN u, Grid3D grid,
        BoundaryCondition bc = BoundaryCondition.Dirichlet)
    {
        int nx = grid.Nx, ny = grid.Ny, nz = grid.Nz;
        var dux = new double[grid.Length];
        var duy = new double[grid.Length];
        var duz = new double[grid.Length];
        double inv2dx = 1.0 / (2.0 * grid.Dx);
        double inv2dy = 1.0 / (2.0 * grid.Dy);
        double inv2dz = 1.0 / (2.0 * grid.Dz);

        for (int iz = 0; iz < nz; iz++)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                for (int ix = 0; ix < nx; ix++)
                {
                    int idx = grid.Index(ix, iy, iz);

                    dux[idx] = (GetValue(u, ix + 1, iy, iz, grid, bc)
                              - GetValue(u, ix - 1, iy, iz, grid, bc)) * inv2dx;

                    duy[idx] = (GetValue(u, ix, iy + 1, iz, grid, bc)
                              - GetValue(u, ix, iy - 1, iz, grid, bc)) * inv2dy;

                    duz[idx] = (GetValue(u, ix, iy, iz + 1, grid, bc)
                              - GetValue(u, ix, iy, iz - 1, grid, bc)) * inv2dz;
                }
            }
        }

        return (new VectorN(dux), new VectorN(duy), new VectorN(duz));
    }

    // ──────────────────────────────────────────────
    //  Divergence: ∂Fx/∂x + ∂Fy/∂y + ∂Fz/∂z
    // ──────────────────────────────────────────────

    /// <summary>
    /// 3D divergence of a vector field (Fx, Fy, Fz) using central differences.
    /// </summary>
    public static VectorN Divergence3D(VectorN fx, VectorN fy, VectorN fz, Grid3D grid,
        BoundaryCondition bc = BoundaryCondition.Dirichlet)
    {
        int nx = grid.Nx, ny = grid.Ny, nz = grid.Nz;
        var result = new double[grid.Length];
        double inv2dx = 1.0 / (2.0 * grid.Dx);
        double inv2dy = 1.0 / (2.0 * grid.Dy);
        double inv2dz = 1.0 / (2.0 * grid.Dz);

        for (int iz = 0; iz < nz; iz++)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                for (int ix = 0; ix < nx; ix++)
                {
                    int idx = grid.Index(ix, iy, iz);

                    double dfxdx = (GetValue(fx, ix + 1, iy, iz, grid, bc)
                                  - GetValue(fx, ix - 1, iy, iz, grid, bc)) * inv2dx;

                    double dfydy = (GetValue(fy, ix, iy + 1, iz, grid, bc)
                                  - GetValue(fy, ix, iy - 1, iz, grid, bc)) * inv2dy;

                    double dfzdz = (GetValue(fz, ix, iy, iz + 1, grid, bc)
                                  - GetValue(fz, ix, iy, iz - 1, grid, bc)) * inv2dz;

                    result[idx] = dfxdx + dfydy + dfzdz;
                }
            }
        }

        return new VectorN(result);
    }

    // ──────────────────────────────────────────────
    //  Advection: v · ∇u  (first-order upwind)
    // ──────────────────────────────────────────────

    /// <summary>
    /// Advection term v·∇u using first-order upwind differencing.
    /// </summary>
    public static VectorN Advection3D(VectorN u, VectorN vx, VectorN vy, VectorN vz, Grid3D grid,
        BoundaryCondition bc = BoundaryCondition.Dirichlet)
    {
        int nx = grid.Nx, ny = grid.Ny, nz = grid.Nz;
        var result = new double[grid.Length];
        double invDx = 1.0 / grid.Dx;
        double invDy = 1.0 / grid.Dy;
        double invDz = 1.0 / grid.Dz;

        for (int iz = 0; iz < nz; iz++)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                for (int ix = 0; ix < nx; ix++)
                {
                    int idx = grid.Index(ix, iy, iz);
                    double center = u[idx];

                    // Upwind in x
                    double dudx = vx[idx] >= 0
                        ? (center - GetValue(u, ix - 1, iy, iz, grid, bc)) * invDx
                        : (GetValue(u, ix + 1, iy, iz, grid, bc) - center) * invDx;

                    // Upwind in y
                    double dudy = vy[idx] >= 0
                        ? (center - GetValue(u, ix, iy - 1, iz, grid, bc)) * invDy
                        : (GetValue(u, ix, iy + 1, iz, grid, bc) - center) * invDy;

                    // Upwind in z
                    double dudz = vz[idx] >= 0
                        ? (center - GetValue(u, ix, iy, iz - 1, grid, bc)) * invDz
                        : (GetValue(u, ix, iy, iz + 1, grid, bc) - center) * invDz;

                    result[idx] = vx[idx] * dudx + vy[idx] * dudy + vz[idx] * dudz;
                }
            }
        }

        return new VectorN(result);
    }

    // ──────────────────────────────────────────────
    //  Poisson solver: ∇²u = f (Gauss-Seidel)
    // ──────────────────────────────────────────────

    /// <summary>
    /// Iterative Gauss-Seidel solver for the 3D Poisson equation ∇²u = f
    /// with Dirichlet boundary conditions supplied via mask and values.
    /// </summary>
    public static (VectorN solution, int iterations) SolvePoisson3D(
        VectorN rhs,
        Grid3D grid,
        bool[] boundaryMask,
        double[] boundaryValues,
        VectorN? initialGuess = null,
        double tolerance = 1e-6,
        int maxIterations = 10000)
    {
        int nx = grid.Nx, ny = grid.Ny, nz = grid.Nz, len = grid.Length;
        double dx2 = grid.Dx * grid.Dx;
        double dy2 = grid.Dy * grid.Dy;
        double dz2 = grid.Dz * grid.Dz;
        double factor = 2.0 / dx2 + 2.0 / dy2 + 2.0 / dz2;

        double[] u = new double[len];
        if (initialGuess.HasValue)
            for (int i = 0; i < len; i++) u[i] = initialGuess.Value[i];

        for (int i = 0; i < len; i++)
            if (boundaryMask[i]) u[i] = boundaryValues[i];

        int iter;
        for (iter = 0; iter < maxIterations; iter++)
        {
            double maxChange = 0;

            for (int iz = 0; iz < nz; iz++)
            {
                for (int iy = 0; iy < ny; iy++)
                {
                    for (int ix = 0; ix < nx; ix++)
                    {
                        int idx = grid.Index(ix, iy, iz);
                        if (boundaryMask[idx]) continue;

                        double xm = ix > 0 ? u[grid.Index(ix - 1, iy, iz)] : 0;
                        double xp = ix < nx - 1 ? u[grid.Index(ix + 1, iy, iz)] : 0;
                        double ym = iy > 0 ? u[grid.Index(ix, iy - 1, iz)] : 0;
                        double yp = iy < ny - 1 ? u[grid.Index(ix, iy + 1, iz)] : 0;
                        double zm = iz > 0 ? u[grid.Index(ix, iy, iz - 1)] : 0;
                        double zp = iz < nz - 1 ? u[grid.Index(ix, iy, iz + 1)] : 0;

                        double newVal = ((xm + xp) / dx2
                                       + (ym + yp) / dy2
                                       + (zm + zp) / dz2
                                       - rhs[idx]) / factor;

                        double change = Math.Abs(newVal - u[idx]);
                        if (change > maxChange) maxChange = change;
                        u[idx] = newVal;
                    }
                }
            }

            if (maxChange < tolerance) { iter++; break; }
        }

        return (new VectorN(u), iter);
    }

    // ──────────────────────────────────────────────
    //  Helper: boundary-aware value lookup
    // ──────────────────────────────────────────────

    private static double GetValue(VectorN u, int ix, int iy, int iz, Grid3D grid, BoundaryCondition bc)
    {
        int nx = grid.Nx, ny = grid.Ny, nz = grid.Nz;

        bool outside = ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz;

        if (!outside)
            return u[grid.Index(ix, iy, iz)];

        switch (bc)
        {
            case BoundaryCondition.Dirichlet:
                return 0.0;

            case BoundaryCondition.Neumann:
                int cx = Math.Max(0, Math.Min(nx - 1, ix));
                int cy = Math.Max(0, Math.Min(ny - 1, iy));
                int cz = Math.Max(0, Math.Min(nz - 1, iz));
                return u[grid.Index(cx, cy, cz)];

            case BoundaryCondition.Periodic:
                int wx = ((ix % nx) + nx) % nx;
                int wy = ((iy % ny) + ny) % ny;
                int wz = ((iz % nz) + nz) % nz;
                return u[grid.Index(wx, wy, wz)];

            default:
                return 0.0;
        }
    }
}
