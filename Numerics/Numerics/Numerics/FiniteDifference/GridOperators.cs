using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Numerics.FiniteDifference
{
    /// <summary>
    /// Discrete finite-difference operators that act on a flat <see cref="VectorN"/>
    /// representing values on a <see cref="Grid2D"/>. These produce the right-hand side
    /// of the semi-discrete ODE system obtained via Method of Lines, ready to be fed to
    /// the existing ODE solvers (RungeKutta, EulerMethod, etc.).
    /// </summary>
    public static class GridOperators
    {
        // ──────────────────────────────────────────────
        //  1-D operators (uniform spacing dx, N points)
        // ──────────────────────────────────────────────

        /// <summary>
        /// First derivative du/dx using central differences: (u[i+1] − u[i−1]) / (2·dx).
        /// </summary>
        public static VectorN Gradient1D(VectorN u, double dx,
            BoundaryCondition bc = BoundaryCondition.Dirichlet)
        {
            int n = u.Length;
            var result = new double[n];
            double inv2dx = 1.0 / (2.0 * dx);

            for (int i = 1; i < n - 1; i++)
                result[i] = (u[i + 1] - u[i - 1]) * inv2dx;

            // Boundaries
            switch (bc)
            {
                case BoundaryCondition.Dirichlet:
                    result[0] = (u[1] - 0.0) * inv2dx;
                    result[n - 1] = (0.0 - u[n - 2]) * inv2dx;
                    break;
                case BoundaryCondition.Neumann:
                    result[0] = (u[1] - u[0]) * inv2dx;
                    result[n - 1] = (u[n - 1] - u[n - 2]) * inv2dx;
                    break;
                case BoundaryCondition.Periodic:
                    result[0] = (u[1] - u[n - 1]) * inv2dx;
                    result[n - 1] = (u[0] - u[n - 2]) * inv2dx;
                    break;
            }

            return new VectorN(result);
        }

        /// <summary>
        /// Second derivative d²u/dx² using the 3-point stencil: (u[i+1] − 2u[i] + u[i−1]) / dx².
        /// </summary>
        public static VectorN Laplacian1D(VectorN u, double dx,
            BoundaryCondition bc = BoundaryCondition.Dirichlet)
        {
            int n = u.Length;
            var result = new double[n];
            double invDx2 = 1.0 / (dx * dx);

            for (int i = 1; i < n - 1; i++)
                result[i] = (u[i + 1] - 2.0 * u[i] + u[i - 1]) * invDx2;

            double uLeft, uRight;
            switch (bc)
            {
                case BoundaryCondition.Dirichlet:
                    uLeft = 0.0;
                    uRight = 0.0;
                    break;
                case BoundaryCondition.Neumann:
                    uLeft = u[0];       // mirror: u[-1] = u[0]
                    uRight = u[n - 1];  // mirror: u[N] = u[N-1]
                    break;
                case BoundaryCondition.Periodic:
                    uLeft = u[n - 1];
                    uRight = u[0];
                    break;
                default:
                    uLeft = 0.0;
                    uRight = 0.0;
                    break;
            }

            result[0] = (u[1] - 2.0 * u[0] + uLeft) * invDx2;
            result[n - 1] = (uRight - 2.0 * u[n - 1] + u[n - 2]) * invDx2;

            return new VectorN(result);
        }

        /// <summary>
        /// Fourth derivative d⁴u/dx⁴ using the 5-point central stencil:
        /// (u[i−2] − 4u[i−1] + 6u[i] − 4u[i+1] + u[i+2]) / dx⁴.
        /// Used for Euler-Bernoulli beam equation EIu⁗ = q.
        /// </summary>
        public static VectorN Biharmonic1D(VectorN u, double dx,
            BoundaryCondition bc = BoundaryCondition.Dirichlet)
        {
            int n = u.Length;
            var result = new double[n];
            double invDx4 = 1.0 / (dx * dx * dx * dx);

            for (int i = 0; i < n; i++)
            {
                double um2 = GetValue1D(u, i - 2, n, bc);
                double um1 = GetValue1D(u, i - 1, n, bc);
                double u0 = u[i];
                double up1 = GetValue1D(u, i + 1, n, bc);
                double up2 = GetValue1D(u, i + 2, n, bc);

                result[i] = (um2 - 4.0 * um1 + 6.0 * u0 - 4.0 * up1 + up2) * invDx4;
            }

            return new VectorN(result);
        }

        private static double GetValue1D(VectorN u, int i, int n, BoundaryCondition bc)
        {
            if (i >= 0 && i < n) return u[i];

            switch (bc)
            {
                case BoundaryCondition.Dirichlet:
                    return 0.0;
                case BoundaryCondition.Neumann:
                    return u[Math.Max(0, Math.Min(n - 1, i))];
                case BoundaryCondition.Periodic:
                    return u[((i % n) + n) % n];
                default:
                    return 0.0;
            }
        }

        // ──────────────────────────────────────────────
        //  2-D operators (Grid2D, row-major VectorN)
        // ──────────────────────────────────────────────

        /// <summary>
        /// 2D gradient (∂u/∂x, ∂u/∂y) using central differences.
        /// Returns two <see cref="VectorN"/>: one for each partial derivative.
        /// </summary>
        public static (VectorN dux, VectorN duy) Gradient2D(VectorN u, Grid2D grid,
            BoundaryCondition bc = BoundaryCondition.Dirichlet)
        {
            int nx = grid.Nx, ny = grid.Ny;
            var dux = new double[grid.Length];
            var duy = new double[grid.Length];
            double inv2dx = 1.0 / (2.0 * grid.Dx);
            double inv2dy = 1.0 / (2.0 * grid.Dy);

            for (int iy = 0; iy < ny; iy++)
            {
                for (int ix = 0; ix < nx; ix++)
                {
                    int idx = grid.Index(ix, iy);

                    // ∂u/∂x
                    double uRight = GetValue(u, ix + 1, iy, grid, bc);
                    double uLeft = GetValue(u, ix - 1, iy, grid, bc);
                    dux[idx] = (uRight - uLeft) * inv2dx;

                    // ∂u/∂y
                    double uUp = GetValue(u, ix, iy + 1, grid, bc);
                    double uDown = GetValue(u, ix, iy - 1, grid, bc);
                    duy[idx] = (uUp - uDown) * inv2dy;
                }
            }

            return (new VectorN(dux), new VectorN(duy));
        }

        /// <summary>
        /// 2D Laplacian ∇²u = ∂²u/∂x² + ∂²u/∂y² using the 5-point stencil.
        /// </summary>
        public static VectorN Laplacian2D(VectorN u, Grid2D grid,
            BoundaryCondition bc = BoundaryCondition.Dirichlet)
        {
            int nx = grid.Nx, ny = grid.Ny;
            var result = new double[grid.Length];
            double invDx2 = 1.0 / (grid.Dx * grid.Dx);
            double invDy2 = 1.0 / (grid.Dy * grid.Dy);

            for (int iy = 0; iy < ny; iy++)
            {
                for (int ix = 0; ix < nx; ix++)
                {
                    int idx = grid.Index(ix, iy);
                    double center = u[idx];

                    double uRight = GetValue(u, ix + 1, iy, grid, bc);
                    double uLeft = GetValue(u, ix - 1, iy, grid, bc);
                    double uUp = GetValue(u, ix, iy + 1, grid, bc);
                    double uDown = GetValue(u, ix, iy - 1, grid, bc);

                    result[idx] = (uRight - 2.0 * center + uLeft) * invDx2
                                + (uUp - 2.0 * center + uDown) * invDy2;
                }
            }

            return new VectorN(result);
        }

        /// <summary>
        /// 2D divergence of a vector field (Fx, Fy): ∂Fx/∂x + ∂Fy/∂y.
        /// Both components are stored as flat <see cref="VectorN"/> on the same <see cref="Grid2D"/>.
        /// </summary>
        public static VectorN Divergence2D(VectorN fx, VectorN fy, Grid2D grid,
            BoundaryCondition bc = BoundaryCondition.Dirichlet)
        {
            int nx = grid.Nx, ny = grid.Ny;
            var result = new double[grid.Length];
            double inv2dx = 1.0 / (2.0 * grid.Dx);
            double inv2dy = 1.0 / (2.0 * grid.Dy);

            for (int iy = 0; iy < ny; iy++)
            {
                for (int ix = 0; ix < nx; ix++)
                {
                    int idx = grid.Index(ix, iy);

                    double fxRight = GetValue(fx, ix + 1, iy, grid, bc);
                    double fxLeft = GetValue(fx, ix - 1, iy, grid, bc);
                    double fyUp = GetValue(fy, ix, iy + 1, grid, bc);
                    double fyDown = GetValue(fy, ix, iy - 1, grid, bc);

                    result[idx] = (fxRight - fxLeft) * inv2dx
                                + (fyUp - fyDown) * inv2dy;
                }
            }

            return new VectorN(result);
        }

        // ──────────────────────────────────────────────
        //  Advection: v · ∇u  (upwind scheme)
        // ──────────────────────────────────────────────

        /// <summary>
        /// Advection term v·∇u using first-order upwind differencing for stability.
        /// <paramref name="vx"/> and <paramref name="vy"/> are velocity components on the grid.
        /// </summary>
        public static VectorN Advection2D(VectorN u, VectorN vx, VectorN vy, Grid2D grid,
            BoundaryCondition bc = BoundaryCondition.Dirichlet)
        {
            int nx = grid.Nx, ny = grid.Ny;
            var result = new double[grid.Length];
            double invDx = 1.0 / grid.Dx;
            double invDy = 1.0 / grid.Dy;

            for (int iy = 0; iy < ny; iy++)
            {
                for (int ix = 0; ix < nx; ix++)
                {
                    int idx = grid.Index(ix, iy);
                    double center = u[idx];

                    // Upwind in x
                    double dudx;
                    if (vx[idx] >= 0)
                        dudx = (center - GetValue(u, ix - 1, iy, grid, bc)) * invDx;
                    else
                        dudx = (GetValue(u, ix + 1, iy, grid, bc) - center) * invDx;

                    // Upwind in y
                    double dudy;
                    if (vy[idx] >= 0)
                        dudy = (center - GetValue(u, ix, iy - 1, grid, bc)) * invDy;
                    else
                        dudy = (GetValue(u, ix, iy + 1, grid, bc) - center) * invDy;

                    result[idx] = vx[idx] * dudx + vy[idx] * dudy;
                }
            }

            return new VectorN(result);
        }

        // ──────────────────────────────────────────────
        //  Helper: boundary-aware value lookup
        // ──────────────────────────────────────────────

        private static double GetValue(VectorN u, int ix, int iy, Grid2D grid, BoundaryCondition bc)
        {
            int nx = grid.Nx, ny = grid.Ny;

            bool outsideX = ix < 0 || ix >= nx;
            bool outsideY = iy < 0 || iy >= ny;

            if (!outsideX && !outsideY)
                return u[grid.Index(ix, iy)];

            switch (bc)
            {
                case BoundaryCondition.Dirichlet:
                    return 0.0;

                case BoundaryCondition.Neumann:
                    int clampedX = Math.Max(0, Math.Min(nx - 1, ix));
                    int clampedY = Math.Max(0, Math.Min(ny - 1, iy));
                    return u[grid.Index(clampedX, clampedY)];

                case BoundaryCondition.Periodic:
                    int wrappedX = ((ix % nx) + nx) % nx;
                    int wrappedY = ((iy % ny) + ny) % ny;
                    return u[grid.Index(wrappedX, wrappedY)];

                default:
                    return 0.0;
            }
        }

        // ──────────────────────────────────────────────
        //  Poisson solver: ∇²u = f (Gauss-Seidel)
        // ──────────────────────────────────────────────

        /// <summary>
        /// Iterative Gauss-Seidel solver for the 2D Poisson equation ∇²u = f
        /// on a <see cref="Grid2D"/> with Dirichlet boundary conditions.
        /// Boundary values are supplied via <paramref name="boundaryMask"/>
        /// and <paramref name="boundaryValues"/>; interior nodes are iterated.
        /// </summary>
        /// <param name="rhs">Right-hand side f at each grid point (flat VectorN).</param>
        /// <param name="grid">The computational grid.</param>
        /// <param name="boundaryMask">True for nodes that are fixed (Dirichlet).</param>
        /// <param name="boundaryValues">Fixed values at Dirichlet nodes.</param>
        /// <param name="initialGuess">Starting solution (or null for zeros).</param>
        /// <param name="tolerance">Convergence tolerance on max absolute residual change.</param>
        /// <param name="maxIterations">Maximum iteration count.</param>
        /// <returns>
        /// The converged solution u and the number of iterations used.
        /// </returns>
        public static (VectorN solution, int iterations) SolvePoisson2D(
            VectorN rhs,
            Grid2D grid,
            bool[] boundaryMask,
            double[] boundaryValues,
            VectorN? initialGuess = null,
            double tolerance = 1e-6,
            int maxIterations = 10000)
        {
            int nx = grid.Nx, ny = grid.Ny, len = grid.Length;
            double dx2 = grid.Dx * grid.Dx;
            double dy2 = grid.Dy * grid.Dy;
            double factor = 2.0 / dx2 + 2.0 / dy2;

            // Working solution — mutable copy
            double[] u = new double[len];
            if (initialGuess.HasValue)
                for (int i = 0; i < len; i++) u[i] = initialGuess.Value[i];

            // Apply boundary values
            for (int i = 0; i < len; i++)
                if (boundaryMask[i]) u[i] = boundaryValues[i];

            int iter;
            for (iter = 0; iter < maxIterations; iter++)
            {
                double maxChange = 0;

                for (int iy = 0; iy < ny; iy++)
                {
                    for (int ix = 0; ix < nx; ix++)
                    {
                        int idx = grid.Index(ix, iy);
                        if (boundaryMask[idx]) continue;

                        double uLeft = ix > 0 ? u[grid.Index(ix - 1, iy)] : boundaryValues[idx];
                        double uRight = ix < nx - 1 ? u[grid.Index(ix + 1, iy)] : boundaryValues[idx];
                        double uDown = iy > 0 ? u[grid.Index(ix, iy - 1)] : boundaryValues[idx];
                        double uUp = iy < ny - 1 ? u[grid.Index(ix, iy + 1)] : boundaryValues[idx];

                        double newVal = ((uLeft + uRight) / dx2 + (uDown + uUp) / dy2 - rhs[idx]) / factor;

                        double change = Math.Abs(newVal - u[idx]);
                        if (change > maxChange) maxChange = change;
                        u[idx] = newVal;
                    }
                }

                if (maxChange < tolerance) { iter++; break; }
            }

            return (new VectorN(u), iter);
        }
    }
}
