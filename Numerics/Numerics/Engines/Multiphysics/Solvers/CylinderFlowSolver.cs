using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Numerics.FiniteDifference;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Multiphysics.Solvers;

/// <summary>
/// Solves 2D incompressible Navier–Stokes around a circular cylinder using
/// Chorin's projection (fractional-step) method on a uniform <see cref="Grid2D"/>.
///
/// Algorithm per time step:
///   1. Predict: v* = vⁿ + Δt·(−(v·∇)v + ν∇²v)
///   2. Pressure Poisson: ∇²p = (∇·v*)/Δt
///   3. Correct: vⁿ⁺¹ = v* − Δt·∇p
///   4. Enforce BCs: inlet, outlet, walls, cylinder no-slip.
/// </summary>
internal class CylinderFlowSolver : IMultiphysicsSolver
{
    public SimulationResult Solve(SimulationBuilder cfg)
    {
        var mat = cfg.Material.Value;
        double nu = mat.KinematicViscosity;
        double rho = mat.Density;
        double uInlet = cfg.InletVelocity;

        double dx = cfg.GeomWidth / cfg.Nx;
        double dy = cfg.GeomHeight / cfg.Ny;
        var grid = new Grid2D(cfg.Nx, cfg.Ny, dx, dy);
        int len = grid.Length;

        // ── Cylinder mask ────────────────────────────────────────
        double cx = cfg.CylinderCenterX;
        double cy = cfg.CylinderCenterY;
        double cr = cfg.CylinderRadius;
        var mask = new bool[len];
        var mask2D = new bool[cfg.Nx, cfg.Ny];
        for (int iy = 0; iy < cfg.Ny; iy++)
        {
            for (int ix = 0; ix < cfg.Nx; ix++)
            {
                double x = ix * dx;
                double y = iy * dy;
                double dist = Math.Sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
                bool inside = dist <= cr;
                int idx = grid.Index(ix, iy);
                mask[idx] = inside;
                mask2D[ix, iy] = inside;
            }
        }

        // ── Pressure Poisson boundary mask ───────────────────────
        // Fixed pressure at outlet (right edge) + inside cylinder
        var pMask = new bool[len];
        var pBC = new double[len];
        for (int iy = 0; iy < cfg.Ny; iy++)
        {
            int outIdx = grid.Index(cfg.Nx - 1, iy);
            pMask[outIdx] = true;
            pBC[outIdx] = 0.0;
        }
        for (int i = 0; i < len; i++)
        {
            if (mask[i])
            {
                pMask[i] = true;
                pBC[i] = 0.0;
            }
        }

        // ── Initial velocity: uniform inlet ──────────────────────
        var vxArr = new double[len];
        var vyArr = new double[len];
        for (int i = 0; i < len; i++)
        {
            if (!mask[i])
                vxArr[i] = uInlet;
        }

        var vx = new VectorN(vxArr);
        var vy = new VectorN(vyArr);
        var p = new VectorN(new double[len]);

        // ── Lift history for Strouhal detection ──────────────────
        var liftHistory = new List<double>();

        // ── Time-step loop ───────────────────────────────────────
        var timeline = new List<double[,]>();
        double time = 0;

        for (int step = 0; step < cfg.Steps; step++)
        {
            double dt = cfg.Dt;

            // CFL clamp: dt ≤ min(dx/|u_max|, dx²/(4ν))
            double uMax = MaxAbs(vx, vy, len);
            if (uMax > 0)
            {
                double dtAdvection = Math.Min(dx, dy) / uMax;
                double dtDiffusion = Math.Min(dx * dx, dy * dy) / (4.0 * nu);
                double dtCFL = 0.4 * Math.Min(dtAdvection, dtDiffusion);
                if (dt > dtCFL) dt = dtCFL;
            }

            // ── Step 1: Advection + Diffusion → v* ──────────────
            var advX = GridOperators.Advection2D(vx, vx, vy, grid, BoundaryCondition.Dirichlet);
            var advY = GridOperators.Advection2D(vy, vx, vy, grid, BoundaryCondition.Dirichlet);
            var lapX = GridOperators.Laplacian2D(vx, grid, BoundaryCondition.Dirichlet);
            var lapY = GridOperators.Laplacian2D(vy, grid, BoundaryCondition.Dirichlet);

            var vxStar = new double[len];
            var vyStar = new double[len];
            for (int i = 0; i < len; i++)
            {
                vxStar[i] = vx[i] + dt * (-advX[i] + nu * lapX[i]);
                vyStar[i] = vy[i] + dt * (-advY[i] + nu * lapY[i]);
            }

            var vxS = new VectorN(vxStar);
            var vyS = new VectorN(vyStar);

            // Apply BCs to intermediate velocity
            ApplyVelocityBC(vxS, vyS, grid, cfg, mask, uInlet);

            // ── Step 2: Pressure Poisson ─────────────────────────
            var div = GridOperators.Divergence2D(vxS, vyS, grid, BoundaryCondition.Dirichlet);
            var rhs = new double[len];
            for (int i = 0; i < len; i++)
                rhs[i] = div[i] / dt;

            var (pSol, _) = GridOperators.SolvePoisson2D(
                new VectorN(rhs), grid, pMask, pBC,
                initialGuess: p,
                tolerance: cfg.Tolerance,
                maxIterations: cfg.MaxIterations);
            p = pSol;

            // ── Step 3: Velocity correction ──────────────────────
            var (dpx, dpy) = GridOperators.Gradient2D(p, grid, BoundaryCondition.Dirichlet);

            var vxNew = new double[len];
            var vyNew = new double[len];
            for (int i = 0; i < len; i++)
            {
                vxNew[i] = vxS[i] - dt * dpx[i];
                vyNew[i] = vyS[i] - dt * dpy[i];
            }

            vx = new VectorN(vxNew);
            vy = new VectorN(vyNew);

            // ── Step 4: Enforce BCs ──────────────────────────────
            ApplyVelocityBC(vx, vy, grid, cfg, mask, uInlet);

            time += dt;

            // Record lift for Strouhal estimation (sum vy on cylinder boundary)
            liftHistory.Add(ComputeLiftProxy(vy, mask, grid, cfg.Nx, cfg.Ny));

            // Store vorticity snapshot every 10 steps (or last step)
            if (step % 10 == 0 || step == cfg.Steps - 1)
                timeline.Add(ComputeVorticityField(vx, vy, grid, cfg.Nx, cfg.Ny));
        }

        // ── Post-processing ─────────────────────────────────────
        var vorticityField = ComputeVorticityField(vx, vy, grid, cfg.Nx, cfg.Ny);

        double cd = ComputeDragCoefficient(p, vx, vy, mask, grid, cfg, nu, rho, uInlet);
        double cl = ComputeLiftCoefficient(p, vx, vy, mask, grid, cfg, nu, rho, uInlet);
        double st = EstimateStrouhal(liftHistory, cfg.Dt, cfg.CylinderRadius * 2.0, uInlet);

        double[,] vxField = grid.ToArray(vx);
        double[,] vyField = grid.ToArray(vy);
        double[,] pField = grid.ToArray(p);
        (double minV, double maxV) = FieldMinMax(vorticityField, cfg.Nx, cfg.Ny);

        return new SimulationResult
        {
            Type = MultiphysicsType.CylinderFlow,
            Field = vorticityField,
            Timeline = timeline,
            Vx = vxField,
            Vy = vyField,
            Pressure = pField,
            Vorticity = vorticityField,
            CylinderMask = mask2D,
            DragCoefficient = cd,
            LiftCoefficient = cl,
            StrouhalNumber = st,
            MaxValue = maxV,
            MinValue = minV,
            Iterations = cfg.Steps,
            FinalTime = time
        };
    }

    // ═══════════════════════════════════════════════════════════════
    //  Boundary conditions
    // ═══════════════════════════════════════════════════════════════

    private static void ApplyVelocityBC(
        VectorN vx, VectorN vy, Grid2D grid,
        SimulationBuilder cfg, bool[] mask, double uInlet)
    {
        int nx = cfg.Nx, ny = cfg.Ny;

        for (int iy = 0; iy < ny; iy++)
        {
            // Inlet (left): Dirichlet u = U∞, v = 0
            int left = grid.Index(0, iy);
            vx[left] = uInlet;
            vy[left] = 0;

            // Outlet (right): Neumann — copy from interior
            int right = grid.Index(nx - 1, iy);
            int inner = grid.Index(nx - 2, iy);
            vx[right] = vx[inner];
            vy[right] = vy[inner];
        }

        for (int ix = 0; ix < nx; ix++)
        {
            // Bottom wall: free-slip (vy = 0, vx copied)
            int bot = grid.Index(ix, 0);
            vy[bot] = 0;
            if (ny > 1) vx[bot] = vx[grid.Index(ix, 1)];

            // Top wall: free-slip (vy = 0, vx copied)
            int top = grid.Index(ix, ny - 1);
            vy[top] = 0;
            if (ny > 1) vx[top] = vx[grid.Index(ix, ny - 2)];
        }

        // Cylinder: no-slip
        for (int i = 0; i < grid.Length; i++)
        {
            if (mask[i])
            {
                vx[i] = 0;
                vy[i] = 0;
            }
        }
    }

    // ═══════════════════════════════════════════════════════════════
    //  Vorticity: ω = ∂vy/∂x − ∂vx/∂y
    // ═══════════════════════════════════════════════════════════════

    private static double[,] ComputeVorticityField(
        VectorN vx, VectorN vy, Grid2D grid, int nx, int ny)
    {
        var (dvydx, _) = GridOperators.Gradient2D(vy, grid, BoundaryCondition.Dirichlet);
        var (_, dvxdy) = GridOperators.Gradient2D(vx, grid, BoundaryCondition.Dirichlet);

        var omega = new double[nx, ny];
        for (int iy = 0; iy < ny; iy++)
            for (int ix = 0; ix < nx; ix++)
            {
                int idx = grid.Index(ix, iy);
                omega[ix, iy] = dvydx[idx] - dvxdy[idx];
            }
        return omega;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Force coefficients (pressure-based approximation)
    // ═══════════════════════════════════════════════════════════════

    private static double ComputeDragCoefficient(
        VectorN p, VectorN vx, VectorN vy, bool[] mask,
        Grid2D grid, SimulationBuilder cfg, double nu, double rho, double uInlet)
    {
        double fx = 0;
        double dy = grid.Dy, dx = grid.Dx;
        int nx = cfg.Nx, ny = cfg.Ny;
        double cx = cfg.CylinderCenterX, cy = cfg.CylinderCenterY;

        // Sum pressure forces on cylinder boundary cells
        for (int iy = 1; iy < ny - 1; iy++)
        {
            for (int ix = 1; ix < nx - 1; ix++)
            {
                if (!mask[grid.Index(ix, iy)]) continue;

                // Check if this cell has at least one fluid neighbour (boundary cell)
                bool hasFluidNeighbour =
                    !mask[grid.Index(ix - 1, iy)] || !mask[grid.Index(ix + 1, iy)] ||
                    !mask[grid.Index(ix, iy - 1)] || !mask[grid.Index(ix, iy + 1)];

                if (!hasFluidNeighbour) continue;

                double x = ix * dx;
                double y = iy * dy;
                double dist = Math.Sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
                if (dist < 1e-15) continue;

                // Outward normal from cylinder centre
                double nxDir = (x - cx) / dist;

                // Pressure contribution to drag
                double pLocal = AverageFluidNeighbourPressure(p, mask, grid, ix, iy, nx, ny);
                fx += pLocal * nxDir * dy;
            }
        }

        double diameter = cfg.CylinderRadius * 2.0;
        double qInf = 0.5 * rho * uInlet * uInlet * diameter;
        return qInf > 0 ? fx / qInf : 0;
    }

    private static double ComputeLiftCoefficient(
        VectorN p, VectorN vx, VectorN vy, bool[] mask,
        Grid2D grid, SimulationBuilder cfg, double nu, double rho, double uInlet)
    {
        double fy = 0;
        double dy = grid.Dy, dx = grid.Dx;
        int nx = cfg.Nx, ny = cfg.Ny;
        double cx = cfg.CylinderCenterX, cy = cfg.CylinderCenterY;

        for (int iy = 1; iy < ny - 1; iy++)
        {
            for (int ix = 1; ix < nx - 1; ix++)
            {
                if (!mask[grid.Index(ix, iy)]) continue;

                bool hasFluidNeighbour =
                    !mask[grid.Index(ix - 1, iy)] || !mask[grid.Index(ix + 1, iy)] ||
                    !mask[grid.Index(ix, iy - 1)] || !mask[grid.Index(ix, iy + 1)];

                if (!hasFluidNeighbour) continue;

                double x = ix * dx;
                double y = iy * dy;
                double dist = Math.Sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
                if (dist < 1e-15) continue;

                double nyDir = (y - cy) / dist;
                double pLocal = AverageFluidNeighbourPressure(p, mask, grid, ix, iy, nx, ny);
                fy += pLocal * nyDir * dx;
            }
        }

        double diameter = cfg.CylinderRadius * 2.0;
        double qInf = 0.5 * rho * uInlet * uInlet * diameter;
        return qInf > 0 ? fy / qInf : 0;
    }

    private static double AverageFluidNeighbourPressure(
        VectorN p, bool[] mask, Grid2D grid, int ix, int iy, int nx, int ny)
    {
        double sum = 0;
        int count = 0;
        int[] dxs = { -1, 1, 0, 0 };
        int[] dys = { 0, 0, -1, 1 };
        for (int d = 0; d < 4; d++)
        {
            int jx = ix + dxs[d], jy = iy + dys[d];
            if (jx >= 0 && jx < nx && jy >= 0 && jy < ny)
            {
                int jIdx = grid.Index(jx, jy);
                if (!mask[jIdx])
                {
                    sum += p[jIdx];
                    count++;
                }
            }
        }
        return count > 0 ? sum / count : 0;
    }

    /// <summary>
    /// Proxy for instantaneous lift: sum of vy just outside cylinder boundary.
    /// </summary>
    private static double ComputeLiftProxy(
        VectorN vy, bool[] mask, Grid2D grid, int nx, int ny)
    {
        double sum = 0;
        for (int iy = 1; iy < ny - 1; iy++)
        {
            for (int ix = 1; ix < nx - 1; ix++)
            {
                int idx = grid.Index(ix, iy);
                if (mask[idx]) continue;

                // Is this fluid cell adjacent to cylinder?
                bool adj = (ix > 0 && mask[grid.Index(ix - 1, iy)]) ||
                           (ix < nx - 1 && mask[grid.Index(ix + 1, iy)]) ||
                           (iy > 0 && mask[grid.Index(ix, iy - 1)]) ||
                           (iy < ny - 1 && mask[grid.Index(ix, iy + 1)]);
                if (adj) sum += vy[idx];
            }
        }
        return sum;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Strouhal estimation from lift history
    // ═══════════════════════════════════════════════════════════════

    private static double EstimateStrouhal(
        List<double> liftHistory, double dt, double diameter, double uInlet)
    {
        if (liftHistory.Count < 20 || diameter <= 0 || uInlet <= 0)
            return 0;

        // Count zero-crossings in the latter half (after transient)
        int start = liftHistory.Count / 2;
        int crossings = 0;
        for (int i = start + 1; i < liftHistory.Count; i++)
        {
            if ((liftHistory[i] >= 0 && liftHistory[i - 1] < 0) ||
                (liftHistory[i] < 0 && liftHistory[i - 1] >= 0))
                crossings++;
        }

        if (crossings < 2) return 0;

        // Each full cycle has 2 zero-crossings
        double cycles = crossings / 2.0;
        double duration = (liftHistory.Count - start) * dt;
        double frequency = cycles / duration;
        return frequency * diameter / uInlet;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Helpers
    // ═══════════════════════════════════════════════════════════════

    private static double MaxAbs(VectorN vx, VectorN vy, int len)
    {
        double max = 0;
        for (int i = 0; i < len; i++)
        {
            double speed = Math.Abs(vx[i]) + Math.Abs(vy[i]);
            if (speed > max) max = speed;
        }
        return max;
    }

    private static (double min, double max) FieldMinMax(double[,] field, int nx, int ny)
    {
        double min = double.MaxValue, max = double.MinValue;
        for (int iy = 0; iy < ny; iy++)
            for (int ix = 0; ix < nx; ix++)
            {
                double v = field[ix, iy];
                if (v < min) min = v;
                if (v > max) max = v;
            }
        return (min, max);
    }
}
