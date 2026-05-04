using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Numerics.FiniteDifference;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Multiphysics.Solvers;

/// <summary>
/// Solves 3D incompressible Navier–Stokes around a circular cylinder using
/// Chorin's projection method on a uniform <see cref="Grid3D"/>.
///
/// The cylinder axis is aligned with z. The cylinder mask is constant across all z-slices.
/// Periodic boundary conditions are applied in the z-direction.
///
/// Algorithm per time step:
///   1. Predict: v* = vⁿ + Δt·(−(v·∇)v + ν∇²v)
///   2. Pressure Poisson: ∇²p = (∇·v*)/Δt
///   3. Correct: vⁿ⁺¹ = v* − Δt·∇p
///   4. Enforce BCs: inlet, outlet, top/bottom walls, cylinder no-slip, periodic z.
/// </summary>
internal class CylinderFlow3DSolver : IMultiphysicsSolver
{
    public SimulationResult Solve(SimulationBuilder cfg)
    {
        var mat = cfg.Material.Value;
        double nu = mat.KinematicViscosity;
        double rho = mat.Density;
        double uInlet = cfg.InletVelocity;

        int nx = cfg.Nx, ny = cfg.Ny, nz = cfg.Nz;
        double dx = cfg.GeomWidth / nx;
        double dy = cfg.GeomHeight / ny;
        double dz = cfg.GeomDepth / nz;
        var grid = new Grid3D(nx, ny, nz, dx, dy, dz);
        int len = grid.Length;

        // ── Cylinder mask (constant along z-axis) ────────────────
        double cx = cfg.CylinderCenterX;
        double cy = cfg.CylinderCenterY;
        double cr = cfg.CylinderRadius;
        var mask = new bool[len];
        var mask2D = new bool[nx, ny];
        for (int iy = 0; iy < ny; iy++)
        {
            for (int ix = 0; ix < nx; ix++)
            {
                double x = ix * dx;
                double y = iy * dy;
                double dist = Math.Sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
                bool inside = dist <= cr;
                mask2D[ix, iy] = inside;
                for (int iz = 0; iz < nz; iz++)
                    mask[grid.Index(ix, iy, iz)] = inside;
            }
        }

        // ── Pressure Poisson boundary mask ───────────────────────
        // Outlet face (ix = nx-1) is p=0 reference, plus cylinder cells.
        var pMask = new bool[len];
        var pBC = new double[len];
        for (int iz = 0; iz < nz; iz++)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                int outIdx = grid.Index(nx - 1, iy, iz);
                pMask[outIdx] = true;
                pBC[outIdx] = 0.0;
            }
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
        var vzArr = new double[len];
        for (int i = 0; i < len; i++)
        {
            if (!mask[i])
                vxArr[i] = uInlet;
        }

        var vx = new VectorN(vxArr);
        var vy = new VectorN(vyArr);
        var vz = new VectorN(vzArr);
        var p = new VectorN(new double[len]);

        // ── Lift history for Strouhal detection ──────────────────
        var liftHistory = new List<double>();

        // ── Time-step loop ───────────────────────────────────────
        var timeline = new List<double[,,]>();
        double time = 0;
        bool cflClamped = false;

        for (int step = 0; step < cfg.Steps; step++)
        {
            double dt = cfg.Dt;

            // CFL clamp
            double uMax = MaxAbs(vx, vy, vz, len);
            if (uMax > 0)
            {
                double minD = Math.Min(dx, Math.Min(dy, dz));
                double dtAdv = minD / uMax;
                double dtDiff = minD * minD / (6.0 * nu);
                double dtCFL = 0.4 * Math.Min(dtAdv, dtDiff);
                if (dt > dtCFL) { dt = dtCFL; cflClamped = true; }
            }

            // ── Step 1: Advection + Diffusion → v* ──────────────
            var advX = GridOperators3D.Advection3D(vx, vx, vy, vz, grid, BoundaryCondition.Periodic);
            var advY = GridOperators3D.Advection3D(vy, vx, vy, vz, grid, BoundaryCondition.Periodic);
            var advZ = GridOperators3D.Advection3D(vz, vx, vy, vz, grid, BoundaryCondition.Periodic);
            var lapX = GridOperators3D.Laplacian3D(vx, grid, BoundaryCondition.Periodic);
            var lapY = GridOperators3D.Laplacian3D(vy, grid, BoundaryCondition.Periodic);
            var lapZ = GridOperators3D.Laplacian3D(vz, grid, BoundaryCondition.Periodic);

            var vxStar = new double[len];
            var vyStar = new double[len];
            var vzStar = new double[len];
            for (int i = 0; i < len; i++)
            {
                vxStar[i] = vx[i] + dt * (-advX[i] + nu * lapX[i]);
                vyStar[i] = vy[i] + dt * (-advY[i] + nu * lapY[i]);
                vzStar[i] = vz[i] + dt * (-advZ[i] + nu * lapZ[i]);
            }

            var vxS = new VectorN(vxStar);
            var vyS = new VectorN(vyStar);
            var vzS = new VectorN(vzStar);

            ApplyVelocityBC(vxS, vyS, vzS, grid, cfg, mask, uInlet);

            // ── Step 2: Pressure Poisson ─────────────────────────
            var div = GridOperators3D.Divergence3D(vxS, vyS, vzS, grid, BoundaryCondition.Periodic);
            var rhs = new double[len];
            for (int i = 0; i < len; i++)
                rhs[i] = div[i] / dt;

            var (pSol, _) = GridOperators3D.SolvePoisson3D(
                new VectorN(rhs), grid, pMask, pBC,
                initialGuess: p,
                tolerance: cfg.Tolerance,
                maxIterations: cfg.MaxIterations);
            p = pSol;

            // ── Step 3: Velocity correction ──────────────────────
            var (dpx, dpy, dpz) = GridOperators3D.Gradient3D(p, grid, BoundaryCondition.Periodic);

            var vxNew = new double[len];
            var vyNew = new double[len];
            var vzNew = new double[len];
            for (int i = 0; i < len; i++)
            {
                vxNew[i] = vxS[i] - dt * dpx[i];
                vyNew[i] = vyS[i] - dt * dpy[i];
                vzNew[i] = vzS[i] - dt * dpz[i];
            }

            vx = new VectorN(vxNew);
            vy = new VectorN(vyNew);
            vz = new VectorN(vzNew);

            // ── Step 4: Enforce BCs ──────────────────────────────
            ApplyVelocityBC(vx, vy, vz, grid, cfg, mask, uInlet);

            time += dt;

            liftHistory.Add(ComputeLiftProxy(vy, mask, len));

            // Store timeline snapshot periodically (vorticity z-component at mid-z)
            if (step % 10 == 0 || step == cfg.Steps - 1)
                timeline.Add(grid.ToArray(vx)); // store vx as representative timeline field
        }

        // ── Post-processing ─────────────────────────────────────
        double cd = ComputeDragCoefficient(p, mask, mask2D, grid, cfg, rho, uInlet);
        double cl = ComputeLiftCoefficient(p, mask, mask2D, grid, cfg, rho, uInlet);
        double st = EstimateStrouhal(liftHistory, cfg.Dt, cfg.CylinderRadius * 2.0, uInlet);

        double[,,] vxField = grid.ToArray(vx);
        double[,,] vyField = grid.ToArray(vy);
        double[,,] vzField = grid.ToArray(vz);
        double[,,] pField = grid.ToArray(p);

        // Velocity magnitude as primary Field3D
        double[,,] magField = new double[nx, ny, nz];
        double minVal = double.MaxValue, maxVal = double.MinValue;
        for (int iz = 0; iz < nz; iz++)
            for (int iy = 0; iy < ny; iy++)
                for (int ix = 0; ix < nx; ix++)
                {
                    double speed = Math.Sqrt(
                        vxField[ix, iy, iz] * vxField[ix, iy, iz] +
                        vyField[ix, iy, iz] * vyField[ix, iy, iz] +
                        vzField[ix, iy, iz] * vzField[ix, iy, iz]);
                    magField[ix, iy, iz] = speed;
                    if (speed < minVal) minVal = speed;
                    if (speed > maxVal) maxVal = speed;
                }

        return new SimulationResult
        {
            Type = MultiphysicsType.CylinderFlow3D,
            Field3D = magField,
            Timeline3D = timeline,
            Vx3D = vxField,
            Vy3D = vyField,
            Vz3D = vzField,
            Pressure3D = pField,
            CylinderMask3D = mask2D,
            DragCoefficient = cd,
            LiftCoefficient = cl,
            StrouhalNumber = st,
            CflClamped = cflClamped,
            Nx3D = nx,
            Ny3D = ny,
            Nz3D = nz,
            MaxValue = maxVal,
            MinValue = minVal,
            Iterations = cfg.Steps,
            FinalTime = time
        };
    }

    // ═════════════════════════════════════════════════════════════
    //  Boundary conditions
    // ═════════════════════════════════════════════════════════════

    private static void ApplyVelocityBC(
        VectorN vx, VectorN vy, VectorN vz,
        Grid3D grid, SimulationBuilder cfg, bool[] mask, double uInlet)
    {
        int nx = grid.Nx, ny = grid.Ny, nz = grid.Nz;

        for (int iz = 0; iz < nz; iz++)
        {
            // Inlet (ix=0): uniform flow
            for (int iy = 0; iy < ny; iy++)
            {
                int idx = grid.Index(0, iy, iz);
                vx[idx] = uInlet;
                vy[idx] = 0;
                vz[idx] = 0;
            }

            // Outlet (ix=nx-1): Neumann (copy from interior)
            for (int iy = 0; iy < ny; iy++)
            {
                int idx = grid.Index(nx - 1, iy, iz);
                int inner = grid.Index(nx - 2, iy, iz);
                vx[idx] = vx[inner];
                vy[idx] = vy[inner];
                vz[idx] = vz[inner];
            }

            // Top/bottom walls: free-slip (vy=0, copy vx/vz from interior)
            for (int ix = 0; ix < nx; ix++)
            {
                int bot = grid.Index(ix, 0, iz);
                vy[bot] = 0;
                if (ny > 1)
                {
                    vx[bot] = vx[grid.Index(ix, 1, iz)];
                    vz[bot] = vz[grid.Index(ix, 1, iz)];
                }

                int top = grid.Index(ix, ny - 1, iz);
                vy[top] = 0;
                if (ny > 1)
                {
                    vx[top] = vx[grid.Index(ix, ny - 2, iz)];
                    vz[top] = vz[grid.Index(ix, ny - 2, iz)];
                }
            }
        }

        // Periodic z: copy iz=0 ← iz=nz-2, iz=nz-1 ← iz=1
        // (handled implicitly by the Periodic BC in operators, but also
        //  enforce explicitly for the velocity correction step)
        for (int iy = 0; iy < ny; iy++)
            for (int ix = 0; ix < nx; ix++)
            {
                int front = grid.Index(ix, iy, 0);
                int backInner = grid.Index(ix, iy, nz - 2);
                int back = grid.Index(ix, iy, nz - 1);
                int frontInner = grid.Index(ix, iy, 1);

                vx[front] = vx[backInner];
                vy[front] = vy[backInner];
                vz[front] = vz[backInner];

                vx[back] = vx[frontInner];
                vy[back] = vy[frontInner];
                vz[back] = vz[frontInner];
            }

        // Cylinder no-slip
        for (int i = 0; i < grid.Length; i++)
        {
            if (mask[i])
            {
                vx[i] = 0;
                vy[i] = 0;
                vz[i] = 0;
            }
        }
    }

    // ═════════════════════════════════════════════════════════════
    //  Force coefficients (pressure-based approximation)
    // ═════════════════════════════════════════════════════════════

    private static double ComputeDragCoefficient(
        VectorN p, bool[] mask, bool[,] mask2D,
        Grid3D grid, SimulationBuilder cfg, double rho, double uInlet)
    {
        double fx = 0;
        int nx = grid.Nx, ny = grid.Ny, nz = grid.Nz;
        double dx = grid.Dx, dy = grid.Dy, dz = grid.Dz;
        double cx = cfg.CylinderCenterX, cy = cfg.CylinderCenterY;

        for (int iz = 0; iz < nz; iz++)
        {
            for (int iy = 1; iy < ny - 1; iy++)
            {
                for (int ix = 1; ix < nx - 1; ix++)
                {
                    int idx = grid.Index(ix, iy, iz);
                    if (!mask[idx]) continue;

                    bool hasFluid =
                        !mask[grid.Index(ix - 1, iy, iz)] || !mask[grid.Index(ix + 1, iy, iz)] ||
                        !mask[grid.Index(ix, iy - 1, iz)] || !mask[grid.Index(ix, iy + 1, iz)];

                    if (!hasFluid) continue;

                    double x = ix * dx, y = iy * dy;
                    double dist = Math.Sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
                    if (dist < 1e-15) continue;

                    double nxDir = (x - cx) / dist;
                    double pLocal = AverageFluidPressure(p, mask, grid, ix, iy, iz, nx, ny, nz);
                    fx += pLocal * nxDir * dy * dz;
                }
            }
        }

        double diameter = cfg.CylinderRadius * 2.0;
        double depth = cfg.GeomDepth;
        double qInf = 0.5 * rho * uInlet * uInlet * diameter * depth;
        return qInf > 0 ? fx / qInf : 0;
    }

    private static double ComputeLiftCoefficient(
        VectorN p, bool[] mask, bool[,] mask2D,
        Grid3D grid, SimulationBuilder cfg, double rho, double uInlet)
    {
        double fy = 0;
        int nx = grid.Nx, ny = grid.Ny, nz = grid.Nz;
        double dx = grid.Dx, dy = grid.Dy, dz = grid.Dz;
        double cx = cfg.CylinderCenterX, cy = cfg.CylinderCenterY;

        for (int iz = 0; iz < nz; iz++)
        {
            for (int iy = 1; iy < ny - 1; iy++)
            {
                for (int ix = 1; ix < nx - 1; ix++)
                {
                    int idx = grid.Index(ix, iy, iz);
                    if (!mask[idx]) continue;

                    bool hasFluid =
                        !mask[grid.Index(ix - 1, iy, iz)] || !mask[grid.Index(ix + 1, iy, iz)] ||
                        !mask[grid.Index(ix, iy - 1, iz)] || !mask[grid.Index(ix, iy + 1, iz)];

                    if (!hasFluid) continue;

                    double x = ix * dx, y = iy * dy;
                    double dist = Math.Sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
                    if (dist < 1e-15) continue;

                    double nyDir = (y - cy) / dist;
                    double pLocal = AverageFluidPressure(p, mask, grid, ix, iy, iz, nx, ny, nz);
                    fy += pLocal * nyDir * dx * dz;
                }
            }
        }

        double diameter = cfg.CylinderRadius * 2.0;
        double depth = cfg.GeomDepth;
        double qInf = 0.5 * rho * uInlet * uInlet * diameter * depth;
        return qInf > 0 ? fy / qInf : 0;
    }

    private static double AverageFluidPressure(
        VectorN p, bool[] mask, Grid3D grid, int ix, int iy, int iz, int nx, int ny, int nz)
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
                int jIdx = grid.Index(jx, jy, iz);
                if (!mask[jIdx])
                {
                    sum += p[jIdx];
                    count++;
                }
            }
        }
        return count > 0 ? sum / count : 0;
    }

    private static double ComputeLiftProxy(VectorN vy, bool[] mask, int len)
    {
        double sum = 0;
        for (int i = 0; i < len; i++)
        {
            if (!mask[i]) continue;
            sum += vy[i];
        }
        return sum;
    }

    // ═════════════════════════════════════════════════════════════
    //  Strouhal estimation
    // ═════════════════════════════════════════════════════════════

    private static double EstimateStrouhal(
        List<double> liftHistory, double dt, double diameter, double uInlet)
    {
        if (liftHistory.Count < 20 || diameter <= 0 || uInlet <= 0)
            return 0;

        int start = liftHistory.Count / 2;
        int crossings = 0;
        for (int i = start + 1; i < liftHistory.Count; i++)
        {
            if ((liftHistory[i] >= 0 && liftHistory[i - 1] < 0) ||
                (liftHistory[i] < 0 && liftHistory[i - 1] >= 0))
                crossings++;
        }

        if (crossings < 2) return 0;

        double cycles = crossings / 2.0;
        double duration = (liftHistory.Count - start) * dt;
        double frequency = cycles / duration;
        return frequency * diameter / uInlet;
    }

    // ═════════════════════════════════════════════════════════════
    //  Helpers
    // ═════════════════════════════════════════════════════════════

    private static double MaxAbs(VectorN vx, VectorN vy, VectorN vz, int len)
    {
        double max = 0;
        for (int i = 0; i < len; i++)
        {
            double speed = Math.Abs(vx[i]) + Math.Abs(vy[i]) + Math.Abs(vz[i]);
            if (speed > max) max = speed;
        }
        return max;
    }
}
