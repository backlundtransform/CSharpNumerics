using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Numerics.FiniteDifference;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Multiphysics.Solvers;

/// <summary>
/// Solves the 2D incompressible Navier–Stokes equations on a rectangular domain
/// using the Chorin projection method (fractional step):
///   1. Advection + diffusion (explicit forward Euler)
///   2. Pressure Poisson solve (enforce ∇·v = 0)
///   3. Velocity correction: v = v* − (dt/ρ)∇p
/// Produces velocity fields (Vx, Vy), pressure, and timeline snapshots.
/// </summary>
internal class FluidFlow2DSolver : IMultiphysicsSolver
{
    public SimulationResult Solve(SimulationBuilder cfg)
    {
        var mat = cfg.Material.Value;
        double nu = mat.KinematicViscosity;
        double rho = mat.Density;

        double dx = cfg.GeomWidth / cfg.Nx;
        double dy = cfg.GeomHeight / cfg.Ny;
        var grid = new Grid2D(cfg.Nx, cfg.Ny, dx, dy);
        int N = grid.Length;

        // Initialize velocity and pressure fields
        var vx = new double[N];
        var vy = new double[N];
        var p = new double[N];

        // Apply initial conditions
        if (cfg.HasInitialFunc)
        {
            var init = grid.Initialize(cfg.InitialFunc);
            for (int i = 0; i < N; i++) vx[i] = init[i];
        }
        else if (cfg.InitialValue != 0)
        {
            for (int i = 0; i < N; i++) vx[i] = cfg.InitialValue;
        }

        // Inlet BC (left boundary): use InletVelocity if set, else LeftBC
        double inletU = cfg.InletVelocity > 0 ? cfg.InletVelocity : cfg.LeftBC;

        var timeline = new List<double[,]>();
        double time = 0;

        // Store initial snapshot
        timeline.Add(ToVelocityMagnitude(vx, vy, grid));

        for (int step = 0; step < cfg.Steps; step++)
        {
            // CFL-limited time step
            double maxV = 1e-10;
            for (int i = 0; i < N; i++)
            {
                double speed = Math.Sqrt(vx[i] * vx[i] + vy[i] * vy[i]);
                if (speed > maxV) maxV = speed;
            }
            double dtCFL = 0.25 * Math.Min(dx, dy) / maxV;
            double dtDiff = 0.2 * dx * dx / (nu > 0 ? nu : 1e-10);
            double dt = Math.Min(cfg.Dt, Math.Min(dtCFL, dtDiff));

            // ─── Step 1: Advection + Diffusion (Forward Euler) ───
            var vxStar = new double[N];
            var vyStar = new double[N];

            var vxVec = new VectorN(vx);
            var vyVec = new VectorN(vy);

            // Diffusion: ν∇²v
            var lapVx = GridOperators.Laplacian2D(vxVec, grid, BoundaryCondition.Neumann);
            var lapVy = GridOperators.Laplacian2D(vyVec, grid, BoundaryCondition.Neumann);

            // Advection: -(v·∇)v using upwind
            var advVx = GridOperators.Advection2D(vxVec, vxVec, vyVec, grid, BoundaryCondition.Neumann);
            var advVy = GridOperators.Advection2D(vyVec, vxVec, vyVec, grid, BoundaryCondition.Neumann);

            for (int i = 0; i < N; i++)
            {
                vxStar[i] = vx[i] + dt * (nu * lapVx[i] - advVx[i]);
                vyStar[i] = vy[i] + dt * (nu * lapVy[i] - advVy[i]);
            }

            // ─── Apply BCs to intermediate velocity ───
            ApplyVelocityBC(vxStar, vyStar, grid, inletU, cfg.TopBC, cfg.BottomBC);

            // ─── Step 2: Pressure Poisson solve ───
            // ∇²p = (ρ/dt)∇·v*
            var divStar = GridOperators.Divergence2D(
                new VectorN(vxStar), new VectorN(vyStar), grid, BoundaryCondition.Neumann);

            var rhs = new double[N];
            double rhoDivDt = rho / dt;
            for (int i = 0; i < N; i++)
                rhs[i] = rhoDivDt * divStar[i];

            // Boundary mask for pressure (Neumann everywhere, reference at corner)
            var mask = new bool[N];
            var bcValues = new double[N];
            mask[0] = true; // Fix pressure at one point to remove null space
            bcValues[0] = 0.0;

            var (pSol, _) = GridOperators.SolvePoisson2D(
                new VectorN(rhs), grid, mask, bcValues,
                initialGuess: new VectorN(p),
                tolerance: cfg.Tolerance, maxIterations: cfg.MaxIterations);

            for (int i = 0; i < N; i++) p[i] = pSol[i];

            // ─── Step 3: Velocity correction ───
            var (dpx, dpy) = GridOperators.Gradient2D(pSol, grid, BoundaryCondition.Neumann);
            double dtOverRho = dt / rho;
            for (int i = 0; i < N; i++)
            {
                vx[i] = vxStar[i] - dtOverRho * dpx[i];
                vy[i] = vyStar[i] - dtOverRho * dpy[i];
            }

            // Re-apply BCs after correction
            ApplyVelocityBC(vx, vy, grid, inletU, cfg.TopBC, cfg.BottomBC);

            time += dt;
            timeline.Add(ToVelocityMagnitude(vx, vy, grid));
        }

        // Build result
        double[,] vxArr = grid.ToArray(new VectorN(vx));
        double[,] vyArr = grid.ToArray(new VectorN(vy));
        double[,] pArr = grid.ToArray(new VectorN(p));
        double[,] magField = ToVelocityMagnitude(vx, vy, grid);

        (double min, double max) = FieldMinMax(magField, cfg.Nx, cfg.Ny);

        return new SimulationResult
        {
            Type = MultiphysicsType.FluidFlow2D,
            Field = magField,
            Vx = vxArr,
            Vy = vyArr,
            Pressure = pArr,
            Timeline = timeline,
            MaxValue = max,
            MinValue = min,
            Iterations = cfg.Steps,
            FinalTime = time
        };
    }

    private static void ApplyVelocityBC(double[] vx, double[] vy, Grid2D grid, double inletU, double topBC, double bottomBC)
    {
        int nx = grid.Nx, ny = grid.Ny;

        // Left (inlet): Dirichlet
        for (int iy = 0; iy < ny; iy++)
        {
            int idx = grid.Index(0, iy);
            vx[idx] = inletU;
            vy[idx] = 0;
        }

        // Right (outlet): Neumann (copy from interior)
        for (int iy = 0; iy < ny; iy++)
        {
            int idx = grid.Index(nx - 1, iy);
            int inner = grid.Index(nx - 2, iy);
            vx[idx] = vx[inner];
            vy[idx] = vy[inner];
        }

        // Bottom wall: no-slip
        for (int ix = 0; ix < nx; ix++)
        {
            int idx = grid.Index(ix, 0);
            vx[idx] = bottomBC;
            vy[idx] = 0;
        }

        // Top wall: no-slip
        for (int ix = 0; ix < nx; ix++)
        {
            int idx = grid.Index(ix, ny - 1);
            vx[idx] = topBC;
            vy[idx] = 0;
        }
    }

    private static double[,] ToVelocityMagnitude(double[] vx, double[] vy, Grid2D grid)
    {
        var mag = new double[grid.Nx, grid.Ny];
        for (int iy = 0; iy < grid.Ny; iy++)
            for (int ix = 0; ix < grid.Nx; ix++)
            {
                int idx = grid.Index(ix, iy);
                mag[ix, iy] = Math.Sqrt(vx[idx] * vx[idx] + vy[idx] * vy[idx]);
            }
        return mag;
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
