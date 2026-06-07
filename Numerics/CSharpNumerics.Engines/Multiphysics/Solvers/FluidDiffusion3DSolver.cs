using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Numerics.FiniteDifference;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Multiphysics.Solvers;

/// <summary>
/// Solves the 3D advection-diffusion equation ∂c/∂t = D∇²c − v·∇c + S
/// on a uniform <see cref="Grid3D"/> using forward Euler time stepping
/// with Dirichlet boundary conditions on all six faces.
/// The velocity field v = (vx, vy, vz) is prescribed (not solved).
/// </summary>
internal class FluidDiffusion3DSolver : IMultiphysicsSolver
{
    public SimulationResult Solve(SimulationBuilder cfg)
    {
        double D = cfg.DiffusionCoefficient;
        double dx = cfg.GeomWidth / cfg.Nx;
        double dy = cfg.GeomHeight / cfg.Ny;
        double dz = cfg.GeomDepth / cfg.Nz;
        var grid = new Grid3D(cfg.Nx, cfg.Ny, cfg.Nz, dx, dy, dz);
        int len = grid.Length;

        // Build prescribed velocity field on grid
        VectorN vx, vy, vz;
        if (cfg.HasVelocityField3D)
        {
            var vxArr = new double[len];
            var vyArr = new double[len];
            var vzArr = new double[len];
            for (int iz = 0; iz < cfg.Nz; iz++)
                for (int iy = 0; iy < cfg.Ny; iy++)
                    for (int ix = 0; ix < cfg.Nx; ix++)
                    {
                        double x = ix * dx, y = iy * dy, z = iz * dz;
                        var (vxVal, vyVal, vzVal) = cfg.VelocityField3D(x, y, z);
                        int idx = grid.Index(ix, iy, iz);
                        vxArr[idx] = vxVal;
                        vyArr[idx] = vyVal;
                        vzArr[idx] = vzVal;
                    }
            vx = new VectorN(vxArr);
            vy = new VectorN(vyArr);
            vz = new VectorN(vzArr);
        }
        else
        {
            vx = grid.Zeros();
            vy = grid.Zeros();
            vz = grid.Zeros();
        }

        // CFL check
        bool cflClamped = false;
        double maxSpeed = MaxAbs(vx, vy, vz, len);
        double minD = Math.Min(dx, Math.Min(dy, dz));
        if (maxSpeed > 0)
        {
            double cfl = maxSpeed * cfg.Dt / minD;
            if (cfl > 1.0) cflClamped = true;
        }

        // Initial condition
        VectorN c;
        if (cfg.HasInitialFunc3D)
            c = grid.Initialize(cfg.InitialFunc3D);
        else
            c = grid.Initialize((x, y, z) => cfg.InitialValue);

        // Source term
        var source = grid.Zeros();
        foreach (var (six, siy, siz, value) in cfg.Sources3D)
        {
            if (six >= 0 && six < cfg.Nx && siy >= 0 && siy < cfg.Ny && siz >= 0 && siz < cfg.Nz)
                source[grid.Index(six, siy, siz)] = value;
        }

        // Apply BCs
        void ApplyBC(VectorN state)
        {
            int nx = grid.Nx, ny = grid.Ny, nz = grid.Nz;
            for (int iz2 = 0; iz2 < nz; iz2++)
            {
                for (int ix2 = 0; ix2 < nx; ix2++)
                {
                    state[grid.Index(ix2, 0, iz2)] = cfg.BottomBC;
                    state[grid.Index(ix2, ny - 1, iz2)] = cfg.TopBC;
                }
                for (int iy2 = 0; iy2 < ny; iy2++)
                {
                    state[grid.Index(0, iy2, iz2)] = cfg.LeftBC;
                    state[grid.Index(nx - 1, iy2, iz2)] = cfg.RightBC;
                }
            }
            for (int iy2 = 0; iy2 < ny; iy2++)
                for (int ix2 = 0; ix2 < nx; ix2++)
                {
                    state[grid.Index(ix2, iy2, 0)] = cfg.FrontBC;
                    state[grid.Index(ix2, iy2, nz - 1)] = cfg.BackBC;
                }
        }

        ApplyBC(c);

        // Time stepping
        var timeline = new List<double[,,]> { grid.ToArray(c) };
        double time = 0;

        for (int step = 0; step < cfg.Steps; step++)
        {
            var lap = GridOperators3D.Laplacian3D(c, grid, BoundaryCondition.Dirichlet);
            var adv = GridOperators3D.Advection3D(c, vx, vy, vz, grid, BoundaryCondition.Dirichlet);

            var next = new double[len];
            for (int i = 0; i < len; i++)
                next[i] = c[i] + cfg.Dt * (D * lap[i] - adv[i] + source[i]);

            c = new VectorN(next);
            ApplyBC(c);
            time += cfg.Dt;
            timeline.Add(grid.ToArray(c));
        }

        double[,,] finalField = grid.ToArray(c);
        (double min, double max) = FieldMinMax3D(finalField, cfg.Nx, cfg.Ny, cfg.Nz);

        return new SimulationResult
        {
            Type = MultiphysicsType.FluidDiffusion3D,
            Field3D = finalField,
            Timeline3D = timeline,
            Nx3D = cfg.Nx,
            Ny3D = cfg.Ny,
            Nz3D = cfg.Nz,
            CflClamped = cflClamped,
            MaxValue = max,
            MinValue = min,
            Iterations = cfg.Steps,
            FinalTime = time
        };
    }

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

    private static (double min, double max) FieldMinMax3D(double[,,] field, int nx, int ny, int nz)
    {
        double min = double.MaxValue, max = double.MinValue;
        for (int iz = 0; iz < nz; iz++)
            for (int iy = 0; iy < ny; iy++)
                for (int ix = 0; ix < nx; ix++)
                {
                    double v = field[ix, iy, iz];
                    if (v < min) min = v;
                    if (v > max) max = v;
                }
        return (min, max);
    }
}
