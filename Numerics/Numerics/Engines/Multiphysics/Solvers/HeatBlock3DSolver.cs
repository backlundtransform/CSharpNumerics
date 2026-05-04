using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Numerics.FiniteDifference;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Thermodynamics.Interfaces;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Multiphysics.Solvers;

/// <summary>
/// Solves the 3D heat equation ∂T/∂t = α∇²T + Q/(ρ·c_p)
/// on a uniform <see cref="Grid3D"/> using forward Euler time stepping
/// with Dirichlet boundary conditions on all six faces.
/// </summary>
internal class HeatBlock3DSolver : IMultiphysicsSolver
{
    private readonly IHeatTransferModel _heat;

    internal HeatBlock3DSolver(IHeatTransferModel heat)
    {
        _heat = heat;
    }

    public SimulationResult Solve(SimulationBuilder cfg)
    {
        var mat = cfg.Material.Value;
        double alpha = _heat.ThermalDiffusivity(mat.ThermalConductivity, mat.Density, mat.SpecificHeat);
        double dx = cfg.GeomWidth / cfg.Nx;
        double dy = cfg.GeomHeight / cfg.Ny;
        double dz = cfg.GeomDepth / cfg.Nz;
        var grid = new Grid3D(cfg.Nx, cfg.Ny, cfg.Nz, dx, dy, dz);

        // Initial condition
        VectorN T;
        if (cfg.HasInitialFunc3D)
            T = grid.Initialize(cfg.InitialFunc3D);
        else
            T = grid.Initialize((x, y, z) => cfg.InitialValue);

        // Source term Q/(ρ·c_p)
        var source = grid.Zeros();
        foreach (var (ix, iy, iz, value) in cfg.Sources3D)
        {
            if (ix >= 0 && ix < cfg.Nx && iy >= 0 && iy < cfg.Ny && iz >= 0 && iz < cfg.Nz)
                source[grid.Index(ix, iy, iz)] = _heat.HeatSourceRate(value, mat.Density, mat.SpecificHeat);
        }

        // Apply Dirichlet BCs on all 6 faces
        void ApplyBC(VectorN state)
        {
            int nx = grid.Nx, ny = grid.Ny, nz = grid.Nz;

            for (int iz = 0; iz < nz; iz++)
            {
                for (int ix = 0; ix < nx; ix++)
                {
                    // Bottom face (iy=0) and top face (iy=ny-1)
                    state[grid.Index(ix, 0, iz)] = cfg.BottomBC;
                    state[grid.Index(ix, ny - 1, iz)] = cfg.TopBC;
                }
                for (int iy = 0; iy < ny; iy++)
                {
                    // Left face (ix=0) and right face (ix=nx-1)
                    state[grid.Index(0, iy, iz)] = cfg.LeftBC;
                    state[grid.Index(nx - 1, iy, iz)] = cfg.RightBC;
                }
            }
            for (int iy = 0; iy < ny; iy++)
            {
                for (int ix = 0; ix < nx; ix++)
                {
                    // Front face (iz=0) and back face (iz=nz-1)
                    state[grid.Index(ix, iy, 0)] = cfg.FrontBC;
                    state[grid.Index(ix, iy, nz - 1)] = cfg.BackBC;
                }
            }
        }

        ApplyBC(T);

        // Time-step loop (forward Euler)
        var timeline = new List<double[,,]> { grid.ToArray(T) };
        double time = 0;

        for (int step = 0; step < cfg.Steps; step++)
        {
            var lap = GridOperators3D.Laplacian3D(T, grid, BoundaryCondition.Dirichlet);

            var next = new double[grid.Length];
            for (int i = 0; i < grid.Length; i++)
                next[i] = T[i] + cfg.Dt * (alpha * lap[i] + source[i]);

            T = new VectorN(next);
            ApplyBC(T);
            time += cfg.Dt;
            timeline.Add(grid.ToArray(T));
        }

        double[,,] finalField = grid.ToArray(T);
        (double min, double max) = FieldMinMax3D(finalField, cfg.Nx, cfg.Ny, cfg.Nz);

        return new SimulationResult
        {
            Type = MultiphysicsType.HeatBlock3D,
            Field3D = finalField,
            Timeline3D = timeline,
            Nx3D = cfg.Nx,
            Ny3D = cfg.Ny,
            Nz3D = cfg.Nz,
            MaxValue = max,
            MinValue = min,
            Iterations = cfg.Steps,
            FinalTime = time
        };
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
