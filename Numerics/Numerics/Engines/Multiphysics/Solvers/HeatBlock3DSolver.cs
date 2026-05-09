using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Numerics.FiniteDifference;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Thermodynamics.Interfaces;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Multiphysics.Solvers;

/// <summary>
/// Solves the 3D heat equation ∂T/∂t = α∇²T + Q/(ρ·c_p)
/// on a uniform <see cref="Grid3D"/> using forward Euler time stepping.
/// Supports Dirichlet (fixed temperature) and Robin/convection
/// (−k ∂T/∂n = h(T − T∞)) boundary conditions on all six faces.
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

        // Resolve per-face BC
        var topBC = cfg.TopFaceBC ?? FaceBoundaryCondition.Dirichlet(cfg.TopBC);
        var bottomBC = cfg.BottomFaceBC ?? FaceBoundaryCondition.Dirichlet(cfg.BottomBC);
        var leftBC = cfg.LeftFaceBC ?? FaceBoundaryCondition.Dirichlet(cfg.LeftBC);
        var rightBC = cfg.RightFaceBC ?? FaceBoundaryCondition.Dirichlet(cfg.RightBC);
        var frontBC = cfg.FrontFaceBC ?? FaceBoundaryCondition.Dirichlet(cfg.FrontBC);
        var backBC = cfg.BackFaceBC ?? FaceBoundaryCondition.Dirichlet(cfg.BackBC);

        bool anyConvection = topBC.Type == FaceBCType.Convection
                          || bottomBC.Type == FaceBCType.Convection
                          || leftBC.Type == FaceBCType.Convection
                          || rightBC.Type == FaceBCType.Convection
                          || frontBC.Type == FaceBCType.Convection
                          || backBC.Type == FaceBCType.Convection;

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

        // Apply Dirichlet BCs on faces that use fixed temperature
        void ApplyBC(VectorN state)
        {
            int nx = grid.Nx, ny = grid.Ny, nz = grid.Nz;

            for (int iz = 0; iz < nz; iz++)
            {
                for (int ix = 0; ix < nx; ix++)
                {
                    if (bottomBC.Type == FaceBCType.Dirichlet)
                        state[grid.Index(ix, 0, iz)] = bottomBC.Temperature;
                    if (topBC.Type == FaceBCType.Dirichlet)
                        state[grid.Index(ix, ny - 1, iz)] = topBC.Temperature;
                }
                for (int iy = 0; iy < ny; iy++)
                {
                    if (leftBC.Type == FaceBCType.Dirichlet)
                        state[grid.Index(0, iy, iz)] = leftBC.Temperature;
                    if (rightBC.Type == FaceBCType.Dirichlet)
                        state[grid.Index(nx - 1, iy, iz)] = rightBC.Temperature;
                }
            }
            for (int iy = 0; iy < ny; iy++)
            {
                for (int ix = 0; ix < nx; ix++)
                {
                    if (frontBC.Type == FaceBCType.Dirichlet)
                        state[grid.Index(ix, iy, 0)] = frontBC.Temperature;
                    if (backBC.Type == FaceBCType.Dirichlet)
                        state[grid.Index(ix, iy, nz - 1)] = backBC.Temperature;
                }
            }
        }

        ApplyBC(T);

        // Time-step loop (forward Euler)
        var timeline = new List<double[,,]> { grid.ToArray(T) };
        double time = 0;

        var lapBC = anyConvection ? BoundaryCondition.Neumann : BoundaryCondition.Dirichlet;

        for (int step = 0; step < cfg.Steps; step++)
        {
            var lap = GridOperators3D.Laplacian3D(T, grid, lapBC);

            var next = new double[grid.Length];
            for (int i = 0; i < grid.Length; i++)
                next[i] = T[i] + cfg.Dt * (alpha * lap[i] + source[i]);

            if (anyConvection)
                ApplyConvectionCorrection(next, T, grid, dx, dy, dz,
                    mat.Density, mat.SpecificHeat,
                    topBC, bottomBC, leftBC, rightBC, frontBC, backBC, cfg.Dt);

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

    /// <summary>
    /// Adds the Robin convection correction −h/(ρ·c_p·dx)·(T−T∞)·dt to boundary cells on all six faces.
    /// </summary>
    private void ApplyConvectionCorrection(
        double[] next, VectorN T, Grid3D grid,
        double dx, double dy, double dz,
        double density, double specificHeat,
        FaceBoundaryCondition topBC, FaceBoundaryCondition bottomBC,
        FaceBoundaryCondition leftBC, FaceBoundaryCondition rightBC,
        FaceBoundaryCondition frontBC, FaceBoundaryCondition backBC,
        double dt)
    {
        int nx = grid.Nx, ny = grid.Ny, nz = grid.Nz;

        for (int iz = 0; iz < nz; iz++)
        {
            for (int ix = 0; ix < nx; ix++)
            {
                // Bottom face (iy=0)
                if (bottomBC.Type == FaceBCType.Convection)
                {
                    int idx = grid.Index(ix, 0, iz);
                    next[idx] += dt * _heat.ConvectiveBoundaryRate(
                        bottomBC.HeatTransferCoefficient, density, specificHeat, dy, T[idx], bottomBC.Temperature);
                }
                // Top face (iy=ny-1)
                if (topBC.Type == FaceBCType.Convection)
                {
                    int idx = grid.Index(ix, ny - 1, iz);
                    next[idx] += dt * _heat.ConvectiveBoundaryRate(
                        topBC.HeatTransferCoefficient, density, specificHeat, dy, T[idx], topBC.Temperature);
                }
            }
            for (int iy = 0; iy < ny; iy++)
            {
                // Left face (ix=0)
                if (leftBC.Type == FaceBCType.Convection)
                {
                    int idx = grid.Index(0, iy, iz);
                    next[idx] += dt * _heat.ConvectiveBoundaryRate(
                        leftBC.HeatTransferCoefficient, density, specificHeat, dx, T[idx], leftBC.Temperature);
                }
                // Right face (ix=nx-1)
                if (rightBC.Type == FaceBCType.Convection)
                {
                    int idx = grid.Index(nx - 1, iy, iz);
                    next[idx] += dt * _heat.ConvectiveBoundaryRate(
                        rightBC.HeatTransferCoefficient, density, specificHeat, dx, T[idx], rightBC.Temperature);
                }
            }
        }
        for (int iy = 0; iy < ny; iy++)
        {
            for (int ix = 0; ix < nx; ix++)
            {
                // Front face (iz=0)
                if (frontBC.Type == FaceBCType.Convection)
                {
                    int idx = grid.Index(ix, iy, 0);
                    next[idx] += dt * _heat.ConvectiveBoundaryRate(
                        frontBC.HeatTransferCoefficient, density, specificHeat, dz, T[idx], frontBC.Temperature);
                }
                // Back face (iz=nz-1)
                if (backBC.Type == FaceBCType.Convection)
                {
                    int idx = grid.Index(ix, iy, nz - 1);
                    next[idx] += dt * _heat.ConvectiveBoundaryRate(
                        backBC.HeatTransferCoefficient, density, specificHeat, dz, T[idx], backBC.Temperature);
                }
            }
        }
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
