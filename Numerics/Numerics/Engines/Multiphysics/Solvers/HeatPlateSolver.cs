using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Numerics.FiniteDifference;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Thermodynamics.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Engines.Multiphysics.Solvers;

/// <summary>
/// Solves the 2D heat equation ∂T/∂t = α∇²T + Q/(ρ·c_p)
/// on a uniform Grid2D using forward Euler time stepping.
/// Supports Dirichlet (fixed temperature) and Robin/convection
/// (−k ∂T/∂n = h(T − T∞)) boundary conditions per face.
/// Physics calculations are delegated to an <see cref="IHeatTransferModel"/>.
/// </summary>
internal class HeatPlateSolver : IMultiphysicsSolver
{
    private readonly IHeatTransferModel _heat;

    internal HeatPlateSolver(IHeatTransferModel heat)
    {
        _heat = heat;
    }

    public SimulationResult Solve(SimulationBuilder cfg)
    {
        var mat = cfg.Material.Value;
        double alpha = _heat.ThermalDiffusivity(mat.ThermalConductivity, mat.Density, mat.SpecificHeat);
        double dx = cfg.GeomWidth / cfg.Nx;
        double dy = cfg.GeomHeight / cfg.Ny;
        var grid = new Grid2D(cfg.Nx, cfg.Ny, dx, dy);

        // Resolve per-face BC: if convection was configured, use it; else Dirichlet from WithBoundary values.
        var topBC = cfg.TopFaceBC ?? FaceBoundaryCondition.Dirichlet(cfg.TopBC);
        var bottomBC = cfg.BottomFaceBC ?? FaceBoundaryCondition.Dirichlet(cfg.BottomBC);
        var leftBC = cfg.LeftFaceBC ?? FaceBoundaryCondition.Dirichlet(cfg.LeftBC);
        var rightBC = cfg.RightFaceBC ?? FaceBoundaryCondition.Dirichlet(cfg.RightBC);

        bool anyConvection = topBC.Type == FaceBCType.Convection
                          || bottomBC.Type == FaceBCType.Convection
                          || leftBC.Type == FaceBCType.Convection
                          || rightBC.Type == FaceBCType.Convection;

        // Initial condition
        VectorN T = cfg.HasInitialFunc
            ? grid.Initialize(cfg.InitialFunc)
            : grid.Initialize((x, y) => cfg.InitialValue);

        // Source term Q/(ρ·c_p) — constant in time
        var source = grid.Zeros();
        foreach (var (ix, iy, value) in cfg.Sources2D)
        {
            if (ix >= 0 && ix < cfg.Nx && iy >= 0 && iy < cfg.Ny)
                source[grid.Index(ix, iy)] = _heat.HeatSourceRate(value, mat.Density, mat.SpecificHeat);
        }

        // Apply Dirichlet BCs to state vector (only on Dirichlet faces)
        void ApplyBC(VectorN state)
        {
            for (int ix = 0; ix < grid.Nx; ix++)
            {
                if (bottomBC.Type == FaceBCType.Dirichlet)
                    state[grid.Index(ix, 0)] = bottomBC.Temperature;
                if (topBC.Type == FaceBCType.Dirichlet)
                    state[grid.Index(ix, grid.Ny - 1)] = topBC.Temperature;
            }
            for (int iy = 0; iy < grid.Ny; iy++)
            {
                if (leftBC.Type == FaceBCType.Dirichlet)
                    state[grid.Index(0, iy)] = leftBC.Temperature;
                if (rightBC.Type == FaceBCType.Dirichlet)
                    state[grid.Index(grid.Nx - 1, iy)] = rightBC.Temperature;
            }
        }

        ApplyBC(T);

        // Time-step loop (forward Euler)
        var timeline = new List<double[,]> { grid.ToArray(T) };
        double time = 0;

        // Choose Laplacian BC: use Neumann on convection faces.
        // When mixing, we use Neumann everywhere and then manually set Dirichlet ghost values
        // through the ApplyBC step. For pure Dirichlet, keep existing behaviour.
        var lapBC = anyConvection ? BoundaryCondition.Neumann : BoundaryCondition.Dirichlet;

        for (int step = 0; step < cfg.Steps; step++)
        {
            var lap = GridOperators.Laplacian2D(T, grid, lapBC);

            var next = new double[grid.Length];
            for (int i = 0; i < grid.Length; i++)
                next[i] = T[i] + cfg.Dt * (alpha * lap[i] + source[i]);

            // Add Robin correction at convection boundary cells
            if (anyConvection)
                ApplyConvectionCorrection(next, T, grid, dx, dy, mat.Density, mat.SpecificHeat,
                    topBC, bottomBC, leftBC, rightBC, cfg.Dt);

            T = new VectorN(next);
            ApplyBC(T);
            time += cfg.Dt;
            timeline.Add(grid.ToArray(T));
        }

        double[,] finalField = grid.ToArray(T);
        (double min, double max) = FieldMinMax(finalField, cfg.Nx, cfg.Ny);

        return new SimulationResult
        {
            Type = MultiphysicsType.HeatPlate,
            Field = finalField,
            Timeline = timeline,
            MaxValue = max,
            MinValue = min,
            Iterations = cfg.Steps,
            FinalTime = time
        };
    }

    /// <summary>
    /// Adds the Robin convection correction −h/(ρ·c_p·dx)·(T−T∞)·dt to boundary cells.
    /// </summary>
    private void ApplyConvectionCorrection(
        double[] next, VectorN T, Grid2D grid,
        double dx, double dy, double density, double specificHeat,
        FaceBoundaryCondition topBC, FaceBoundaryCondition bottomBC,
        FaceBoundaryCondition leftBC, FaceBoundaryCondition rightBC,
        double dt)
    {
        int nx = grid.Nx, ny = grid.Ny;

        // Bottom face (iy=0) — normal spacing is dy
        if (bottomBC.Type == FaceBCType.Convection)
        {
            for (int ix = 0; ix < nx; ix++)
            {
                int idx = grid.Index(ix, 0);
                next[idx] += dt * _heat.ConvectiveBoundaryRate(
                    bottomBC.HeatTransferCoefficient, density, specificHeat, dy, T[idx], bottomBC.Temperature);
            }
        }

        // Top face (iy=ny-1) — normal spacing is dy
        if (topBC.Type == FaceBCType.Convection)
        {
            for (int ix = 0; ix < nx; ix++)
            {
                int idx = grid.Index(ix, ny - 1);
                next[idx] += dt * _heat.ConvectiveBoundaryRate(
                    topBC.HeatTransferCoefficient, density, specificHeat, dy, T[idx], topBC.Temperature);
            }
        }

        // Left face (ix=0) — normal spacing is dx
        if (leftBC.Type == FaceBCType.Convection)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                int idx = grid.Index(0, iy);
                next[idx] += dt * _heat.ConvectiveBoundaryRate(
                    leftBC.HeatTransferCoefficient, density, specificHeat, dx, T[idx], leftBC.Temperature);
            }
        }

        // Right face (ix=nx-1) — normal spacing is dx
        if (rightBC.Type == FaceBCType.Convection)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                int idx = grid.Index(nx - 1, iy);
                next[idx] += dt * _heat.ConvectiveBoundaryRate(
                    rightBC.HeatTransferCoefficient, density, specificHeat, dx, T[idx], rightBC.Temperature);
            }
        }
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
