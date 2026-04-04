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
/// on a uniform Grid2D using forward Euler time stepping
/// with Dirichlet boundary conditions.
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

        // Apply Dirichlet BCs to state vector
        void ApplyBC(VectorN state)
        {
            for (int ix = 0; ix < grid.Nx; ix++)
            {
                state[grid.Index(ix, 0)] = cfg.BottomBC;
                state[grid.Index(ix, grid.Ny - 1)] = cfg.TopBC;
            }
            for (int iy = 0; iy < grid.Ny; iy++)
            {
                state[grid.Index(0, iy)] = cfg.LeftBC;
                state[grid.Index(grid.Nx - 1, iy)] = cfg.RightBC;
            }
        }

        ApplyBC(T);

        // Time-step loop (forward Euler)
        var timeline = new List<double[,]> { grid.ToArray(T) };
        double time = 0;

        for (int step = 0; step < cfg.Steps; step++)
        {
            var lap = GridOperators.Laplacian2D(T, grid, BoundaryCondition.Dirichlet);

            var next = new double[grid.Length];
            for (int i = 0; i < grid.Length; i++)
                next[i] = T[i] + cfg.Dt * (alpha * lap[i] + source[i]);

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
