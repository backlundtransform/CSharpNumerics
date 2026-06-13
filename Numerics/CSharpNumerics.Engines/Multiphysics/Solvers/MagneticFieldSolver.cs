using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Numerics.FiniteDifference;
using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Engines.Multiphysics.Solvers;

/// <summary>
/// Solves 2D magnetostatics via the Poisson equation for the magnetic vector potential:
///   ∇²A = −μ₀μ_r · J
/// where A is the z-component of the vector potential (2D: A = Az ẑ),
/// J is the current density (z-component), and μ_r is relative permeability.
/// After solving for A, computes B = ∇×A → Bx = ∂A/∂y, By = −∂A/∂x.
/// </summary>
internal class MagneticFieldSolver : IMultiphysicsSolver
{
    private const double Mu0 = 4.0 * Math.PI * 1e-7; // H/m

    public SimulationResult Solve(SimulationBuilder cfg)
    {
        var mat = cfg.Material.Value;
        double muR = mat.MagneticPermeability;
        double mu = Mu0 * muR;

        double dx = cfg.GeomWidth / cfg.Nx;
        double dy = cfg.GeomHeight / cfg.Ny;
        var grid = new Grid2D(cfg.Nx, cfg.Ny, dx, dy);

        // Right-hand side: −μ·J at source locations
        // Sources2D values represent current density J (A/m²)
        var rhs = new double[grid.Length];
        foreach (var (ix, iy, currentDensity) in cfg.Sources2D)
        {
            if (ix >= 0 && ix < cfg.Nx && iy >= 0 && iy < cfg.Ny)
                rhs[grid.Index(ix, iy)] = -mu * currentDensity;
        }

        // Boundary mask & values: edges are Dirichlet (A = 0 at infinity approximation)
        var mask = new bool[grid.Length];
        var bcValues = new double[grid.Length];

        for (int ix = 0; ix < cfg.Nx; ix++)
        {
            int bot = grid.Index(ix, 0);
            int top = grid.Index(ix, cfg.Ny - 1);
            mask[bot] = true; bcValues[bot] = cfg.BottomBC;
            mask[top] = true; bcValues[top] = cfg.TopBC;
        }
        for (int iy = 0; iy < cfg.Ny; iy++)
        {
            int left = grid.Index(0, iy);
            int right = grid.Index(cfg.Nx - 1, iy);
            mask[left] = true; bcValues[left] = cfg.LeftBC;
            mask[right] = true; bcValues[right] = cfg.RightBC;
        }

        // Solve Poisson for vector potential A
        var (A, iterations) = GridOperators.SolvePoisson2D(
            new VectorN(rhs), grid, mask, bcValues,
            tolerance: cfg.Tolerance, maxIterations: cfg.MaxIterations);

        // Compute B = curl(A ẑ) → Bx = ∂A/∂y, By = −∂A/∂x
        var (dAdx, dAdy) = GridOperators.Gradient2D(A, grid, BoundaryCondition.Neumann);
        var bxFlat = new double[grid.Length];
        var byFlat = new double[grid.Length];
        for (int i = 0; i < grid.Length; i++)
        {
            bxFlat[i] = dAdy[i];     // Bx = ∂A/∂y
            byFlat[i] = -dAdx[i];    // By = −∂A/∂x
        }

        double[,] fieldA = grid.ToArray(A);
        double[,] bxArr = grid.ToArray(new VectorN(bxFlat));
        double[,] byArr = grid.ToArray(new VectorN(byFlat));

        // Compute |B| as the primary field for min/max
        double[,] bMag = new double[cfg.Nx, cfg.Ny];
        double min = double.MaxValue, max = double.MinValue;
        for (int iy = 0; iy < cfg.Ny; iy++)
            for (int ix = 0; ix < cfg.Nx; ix++)
            {
                double v = Math.Sqrt(bxArr[ix, iy] * bxArr[ix, iy] + byArr[ix, iy] * byArr[ix, iy]);
                bMag[ix, iy] = v;
                if (v < min) min = v;
                if (v > max) max = v;
            }

        return new SimulationResult
        {
            Type = MultiphysicsType.MagneticField,
            Field = bMag,
            VectorPotential = fieldA,
            Bx = bxArr,
            By = byArr,
            MaxValue = max,
            MinValue = min,
            Iterations = iterations,
            FinalTime = 0
        };
    }
}
