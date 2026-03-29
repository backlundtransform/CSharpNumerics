using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Numerics.FiniteDifference;
using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Engines.Multiphysics.Solvers;

/// <summary>
/// Solves 2D electrostatics via the Poisson equation ∇²φ = −ρ/ε₀ε_r.
/// Uses the iterative Gauss-Seidel Poisson solver from <see cref="GridOperators"/>.
/// After solving for potential φ, computes E = −∇φ via <see cref="GridOperators.Gradient2D"/>.
/// </summary>
internal class ElectricFieldSolver : IMultiphysicsSolver
{
    private const double Epsilon0 = 8.854187817e-12; // F/m

    public SimulationResult Solve(SimulationBuilder cfg)
    {
        var mat = cfg.Material.Value;
        double epsR = mat.ElectricPermittivity > 0 ? mat.ElectricPermittivity : 1.0;
        double eps = Epsilon0 * epsR;

        double dx = cfg.GeomWidth / cfg.Nx;
        double dy = cfg.GeomHeight / cfg.Ny;
        var grid = new Grid2D(cfg.Nx, cfg.Ny, dx, dy);

        // Right-hand side: −ρ/ε at source locations
        var rhs = new double[grid.Length];
        foreach (var (ix, iy, chargeDensity) in cfg.Sources2D)
        {
            if (ix >= 0 && ix < cfg.Nx && iy >= 0 && iy < cfg.Ny)
                rhs[grid.Index(ix, iy)] = -chargeDensity / eps;
        }

        // Boundary mask & values: edges are Dirichlet
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

        // Solve Poisson
        var (phi, iterations) = GridOperators.SolvePoisson2D(
            new VectorN(rhs), grid, mask, bcValues,
            tolerance: cfg.Tolerance, maxIterations: cfg.MaxIterations);

        // E = −∇φ
        var (dPhiDx, dPhiDy) = GridOperators.Gradient2D(phi, grid, BoundaryCondition.Neumann);
        var exFlat = new double[grid.Length];
        var eyFlat = new double[grid.Length];
        for (int i = 0; i < grid.Length; i++)
        {
            exFlat[i] = -dPhiDx[i];
            eyFlat[i] = -dPhiDy[i];
        }

        double[,] field = grid.ToArray(phi);
        double[,] exArr = grid.ToArray(new VectorN(exFlat));
        double[,] eyArr = grid.ToArray(new VectorN(eyFlat));

        (double min, double max) = FieldMinMax(field, cfg.Nx, cfg.Ny);

        return new SimulationResult
        {
            Type = MultiphysicsType.ElectricField,
            Field = field,
            Ex = exArr,
            Ey = eyArr,
            MaxValue = max,
            MinValue = min,
            Iterations = iterations,
            FinalTime = 0
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
