using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Numerics.FiniteDifference;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Engines.Multiphysics.Solvers;

/// <summary>
/// Solves 1D transient Hagen–Poiseuille pipe flow via finite differences.
/// Radial momentum equation: ∂v/∂t = ν·(1/r)·∂/∂r(r·∂v/∂r) − (1/ρ)·dP/dx.
/// Simplified to 1D radial diffusion on r ∈ [0, R] with v(R) = 0 (no-slip)
/// and ∂v/∂r(0) = 0 (symmetry).
/// </summary>
internal class PipeFlowSolver : IMultiphysicsSolver
{
    public SimulationResult Solve(SimulationBuilder cfg)
    {
        var mat = cfg.Material.Value;
        double nu = mat.KinematicViscosity;
        double rho = mat.Density;
        double R = cfg.GeomRadius;
        int n = cfg.Nodes;
        double dr = R / (n - 1);
        double dPdx = cfg.PressureGradient;

        // Radial node positions: r = [0, dr, 2dr, ..., R]
        var positions = new double[n];
        for (int i = 0; i < n; i++)
            positions[i] = i * dr;

        // Initial velocity (zero)
        var v = new double[n];

        // Analytical steady-state for reference: v(r) = -(dP/dx)(R² − r²)/(4μ)
        double mu = mat.DynamicViscosity;
        double drivingTerm = -dPdx / rho; // body force per unit mass

        // Forward Euler time stepping with cylindrical Laplacian
        var timeline = new List<double[]> { (double[])v.Clone() };
        double time = 0;

        for (int step = 0; step < cfg.Steps; step++)
        {
            var vNew = new double[n];

            // Interior nodes: ν * [d²v/dr² + (1/r)(dv/dr)] + driving
            for (int i = 1; i < n - 1; i++)
            {
                double r = positions[i];
                double d2vdr2 = (v[i + 1] - 2.0 * v[i] + v[i - 1]) / (dr * dr);
                double dvdr = (v[i + 1] - v[i - 1]) / (2.0 * dr);
                double cylLaplacian = d2vdr2 + dvdr / r;
                vNew[i] = v[i] + cfg.Dt * (nu * cylLaplacian + drivingTerm);
            }

            // r = 0 (symmetry): use L'Hôpital — (1/r)(dv/dr)|_{r→0} = d²v/dr²
            // so Laplacian = 2·d²v/dr²
            {
                double d2vdr2 = 2.0 * (v[1] - v[0]) / (dr * dr);
                vNew[0] = v[0] + cfg.Dt * (nu * d2vdr2 + drivingTerm);
            }

            // r = R (no-slip wall)
            vNew[n - 1] = 0;

            v = vNew;
            time += cfg.Dt;
            timeline.Add((double[])v.Clone());
        }

        return new SimulationResult
        {
            Type = MultiphysicsType.PipeFlow,
            Values = v,
            Positions = positions,
            Timeline1D = timeline,
            MaxValue = v.Max(),
            MinValue = v.Min(),
            Iterations = cfg.Steps,
            FinalTime = time
        };
    }
}
