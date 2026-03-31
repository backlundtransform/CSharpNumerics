using CSharpNumerics.Engines.Multiphysics.Enums;

using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.SolidMechanics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Engines.Multiphysics.Solvers;

/// <summary>
/// Solves a 1D Euler–Bernoulli beam problem using analytical deflection
/// formulas from <see cref="SolidExtensions"/>.
/// Supports cantilever and simply supported beams with point loads
/// and/or uniformly distributed loads.
/// </summary>
internal class BeamStressSolver : IMultiphysicsSolver
{
    public SimulationResult Solve(SimulationBuilder cfg)
    {
        var mat = cfg.Material.Value;
        double E = mat.YoungsModulus;
        double I = cfg.ComputeSecondMoment();
        double EI = E * I;
        double L = cfg.GeomLength;
        int n = cfg.Nodes;
        double dx = L / (n - 1);

        double halfHeight = cfg.SectionHeight > 0
            ? cfg.SectionHeight / 2.0
            : cfg.SectionRadius > 0 ? cfg.SectionRadius : 0;

        // Node positions
        var positions = new double[n];
        for (int i = 0; i < n; i++)
            positions[i] = i * dx;

        // Compute deflection, moment, shear from analytical formulas
        var deflection = new double[n];
        var moment = new double[n];
        var shear = new double[n];
        var stress = new double[n];

        switch (cfg.Support.Value)
        {
            case BeamSupport.Cantilever:
                SolveCantilever(cfg, EI, L, n, positions, deflection, moment, shear);
                break;
            case BeamSupport.SimplySupported:
                SolveSimplySupported(cfg, EI, L, n, positions, deflection, moment, shear);
                break;
            case BeamSupport.FixedFixed:
                SolveFixedFixed(cfg, EI, L, n, positions, deflection, moment, shear);
                break;
        }

        // Bending stress σ = M · y_max / I
        if (halfHeight > 0 && I > 0)
        {
            for (int i = 0; i < n; i++)
                stress[i] = Math.Abs(moment[i]) * halfHeight / I;
        }

        double maxDef = deflection.Max();
        double minDef = deflection.Min();

        return new SimulationResult
        {
            Type = MultiphysicsType.BeamStress,
            Values = deflection,
            Positions = positions,
            BendingMoment = moment,
            ShearForce = shear,
            Stress = stress,
            MaxValue = Math.Max(Math.Abs(maxDef), Math.Abs(minDef)),
            MinValue = Math.Min(maxDef, minDef),
            Iterations = 0,
            FinalTime = 0
        };
    }

    private static void SolveCantilever(
        SimulationBuilder cfg, double EI, double L, int n, double[] positions,
        double[] deflection, double[] moment, double[] shear)
    {
        // Point loads at the free end or at specified positions
        foreach (var (pos, load) in cfg.Sources1D)
        {
            double a = pos; // distance from fixed end to load application point
            for (int i = 0; i < n; i++)
            {
                double x = positions[i];
                if (x <= a)
                {
                    // Deflection: u(x) = P·x²(3a − x) / (6EI) for x ≤ a
                    deflection[i] += load * x * x * (3.0 * a - x) / (6.0 * EI);
                    // Moment: M(x) = P(a − x) for x ≤ a
                    moment[i] += load * (a - x);
                    // Shear: V = P for x < a
                    shear[i] += load;
                }
                else
                {
                    // Deflection: u(x) = P·a²(3x − a) / (6EI) for x > a
                    deflection[i] += load * a * a * (3.0 * x - a) / (6.0 * EI);
                    // Moment = 0 for x > a, Shear = 0 for x > a
                }
            }
        }

        // Uniform distributed load
        if (cfg.UniformLoad != 0)
        {
            double q = cfg.UniformLoad;
            for (int i = 0; i < n; i++)
            {
                double x = positions[i];
                deflection[i] += q.CantileverUniformLoadDeflection(L, EI, x);
                moment[i] += q * (L - x) * (L - x) / 2.0;
                shear[i] += q * (L - x);
            }
        }
    }

    private static void SolveSimplySupported(
        SimulationBuilder cfg, double EI, double L, int n, double[] positions,
        double[] deflection, double[] moment, double[] shear)
    {
        // Point loads (assumed at midspan for analytical formulas)
        foreach (var (pos, load) in cfg.Sources1D)
        {
            double a = pos; // distance from left support
            double b = L - a;
            for (int i = 0; i < n; i++)
            {
                double x = positions[i];
                if (x <= a)
                {
                    deflection[i] += load * b * x * (L * L - b * b - x * x) / (6.0 * L * EI);
                    moment[i] += load * b * x / L;
                    shear[i] += load * b / L;
                }
                else
                {
                    deflection[i] += load * a * (L - x) * (2.0 * L * x - a * a - x * x) / (6.0 * L * EI);
                    moment[i] += load * a * (L - x) / L;
                    shear[i] += -load * a / L;
                }
            }
        }

        // Uniform distributed load
        if (cfg.UniformLoad != 0)
        {
            double q = cfg.UniformLoad;
            for (int i = 0; i < n; i++)
            {
                double x = positions[i];
                deflection[i] += q.SimplySupportedUniformLoadDeflection(L, EI, x);
                moment[i] += q * x * (L - x) / 2.0;
                shear[i] += q * (L / 2.0 - x);
            }
        }
    }

    private static void SolveFixedFixed(
        SimulationBuilder cfg, double EI, double L, int n, double[] positions,
        double[] deflection, double[] moment, double[] shear)
    {
        // Uniform load on fixed-fixed beam: δ(x) = q·x²(L−x)² / (24EI)
        if (cfg.UniformLoad != 0)
        {
            double q = cfg.UniformLoad;
            for (int i = 0; i < n; i++)
            {
                double x = positions[i];
                double Lx = L - x;
                deflection[i] += q * x * x * Lx * Lx / (24.0 * EI);
                moment[i] += q * (L * x - x * x) / 2.0 - q * L * L / 12.0;
                shear[i] += q * (L / 2.0 - x);
            }
        }

        // Point load at position a on fixed-fixed beam
        foreach (var (pos, load) in cfg.Sources1D)
        {
            double a = pos;
            double b = L - a;
            for (int i = 0; i < n; i++)
            {
                double x = positions[i];
                if (x <= a)
                {
                    deflection[i] += load * b * b * x * x * (3.0 * a * L - x * (3.0 * a + b))
                                     / (6.0 * EI * L * L * L);
                }
                else
                {
                    deflection[i] += load * a * a * (L - x) * (L - x) * (3.0 * b * L - (L - x) * (3.0 * b + a))
                                     / (6.0 * EI * L * L * L);
                }
                // Moment and shear for fixed-fixed with point load
                double Ma = -load * a * b * b / (L * L);
                double Mb = -load * a * a * b / (L * L);
                double Ra = load * b * b * (3.0 * a + b) / (L * L * L);
                if (x <= a)
                {
                    shear[i] += Ra;
                    moment[i] += Ma + Ra * x;
                }
                else
                {
                    shear[i] += Ra - load;
                    moment[i] += Ma + Ra * x - load * (x - a);
                }
            }
        }
    }
}
