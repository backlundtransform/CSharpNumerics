using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Physics.SolidMechanics.Enums;
using CSharpNumerics.Physics.SolidMechanics.Interfaces;
using System;
using System.Linq;

namespace CSharpNumerics.Engines.Multiphysics.Solvers;

/// <summary>
/// Solves a 1D Euler–Bernoulli beam problem using analytical deflection
/// formulas provided by an <see cref="IBeamModel"/> from Physics.
/// Supports cantilever and simply supported beams with point loads
/// and/or uniformly distributed loads.
/// </summary>
internal class BeamStressSolver : IMultiphysicsSolver
{
    private readonly IBeamModel _beam;

    internal BeamStressSolver(IBeamModel beam)
    {
        _beam = beam;
    }

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
        for (int i = 0; i < n; i++)
            stress[i] = _beam.BendingStress(moment[i], halfHeight, I);

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

    private void SolveCantilever(
        SimulationBuilder cfg, double EI, double L, int n, double[] positions,
        double[] deflection, double[] moment, double[] shear)
    {
        foreach (var (pos, load) in cfg.Sources1D)
        {
            for (int i = 0; i < n; i++)
            {
                var (d, m, s) = _beam.CantileverPointLoad(load, pos, L, EI, positions[i]);
                deflection[i] += d;
                moment[i] += m;
                shear[i] += s;
            }
        }

        if (cfg.UniformLoad != 0)
        {
            for (int i = 0; i < n; i++)
            {
                var (d, m, s) = _beam.CantileverUniformLoad(cfg.UniformLoad, L, EI, positions[i]);
                deflection[i] += d;
                moment[i] += m;
                shear[i] += s;
            }
        }
    }

    private void SolveSimplySupported(
        SimulationBuilder cfg, double EI, double L, int n, double[] positions,
        double[] deflection, double[] moment, double[] shear)
    {
        foreach (var (pos, load) in cfg.Sources1D)
        {
            for (int i = 0; i < n; i++)
            {
                var (d, m, s) = _beam.SimplySupportedPointLoad(load, pos, L, EI, positions[i]);
                deflection[i] += d;
                moment[i] += m;
                shear[i] += s;
            }
        }

        if (cfg.UniformLoad != 0)
        {
            for (int i = 0; i < n; i++)
            {
                var (d, m, s) = _beam.SimplySupportedUniformLoad(cfg.UniformLoad, L, EI, positions[i]);
                deflection[i] += d;
                moment[i] += m;
                shear[i] += s;
            }
        }
    }

    private void SolveFixedFixed(
        SimulationBuilder cfg, double EI, double L, int n, double[] positions,
        double[] deflection, double[] moment, double[] shear)
    {
        if (cfg.UniformLoad != 0)
        {
            for (int i = 0; i < n; i++)
            {
                var (d, m, s) = _beam.FixedFixedUniformLoad(cfg.UniformLoad, L, EI, positions[i]);
                deflection[i] += d;
                moment[i] += m;
                shear[i] += s;
            }
        }

        foreach (var (pos, load) in cfg.Sources1D)
        {
            for (int i = 0; i < n; i++)
            {
                var (d, m, s) = _beam.FixedFixedPointLoad(load, pos, L, EI, positions[i]);
                deflection[i] += d;
                moment[i] += m;
                shear[i] += s;
            }
        }
    }
}
