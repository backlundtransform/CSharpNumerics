using CSharpNumerics.Physics.SolidMechanics.Enums;
using System;

namespace CSharpNumerics.Engines.Multiphysics.Snapshots;

/// <summary>
/// Immutable snapshot of 1D beam analysis results at a single load case.
/// Stores deflection, bending moment, shear force, and stress arrays.
/// </summary>
public class BeamSnapshot
{
    /// <summary>Node positions along the beam (m).</summary>
    public double[] Positions { get; }

    /// <summary>Deflection u(x) in metres.</summary>
    public double[] Deflection { get; }

    /// <summary>Bending moment M(x) in N·m.</summary>
    public double[] BendingMoment { get; }

    /// <summary>Shear force V(x) in N.</summary>
    public double[] ShearForce { get; }

    /// <summary>Maximum bending stress σ(x) in Pa.</summary>
    public double[] Stress { get; }

    /// <summary>Beam support type used.</summary>
    public BeamSupport Support { get; }

    /// <summary>Beam length in metres.</summary>
    public double Length { get; }

    /// <summary>Number of nodes.</summary>
    public int NodeCount => Positions.Length;

    /// <summary>Maximum absolute deflection.</summary>
    public double MaxDeflection { get; }

    /// <summary>Maximum absolute stress.</summary>
    public double MaxStress { get; }

    public BeamSnapshot(
        double[] positions, double[] deflection,
        double[] bendingMoment, double[] shearForce, double[] stress,
        BeamSupport support, double length)
    {
        Positions = (double[])positions.Clone();
        Deflection = (double[])deflection.Clone();
        BendingMoment = (double[])bendingMoment.Clone();
        ShearForce = (double[])shearForce.Clone();
        Stress = (double[])stress.Clone();
        Support = support;
        Length = length;

        double maxDef = 0, maxStr = 0;
        for (int i = 0; i < deflection.Length; i++)
        {
            double d = Math.Abs(deflection[i]);
            if (d > maxDef) maxDef = d;
            double s = Math.Abs(stress[i]);
            if (s > maxStr) maxStr = s;
        }
        MaxDeflection = maxDef;
        MaxStress = maxStr;
    }

    /// <summary>
    /// Creates a <see cref="BeamSnapshot"/> from a <see cref="SimulationResult"/>
    /// of type <see cref="MultiphysicsType.BeamStress"/>.
    /// </summary>
    public static BeamSnapshot FromResult(SimulationResult result, BeamSupport support, double length)
    {
        return new BeamSnapshot(
            result.Positions, result.Values,
            result.BendingMoment, result.ShearForce, result.Stress,
            support, length);
    }
}
