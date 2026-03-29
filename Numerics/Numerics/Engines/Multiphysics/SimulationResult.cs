using CSharpNumerics.Engines.Multiphysics.Enums;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Multiphysics;

/// <summary>
/// Result of a multiphysics simulation.
/// Populated fields depend on the simulation type.
/// </summary>
public class SimulationResult
{
    /// <summary>Which simulation produced this result.</summary>
    public MultiphysicsType Type { get; internal set; }

    // ── 2D field data (HeatPlate, ElectricField) ─────────────────

    /// <summary>Final scalar field values [ix, iy]. Temperature (K) or potential (V).</summary>
    public double[,] Field { get; internal set; }

    /// <summary>Time series of 2D snapshots (transient simulations).</summary>
    public List<double[,]> Timeline { get; internal set; }

    /// <summary>Time series of 1D snapshots (transient PipeFlow).</summary>
    public List<double[]> Timeline1D { get; internal set; }

    // ── 1D array data (BeamStress, PipeFlow) ─────────────────────

    /// <summary>Primary 1D result: deflection (m) for beams, velocity (m/s) for pipe flow.</summary>
    public double[] Values { get; internal set; }

    /// <summary>Node positions along the domain (m).</summary>
    public double[] Positions { get; internal set; }

    // ── Beam-specific ────────────────────────────────────────────

    /// <summary>Bending moment M(x) in N·m.</summary>
    public double[] BendingMoment { get; internal set; }

    /// <summary>Shear force V(x) in N.</summary>
    public double[] ShearForce { get; internal set; }

    /// <summary>Maximum bending stress σ(x) in Pa.</summary>
    public double[] Stress { get; internal set; }

    // ── Electric field–specific ──────────────────────────────────

    /// <summary>Electric field x-component E_x [ix, iy] in V/m.</summary>
    public double[,] Ex { get; internal set; }

    /// <summary>Electric field y-component E_y [ix, iy] in V/m.</summary>
    public double[,] Ey { get; internal set; }

    // ── Metadata ─────────────────────────────────────────────────

    /// <summary>Maximum value in the primary result field.</summary>
    public double MaxValue { get; internal set; }

    /// <summary>Minimum value in the primary result field.</summary>
    public double MinValue { get; internal set; }

    /// <summary>Number of time steps or solver iterations completed.</summary>
    public int Iterations { get; internal set; }

    /// <summary>Final simulation time (transient) or solver residual (steady-state).</summary>
    public double FinalTime { get; internal set; }
}
