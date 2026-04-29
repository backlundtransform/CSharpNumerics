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

    // ── Cylinder flow–specific ───────────────────────────────────

    /// <summary>X-velocity field vx [ix, iy] in m/s (CylinderFlow).</summary>
    public double[,] Vx { get; internal set; }

    /// <summary>Y-velocity field vy [ix, iy] in m/s (CylinderFlow).</summary>
    public double[,] Vy { get; internal set; }

    /// <summary>Pressure field p [ix, iy] in Pa (CylinderFlow).</summary>
    public double[,] Pressure { get; internal set; }

    /// <summary>Vorticity field ω = ∂vy/∂x − ∂vx/∂y [ix, iy] in 1/s (CylinderFlow).</summary>
    public double[,] Vorticity { get; internal set; }

    /// <summary>Boolean mask [ix, iy] — true for cells inside the cylinder obstacle.</summary>
    public bool[,] CylinderMask { get; internal set; }

    /// <summary>Drag coefficient Cd (CylinderFlow).</summary>
    public double DragCoefficient { get; internal set; }

    /// <summary>Lift coefficient Cl (CylinderFlow).</summary>
    public double LiftCoefficient { get; internal set; }

    /// <summary>Strouhal number St = fD/U∞ (CylinderFlow, 0 if no shedding detected).</summary>
    public double StrouhalNumber { get; internal set; }

    // ── Magnetic field–specific ──────────────────────────────────

    /// <summary>Magnetic vector potential A [ix, iy] in T·m (MagneticField).</summary>
    public double[,] VectorPotential { get; internal set; }

    /// <summary>Magnetic field x-component Bx [ix, iy] in T (MagneticField).</summary>
    public double[,] Bx { get; internal set; }

    /// <summary>Magnetic field y-component By [ix, iy] in T (MagneticField).</summary>
    public double[,] By { get; internal set; }

    // ── Plane stress–specific ────────────────────────────────────

    /// <summary>X-displacement field ux [ix, iy] in m (PlaneStress).</summary>
    public double[,] Ux { get; internal set; }

    /// <summary>Y-displacement field uy [ix, iy] in m (PlaneStress).</summary>
    public double[,] Uy { get; internal set; }

    /// <summary>Normal stress σxx [ix, iy] in Pa (PlaneStress).</summary>
    public double[,] StressXX { get; internal set; }

    /// <summary>Normal stress σyy [ix, iy] in Pa (PlaneStress).</summary>
    public double[,] StressYY { get; internal set; }

    /// <summary>Shear stress τxy [ix, iy] in Pa (PlaneStress).</summary>
    public double[,] StressXY { get; internal set; }

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
