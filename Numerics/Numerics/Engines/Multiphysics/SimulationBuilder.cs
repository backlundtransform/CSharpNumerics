using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Engines.Multiphysics.Solvers;
using CSharpNumerics.Physics.Electromagnetism;
using CSharpNumerics.Physics.FluidDynamics;
using CSharpNumerics.Physics.Materials.Engineering;
using CSharpNumerics.Physics.SolidMechanics;
using CSharpNumerics.Physics.SolidMechanics.Enums;
using CSharpNumerics.Physics.Thermodynamics;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Multiphysics;

/// <summary>
/// Fluent builder that collects all configuration for a multiphysics simulation
/// and dispatches to the appropriate solver on <see cref="Solve"/> / <see cref="Run"/>.
/// </summary>
public class SimulationBuilder
{
    // ── Type ─────────────────────────────────────────────────────
    internal MultiphysicsType Type { get; }

    // ── Material ─────────────────────────────────────────────────
    internal EngineeringMaterial? Material { get; private set; }

    // ── 2D Geometry (HeatPlate / ElectricField) ──────────────────
    internal double GeomWidth { get; private set; }
    internal double GeomHeight { get; private set; }
    internal int Nx { get; private set; }
    internal int Ny { get; private set; }

    // ── 1D Geometry (PipeFlow / BeamStress) ──────────────────────
    internal double GeomLength { get; private set; }
    internal double GeomRadius { get; private set; }
    internal int Nodes { get; private set; }

    // ── Cross-section (BeamStress) ───────────────────────────────
    internal double SectionWidth { get; private set; }
    internal double SectionHeight { get; private set; }
    internal double SectionRadius { get; private set; }
    internal double SecondMomentOverride { get; private set; }

    // ── Cylinder (CylinderFlow) ──────────────────────────────────
    internal double CylinderCenterX { get; private set; }
    internal double CylinderCenterY { get; private set; }
    internal double CylinderRadius { get; private set; }
    internal double InletVelocity { get; private set; }

    // ── Boundary conditions ──────────────────────────────────────
    internal double TopBC { get; private set; }
    internal double BottomBC { get; private set; }
    internal double LeftBC { get; private set; }
    internal double RightBC { get; private set; }
    internal BeamSupport? Support { get; private set; }
    internal double PressureGradient { get; private set; }

    // ── Initial condition ────────────────────────────────────────
    internal double InitialValue { get; private set; }
    internal Func<double, double, double> InitialFunc { get; private set; }
    internal bool HasInitialFunc { get; private set; }

    // ── Sources ──────────────────────────────────────────────────
    internal List<(int ix, int iy, double value)> Sources2D { get; } = new();
    internal List<(double position, double value)> Sources1D { get; } = new();
    internal double UniformLoad { get; private set; }

    // ── Solver parameters ────────────────────────────────────────
    internal double Dt { get; private set; } = 0.01;
    internal int Steps { get; private set; }
    internal int MaxIterations { get; private set; } = 10000;
    internal double Tolerance { get; private set; } = 1e-6;

    internal SimulationBuilder(MultiphysicsType type)
    {
        Type = type;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Material
    // ═══════════════════════════════════════════════════════════════

    /// <summary>Set the engineering material (drives α, ν, ε, E, etc.).</summary>
    public SimulationBuilder WithMaterial(EngineeringMaterial material)
    {
        Material = material;
        return this;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Geometry
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Define a 2D rectangular domain (HeatPlate, ElectricField).
    /// </summary>
    /// <param name="width">Physical width in metres.</param>
    /// <param name="height">Physical height in metres.</param>
    /// <param name="nx">Grid cells in x.</param>
    /// <param name="ny">Grid cells in y.</param>
    public SimulationBuilder WithGeometry(double width, double height, int nx, int ny)
    {
        GeomWidth = width;
        GeomHeight = height;
        Nx = nx;
        Ny = ny;
        return this;
    }

    /// <summary>
    /// Define a 1D domain (BeamStress, PipeFlow without radius).
    /// </summary>
    /// <param name="length">Domain length in metres.</param>
    /// <param name="nodes">Number of nodes.</param>
    public SimulationBuilder WithGeometry(double length, int nodes)
    {
        GeomLength = length;
        Nodes = nodes;
        return this;
    }

    /// <summary>
    /// Define a 1D pipe domain with radius (PipeFlow).
    /// </summary>
    /// <param name="length">Pipe length in metres.</param>
    /// <param name="radius">Pipe radius in metres.</param>
    /// <param name="nodes">Number of radial nodes.</param>
    public SimulationBuilder WithGeometry(double length, double radius, int nodes)
    {
        GeomLength = length;
        GeomRadius = radius;
        Nodes = nodes;
        return this;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Cross-Section (Beam)
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Set a rectangular cross-section for beam analysis. Computes I = bh³/12.
    /// </summary>
    public SimulationBuilder WithCrossSection(double width, double height)
    {
        SectionWidth = width;
        SectionHeight = height;
        return this;
    }

    /// <summary>
    /// Set a circular cross-section for beam analysis. Computes I = πr⁴/4.
    /// </summary>
    public SimulationBuilder WithCrossSection(double radius)
    {
        SectionRadius = radius;
        return this;
    }

    /// <summary>
    /// Directly set the second moment of area I (m⁴) for beam analysis.
    /// </summary>
    public SimulationBuilder WithSecondMoment(double I)
    {
        SecondMomentOverride = I;
        return this;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Cylinder (CylinderFlow)
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Define the cylinder obstacle for a CylinderFlow simulation.
    /// </summary>
    /// <param name="centerX">Centre x-position in metres.</param>
    /// <param name="centerY">Centre y-position in metres.</param>
    /// <param name="radius">Cylinder radius in metres.</param>
    public SimulationBuilder WithCylinder(double centerX, double centerY, double radius)
    {
        CylinderCenterX = centerX;
        CylinderCenterY = centerY;
        CylinderRadius = radius;
        return this;
    }

    /// <summary>
    /// Set the uniform freestream inlet velocity U∞ for CylinderFlow.
    /// </summary>
    /// <param name="u">Inlet velocity in m/s.</param>
    public SimulationBuilder WithInletVelocity(double u)
    {
        InletVelocity = u;
        return this;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Boundary Conditions
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Set Dirichlet boundary values on a 2D domain (HeatPlate, ElectricField).
    /// </summary>
    public SimulationBuilder WithBoundary(double top, double bottom, double left, double right)
    {
        TopBC = top;
        BottomBC = bottom;
        LeftBC = left;
        RightBC = right;
        return this;
    }

    /// <summary>
    /// Set beam support type (BeamStress).
    /// </summary>
    public SimulationBuilder WithBoundary(BeamSupport support)
    {
        Support = support;
        return this;
    }

    /// <summary>
    /// Set driving pressure gradient for pipe flow (PipeFlow).
    /// </summary>
    /// <param name="pressureGradient">dP/dx in Pa/m (negative = flow in +x direction).</param>
    public SimulationBuilder WithBoundary(double pressureGradient)
    {
        PressureGradient = pressureGradient;
        return this;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Initial Conditions
    // ═══════════════════════════════════════════════════════════════

    /// <summary>Set a uniform initial value across the domain.</summary>
    public SimulationBuilder WithInitialCondition(double value)
    {
        InitialValue = value;
        HasInitialFunc = false;
        return this;
    }

    /// <summary>Set a spatially varying initial condition via a function (x, y) → value.</summary>
    public SimulationBuilder WithInitialCondition(Func<double, double, double> func)
    {
        InitialFunc = func ?? throw new ArgumentNullException(nameof(func));
        HasInitialFunc = true;
        return this;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Sources
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Add a 2D point source at grid indices (ix, iy) with given intensity.
    /// For HeatPlate: power in W/m³. For ElectricField: charge density ρ in C/m³.
    /// </summary>
    public SimulationBuilder AddSource(int ix, int iy, double value)
    {
        Sources2D.Add((ix, iy, value));
        return this;
    }

    /// <summary>
    /// Add a 1D point load/source at a given position along the domain.
    /// For BeamStress: load in N. For PipeFlow: localised source.
    /// </summary>
    public SimulationBuilder AddSource(double position, double value)
    {
        Sources1D.Add((position, value));
        return this;
    }

    /// <summary>
    /// Set a uniform distributed load along the beam (BeamStress), in N/m.
    /// </summary>
    public SimulationBuilder WithSource(double uniformLoad)
    {
        UniformLoad = uniformLoad;
        return this;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Terminal: Solve / Run
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Execute the simulation and return results.
    /// </summary>
    /// <param name="dt">Time step for transient problems (HeatPlate, PipeFlow).</param>
    /// <param name="steps">Number of time steps (transient) or 0 for steady-state/static.</param>
    /// <param name="maxIterations">Max iterations for iterative solvers (ElectricField).</param>
    /// <param name="tolerance">Convergence tolerance for iterative solvers.</param>
    public SimulationResult Solve(
        double dt = 0.01, int steps = 0,
        int maxIterations = 10000, double tolerance = 1e-6)
    {
        Dt = dt;
        Steps = steps;
        MaxIterations = maxIterations;
        Tolerance = tolerance;

        Validate();

        IMultiphysicsSolver solver = Type switch
        {
            MultiphysicsType.HeatPlate => new HeatPlateSolver(new HeatTransferModel()),
            MultiphysicsType.PipeFlow => new PipeFlowSolver(new ViscousFlowModel()),
            MultiphysicsType.ElectricField => new ElectricFieldSolver(new ElectrostaticModel()),
            MultiphysicsType.BeamStress => new BeamStressSolver(new BeamModel()),
            MultiphysicsType.CylinderFlow => new CylinderFlowSolver(),
            _ => throw new NotSupportedException($"Unknown simulation type: {Type}")
        };

        return solver.Solve(this);
    }

    /// <summary>Alias for <see cref="Solve"/>.</summary>
    public SimulationResult Run(
        double dt = 0.01, int steps = 0,
        int maxIterations = 10000, double tolerance = 1e-6)
        => Solve(dt, steps, maxIterations, tolerance);

    // ═══════════════════════════════════════════════════════════════
    //  Helpers
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Computes the second moment of area from the configured cross-section.
    /// </summary>
    internal double ComputeSecondMoment()
    {
        if (SecondMomentOverride > 0) return SecondMomentOverride;
        if (SectionRadius > 0) return SectionRadius.CircularSecondMoment();
        if (SectionWidth > 0 && SectionHeight > 0)
            return SectionWidth.RectangularSecondMoment(SectionHeight);
        return 0;
    }

    private void Validate()
    {
        if (!Material.HasValue)
            throw new InvalidOperationException("Material not set. Call WithMaterial().");

        switch (Type)
        {
            case MultiphysicsType.HeatPlate:
            case MultiphysicsType.ElectricField:
            case MultiphysicsType.CylinderFlow:
                if (Nx <= 0 || Ny <= 0)
                    throw new InvalidOperationException("2D geometry not set. Call WithGeometry(width, height, nx, ny).");
                break;
            case MultiphysicsType.PipeFlow:
            case MultiphysicsType.BeamStress:
                if (Nodes <= 0)
                    throw new InvalidOperationException("1D geometry not set. Call WithGeometry(length, nodes).");
                break;
        }

        if (Type == MultiphysicsType.CylinderFlow)
        {
            if (CylinderRadius <= 0)
                throw new InvalidOperationException("Cylinder not set. Call WithCylinder(centerX, centerY, radius).");
            if (InletVelocity <= 0)
                throw new InvalidOperationException("Inlet velocity not set. Call WithInletVelocity(u).");
        }

        if (Type == MultiphysicsType.BeamStress)
        {
            if (!Support.HasValue)
                throw new InvalidOperationException("Beam support not set. Call WithBoundary(BeamSupport).");
            if (ComputeSecondMoment() <= 0)
                throw new InvalidOperationException("Cross-section not set. Call WithCrossSection() or WithSecondMoment().");
        }
    }
}
