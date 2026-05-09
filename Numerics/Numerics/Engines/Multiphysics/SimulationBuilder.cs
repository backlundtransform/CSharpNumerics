using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Numerics.FiniteElement.Enums;
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

    // ── 3D Geometry (HeatBlock3D) ────────────────────────────────
    internal double GeomDepth { get; private set; }
    internal int Nz { get; private set; }

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
    internal double Thickness { get; private set; } = 1.0;
    internal PlaneType PlaneType { get; private set; } = PlaneType.PlaneStress;

    // ── FluidDiffusion3D ─────────────────────────────────────────
    internal double DiffusionCoefficient { get; private set; }
    internal Func<double, double, double, (double, double, double)> VelocityField3D { get; private set; }
    internal bool HasVelocityField3D { get; private set; }

    // ── Boundary conditions ──────────────────────────────────────
    internal double TopBC { get; private set; }
    internal double BottomBC { get; private set; }
    internal double LeftBC { get; private set; }
    internal double RightBC { get; private set; }
    internal double FrontBC { get; private set; }
    internal double BackBC { get; private set; }
    internal BeamSupport? Support { get; private set; }
    internal double PressureGradient { get; private set; }

    // ── Per-face thermal BC overrides (default = Dirichlet) ──────
    internal FaceBoundaryCondition TopFaceBC { get; private set; }
    internal FaceBoundaryCondition BottomFaceBC { get; private set; }
    internal FaceBoundaryCondition LeftFaceBC { get; private set; }
    internal FaceBoundaryCondition RightFaceBC { get; private set; }
    internal FaceBoundaryCondition FrontFaceBC { get; private set; }
    internal FaceBoundaryCondition BackFaceBC { get; private set; }
    internal bool HasConvectionBC { get; private set; }

    // ── Initial condition ────────────────────────────────────────
    internal double InitialValue { get; private set; }
    internal Func<double, double, double> InitialFunc { get; private set; }
    internal bool HasInitialFunc { get; private set; }
    internal Func<double, double, double, double> InitialFunc3D { get; private set; }
    internal bool HasInitialFunc3D { get; private set; }

    // ── Sources ──────────────────────────────────────────────────
    internal List<(int ix, int iy, double value)> Sources2D { get; } = new();
    internal List<(int ix, int iy, int iz, double value)> Sources3D { get; } = new();
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
    /// Define a 3D rectangular domain (HeatBlock3D).
    /// </summary>
    /// <param name="width">Physical width (x) in metres.</param>
    /// <param name="height">Physical height (y) in metres.</param>
    /// <param name="depth">Physical depth (z) in metres.</param>
    /// <param name="nx">Grid cells in x.</param>
    /// <param name="ny">Grid cells in y.</param>
    /// <param name="nz">Grid cells in z.</param>
    public SimulationBuilder WithGeometry3D(double width, double height, double depth, int nx, int ny, int nz)
    {
        GeomWidth = width;
        GeomHeight = height;
        GeomDepth = depth;
        Nx = nx;
        Ny = ny;
        Nz = nz;
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

    /// <summary>
    /// Set a prescribed 3D velocity field for advection-diffusion (FluidDiffusion3D).
    /// The function maps physical coordinates (x, y, z) to velocity (vx, vy, vz).
    /// </summary>
    public SimulationBuilder WithVelocityField3D(Func<double, double, double, (double, double, double)> field)
    {
        VelocityField3D = field ?? throw new ArgumentNullException(nameof(field));
        HasVelocityField3D = true;
        return this;
    }

    /// <summary>
    /// Set the scalar diffusion coefficient D (m²/s) for advection-diffusion (FluidDiffusion3D).
    /// </summary>
    public SimulationBuilder WithDiffusionCoefficient(double D)
    {
        DiffusionCoefficient = D;
        return this;
    }

    /// <summary>
    /// Set the constitutive assumption for plane elasticity problems.
    /// </summary>
    public SimulationBuilder WithPlaneType(PlaneType planeType)
    {
        PlaneType = planeType;
        return this;
    }

    /// <summary>
    /// Set the section thickness used by 2D plane-stress/plane-strain elements.
    /// </summary>
    public SimulationBuilder WithThickness(double thickness)
    {
        Thickness = thickness;
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
    /// Set Dirichlet boundary values on all six faces of a 3D domain.
    /// </summary>
    public SimulationBuilder WithBoundary3D(
        double top, double bottom, double left, double right, double front, double back)
    {
        TopBC = top;
        BottomBC = bottom;
        LeftBC = left;
        RightBC = right;
        FrontBC = front;
        BackBC = back;
        return this;
    }

    /// <summary>
    /// Set convective (Robin) boundary conditions on all four faces of a 2D domain.
    /// −k ∂T/∂n = h(T − T∞) at each boundary face.
    /// </summary>
    /// <param name="heatTransferCoefficient">Convective coefficient h in W/(m²·K).</param>
    /// <param name="ambientTemperature">Far-field fluid temperature T∞ in K.</param>
    public SimulationBuilder WithConvectionBoundary(double heatTransferCoefficient, double ambientTemperature)
    {
        var bc = FaceBoundaryCondition.Convection(heatTransferCoefficient, ambientTemperature);
        TopFaceBC = bc;
        BottomFaceBC = bc;
        LeftFaceBC = bc;
        RightFaceBC = bc;
        HasConvectionBC = true;
        return this;
    }

    /// <summary>
    /// Set convective (Robin) boundary conditions on all six faces of a 3D domain.
    /// −k ∂T/∂n = h(T − T∞) at each boundary face.
    /// </summary>
    /// <param name="heatTransferCoefficient">Convective coefficient h in W/(m²·K).</param>
    /// <param name="ambientTemperature">Far-field fluid temperature T∞ in K.</param>
    public SimulationBuilder WithConvectionBoundary3D(double heatTransferCoefficient, double ambientTemperature)
    {
        var bc = FaceBoundaryCondition.Convection(heatTransferCoefficient, ambientTemperature);
        TopFaceBC = bc;
        BottomFaceBC = bc;
        LeftFaceBC = bc;
        RightFaceBC = bc;
        FrontFaceBC = bc;
        BackFaceBC = bc;
        HasConvectionBC = true;
        return this;
    }

    /// <summary>
    /// Override the boundary condition on a single 2D face. Use after
    /// <see cref="WithBoundary"/> or <see cref="WithConvectionBoundary"/>
    /// to mix Dirichlet and Robin conditions.
    /// </summary>
    /// <param name="face">Which face to set (Top, Bottom, Left, Right).</param>
    /// <param name="bc">The boundary condition for that face.</param>
    public SimulationBuilder WithFaceBoundary(string face, FaceBoundaryCondition bc)
    {
        switch (face.ToLowerInvariant())
        {
            case "top": TopFaceBC = bc; break;
            case "bottom": BottomFaceBC = bc; break;
            case "left": LeftFaceBC = bc; break;
            case "right": RightFaceBC = bc; break;
            case "front": FrontFaceBC = bc; break;
            case "back": BackFaceBC = bc; break;
            default: throw new ArgumentException($"Unknown face: {face}");
        }
        HasConvectionBC = true;
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

    /// <summary>Set a spatially varying 3D initial condition via a function (x, y, z) → value.</summary>
    public SimulationBuilder WithInitialCondition3D(Func<double, double, double, double> func)
    {
        InitialFunc3D = func ?? throw new ArgumentNullException(nameof(func));
        HasInitialFunc3D = true;
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
    /// Add a 3D point source at grid indices (ix, iy, iz).
    /// For HeatBlock3D: volumetric power in W/m³.
    /// </summary>
    public SimulationBuilder AddSource3D(int ix, int iy, int iz, double value)
    {
        Sources3D.Add((ix, iy, iz, value));
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
            MultiphysicsType.FluidFlow2D => new FluidFlow2DSolver(),
            MultiphysicsType.MagneticField => new MagneticFieldSolver(),
            MultiphysicsType.PlaneStress => new PlaneStressSolver(),
            MultiphysicsType.HeatBlock3D => new HeatBlock3DSolver(new HeatTransferModel()),
            MultiphysicsType.FluidDiffusion3D => new FluidDiffusion3DSolver(),
            MultiphysicsType.CylinderFlow3D => new CylinderFlow3DSolver(),
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
            case MultiphysicsType.FluidFlow2D:
            case MultiphysicsType.MagneticField:
            case MultiphysicsType.PlaneStress:
                if (Nx <= 0 || Ny <= 0)
                    throw new InvalidOperationException("2D geometry not set. Call WithGeometry(width, height, nx, ny).");
                break;
            case MultiphysicsType.HeatBlock3D:
            case MultiphysicsType.FluidDiffusion3D:
            case MultiphysicsType.CylinderFlow3D:
                if (Nx <= 0 || Ny <= 0 || Nz <= 0)
                    throw new InvalidOperationException("3D geometry not set. Call WithGeometry3D(width, height, depth, nx, ny, nz).");
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

        if (Type == MultiphysicsType.CylinderFlow3D)
        {
            if (CylinderRadius <= 0)
                throw new InvalidOperationException("Cylinder not set. Call WithCylinder(centerX, centerY, radius).");
            if (InletVelocity <= 0)
                throw new InvalidOperationException("Inlet velocity not set. Call WithInletVelocity(u).");
        }

        if (Type == MultiphysicsType.FluidDiffusion3D && DiffusionCoefficient <= 0)
            throw new InvalidOperationException("Diffusion coefficient not set. Call WithDiffusionCoefficient(D).");

        if (Type == MultiphysicsType.BeamStress)
        {
            if (!Support.HasValue)
                throw new InvalidOperationException("Beam support not set. Call WithBoundary(BeamSupport).");
            if (ComputeSecondMoment() <= 0)
                throw new InvalidOperationException("Cross-section not set. Call WithCrossSection() or WithSecondMoment().");
        }

        if (Type == MultiphysicsType.PlaneStress && Thickness <= 0)
            throw new InvalidOperationException("Thickness must be positive. Call WithThickness(thickness).");
    }
}
