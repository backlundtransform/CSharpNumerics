# Multiphysics Engine Roadmap

## Overview

A new `Engines/Multiphysics/` simulation engine enabling four student-friendly PDE simulations on simple geometries, designed for visualization in Unity or web platforms. Built on existing finite difference infrastructure with new engineering materials, Unity/web export, and optional ML/Monte Carlo integration following the GIS engine pattern.

**Target complexity:** ~20–25 files, ~3.5–5K LOC (comparable to Game/Quantum engines).

### LBM Assessment

Lattice Boltzmann Methods were evaluated and **rejected**. FD (Grid2D, GridOperators, 4 time steppers) already covers all four simulation types. LBM would require ~1–2K LOC of new infrastructure (D2Q9 distributions, BGK collision, streaming) duplicating existing FD capability. LBM strengths (parallelism, complex boundaries) are unnecessary for simple student geometries.

### Simulation Types

| Simulation | PDE | Method | Dimension |
|------------|-----|--------|-----------|
| Heat plate | ∂T/∂t = α∇²T | FD (Grid2D + Laplacian2D + ITimeStepper) | 2D |
| Pipe flow | Navier-Stokes (simplified) | FD (1D transient Hagen-Poiseuille) | 1D (2D stretch) |
| Electric field | ∇²φ = −ρ/ε | FD (Poisson solver + Gradient2D) | 2D |
| Beam stress | EIu⁗ = q | 1D FEM (Euler-Bernoulli) | 1D |

---

## Phase 1 — Foundations & Shared Infrastructure

- [x] Add iterative Poisson solver (Gauss-Seidel) to `GridOperators` for ∇²φ = f on Grid2D with Dirichlet BC — `GridOperators.SolvePoisson2D()`
- [x] Add `SolidExtensions` to `Physics/` — Hooke's law, second moment of area, Euler-Bernoulli beam equation, analytical deflections
- [x] Add minimal 1D FEM primitives in `Numerics/Numerics/Numerics/FiniteElement/`
  - [x] `IElement1D` interface — shape functions, local stiffness, local load
  - [x] `BarElement` (2-node, linear) — axial stress/strain
  - [x] `BeamElement` (2-node, Hermite cubic) — Euler-Bernoulli bending
  - [x] `Assembler1D` — local→global stiffness matrix assembly with Gaussian elimination solver
  - [x] `Mesh1D` — 1D mesh from interval subdivision (nodes + elements)
- [x] Add `EngineeringMaterial` to `Physics/Materials/Engineering/`
  - [x] `EngineeringMaterial` immutable struct: ThermalConductivity, SpecificHeat, Density, DynamicViscosity, ElectricPermittivity, YoungsModulus, PoissonsRatio
  - [x] `EngineeringLibrary` with common materials: steel, aluminum, copper, water, air, concrete, glass
  - [ ] Extend `Materials` factory with `Materials.Engineering("Steel")` — add `EngineeringMaterial?` property to `MaterialDescriptor`

## Phase 2 — Core Engine & Individual Solvers

- [x] Create engine scaffold in `Engines/Multiphysics/` with **fluent API** (*pattern from GIS `RiskScenario`*)
  - [x] `SimulationType` static entry point — `SimulationType.Create(MultiphysicsType)`
  - [x] `SimulationBuilder` — fluent builder: `.WithMaterial()` → `.WithGeometry()` → `.WithBoundary()` → `.WithInitialCondition()` → `.AddSource()` → `.Solve()` / `.Run()`
  - [x] `MultiphysicsType` enum: HeatPlate, PipeFlow, ElectricField, BeamStress
  - [x] `BeamSupport` enum: Cantilever, SimplySupported, FixedFixed
  - [x] `SimulationResult` — unified result type with 2D fields, 1D arrays, beam-specific outputs, E-field vectors
  - [x] `IMultiphysicsSolver` — internal solver interface
- [x] `HeatPlateSolver` — 2D heat equation ∂T/∂t = α∇²T + source (Grid2D + Laplacian2D + forward Euler)
  - [x] Material-driven α = k/(ρ·cp) from `EngineeringMaterial`
  - [x] BCs: fixed temperature (Dirichlet edges)
  - [x] Initial conditions: uniform or function-based
- [x] `PipeFlowSolver` — 1D transient Hagen-Poiseuille via FD (cylindrical Laplacian)
  - [x] Material viscosity ν from `EngineeringMaterial.KinematicViscosity`
  - [x] Symmetry BC at r=0, no-slip at r=R
  - [ ] Optional stretch: 2D lid-driven cavity (projection method + Poisson solver)
- [x] `ElectricFieldSolver` — 2D Poisson/Laplace ∇²φ = −ρ/ε
  - [x] Poisson solver from Phase 1, E = −∇φ via `GridOperators.Gradient2D`
  - [x] Conductor as BC (φ = V on edge cells), charge distributions as source
  - [x] Material permittivity from `EngineeringMaterial`
- [x] `BeamStressSolver` — 1D Euler-Bernoulli via analytical `SolidExtensions` formulas
  - [x] Support: cantilever, simply supported, fixed-fixed
  - [x] Loads: point loads, distributed loads
  - [x] Output: deflection curve, bending moment, shear force, stress distribution
  - [x] EI from `EngineeringMaterial.YoungsModulus` + cross-section (rectangular, circular, or custom I)

## Phase 3 — Export & Visualization Support

- [x] `MultiphysicsBinaryExporter` — Unity-compatible binary format (MPHY magic + header + float layers)
  - [x] 2D timeline export: per-step float[Nx*Ny] layers
  - [x] 1D beam export: positions + 4 layers (deflection, moment, shear, stress)
  - [x] Reader: `ReadHeader()` and `Read()` for round-trip loading
- [x] `MultiphysicsJsonExporter` — JSON format for web visualization
  - [x] `SimulationResult` direct export (all 4 simulation types)
  - [x] `SimulationTimeline` export with per-step snapshots
  - [x] `FieldSnapshot` single-frame export
  - [x] `BeamSnapshot` export with all curves
  - [x] `ExportMetadata` support (simulation, unit, density)
  - [x] `Save()` to file for all types
- [x] Snapshot system
  - [x] `FieldSnapshot` — time-stamped 2D scalar field with flat indexing, Min/Max, ToArray round-trip
  - [x] `BeamSnapshot` — immutable 1D arrays (deflection, moment, shear, stress) with `FromResult()` factory
  - [x] `SimulationTimeline` — ordered collection with `FromResult()`, `InterpolateAt()` linear blending

## Phase 4 — ML & Monte Carlo Integration

- [x] `MultiphysicsMonteCarloModel : IMonteCarloModel` — runs solver N times with sampled parameters (*depends on Phase 2; pattern from GIS `PlumeMonteCarloModel`*)
  - [x] `ParameterVariation` — ranges for material properties, BCs, loads
  - [x] Output: scenario matrix, percentile maps
- [x] `SurrogateTrainer` — trains regression model (Linear/SVR/MLP) on MC scenario data (*depends on MC model*)
  - [x] Predict simulation output without re-running solver
- [x] `MultiphysicsClusterAnalyzer` — wraps `ClusteringExperiment` on scenario matrix (*parallel with surrogate; pattern from GIS `ScenarioClusterAnalyzer`*)
  - [x] Identifies representative scenarios from MC ensemble

## Phase 5 — Tests & Documentation

- [x] Unit tests with analytical validation
  - [x] Heat plate: steady-state convergence vs known Fourier solutions
  - [x] Pipe flow: velocity profile vs analytical Poiseuille
  - [x] Electric field: capacitor/line charge vs analytical E-field
  - [x] Beam stress: cantilever deflection vs PL³/3EI
  - [x] FEM primitives: stiffness matrix assembly, known beam cases (bar axial, cantilever, simply supported)
  - [x] MC integration: scenario matrix dimensions, clustering produces valid results
  - [x] Export: binary format header parsing, dimension checks
- [x] Update `Engines/Multiphysics/README.md` — add Multiphysics engine section
- [x] Update `Numerics/README.md` — add Poisson solver section
- [x] Update `Physics/README.md` — add EngineeringMaterial section

---

## Key Files

### Existing — to modify
- `Numerics/Numerics/Numerics/FiniteDifference/GridOperators.cs` — add Poisson solver
- `Numerics/Numerics/Physics/Materials/Materials.cs` — extend factory with Engineering path

### Existing — to reuse as patterns
- `Numerics/Numerics/Engines/GIS/Simulation/PlumeMonteCarloModel.cs` — MC integration
- `Numerics/Numerics/Engines/GIS/Analysis/ScenarioClusterAnalyzer.cs` — clustering
- `Numerics/Numerics/Engines/GIS/Export/UnityBinaryExporter.cs` — binary export format
- `Numerics/Numerics/Engines/Common/ISimulationEngine.cs` — engine interface
- `Numerics/Numerics/Physics/HeatExtensions.cs` — HeatEquationRate()
- `Numerics/Numerics/Physics/FluidExtensions.cs` — NS residuals, pipe flow
- `Numerics/Numerics/Physics/ElectroMagneticFieldExtensions.cs` — E-field calculations


### New — to create
- `Numerics/Numerics/Numerics/FiniteElement/` — IElement1D, BarElement, BeamElement, Assembler1D, Mesh1D (~6 files)
- `Numerics/Numerics/Physics/Materials/Engineering/` — EngineeringMaterial, EngineeringLibrary (~2 files)
- `Numerics/Numerics/Engines/Multiphysics/` — engine, 4 solvers, config, export, snapshots, MC, ML (~15 files)


## Decisions

- **FD over LBM** — existing infrastructure sufficient; LBM complexity unjustified for simple geometries
- **1D FEM primitives in Numerics** — minimal (Element, Assembly, Mesh), beam solver logic in engine
- **Export like GIS** — UnityBinaryExporter + JSON for web
- **Engineering materials** — new struct + library extending existing MaterialDescriptor
- **Pipe flow scope** — 1D Hagen-Poiseuille primary, 2D lid-driven cavity as stretch goal
- **Beam analysis** — static only initially; dynamic vibration deferred
- **RL environment** — deferred to follow-up roadmap