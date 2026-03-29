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

- [ ] Add iterative Poisson solver (Jacobi/Gauss-Seidel) to `GridOperators` for ∇²φ = f on Grid2D with Dirichlet/Neumann BC — `Numerics/Numerics/Numerics/FiniteDifference/GridOperators.cs`
- [ ] Add `EngineeringMaterial` to `Physics/Materials/Engineering/`
- [ ] Add minimal 1D FEM primitives in `Numerics/Numerics/Numerics/FiniteElement/` (*parallel with Poisson solver*)
  - [ ] `IElement1D` interface — shape functions, local stiffness, local load
  - [ ] `BarElement` (2-node, linear) — axial stress/strain
  - [ ] `BeamElement` (2-node, Hermite cubic) — Euler-Bernoulli bending
  - [ ] `Assembler1D` — local→global stiffness matrix assembly into `Matrix`
  - [ ] `Mesh1D` — 1D mesh from interval subdivision (nodes + elements)
- [ ] Add `EngineeringMaterial` to `Physics/Materials/Engineering/` (*parallel with above*)
  - [ ] `EngineeringMaterial` immutable struct: ThermalConductivity, SpecificHeat, Density, DynamicViscosity, ElectricPermittivity, YoungsModulus, PoissonsRatio
  - [ ] `EngineeringLibrary` with common materials: steel, aluminum, copper, water, air, concrete, glass
  - [ ] Extend `Materials` factory with `Materials.Engineering("Steel")` — add `EngineeringMaterial?` property to `MaterialDescriptor`

## Phase 2 — Core Engine & Individual Solvers

- [ ] Create engine scaffold in `Engines/Multiphysics/` (*depends on Phase 1*)
  - [ ] `MultiphysicsEngine : ISimulationEngine` — orchestrates selected simulation type
  - [ ] `SimulationType` enum: PipeFlow, HeatPlate, ElectricField, BeamStress
  - [ ] `SimulationConfig` — unified config (grid size, material, BCs, time params)
  - [ ] `MultiphysicsEngineBuilder` — fluent builder (pattern from GIS `RiskScenarioBuilder`)
- [ ] `HeatPlateSolver` — 2D heat equation ∂T/∂t = α∇²T + source (*depends on scaffold; uses existing HeatExtensions.HeatEquationRate + Grid2D + Laplacian2D + ITimeStepper*)
  - [ ] Material-driven α = k/(ρ·cp) from `EngineeringMaterial`
  - [ ] BCs: fixed temperature, insulated, mixed
  - [ ] Initial conditions: hot spot, gradient, uniform
- [ ] `PipeFlowSolver` — 1D transient Hagen-Poiseuille via FD (*parallel with HeatPlateSolver*)
  - [ ] Uses `FluidExtensions` NS terms + `GridOperators`
  - [ ] Material viscosity from `EngineeringMaterial`
  - [ ] Optional stretch: 2D lid-driven cavity (projection method + Poisson solver)
- [ ] `ElectricFieldSolver` — 2D Poisson/Laplace ∇²φ = −ρ/ε (*parallel with above*)
  - [ ] Poisson solver from Phase 1, E = −∇φ via `GridOperators.Gradient2D`
  - [ ] Conductor as BC (φ = V on conductor cells), charge distributions as source
  - [ ] Material permittivity from `EngineeringMaterial`
- [ ] `BeamStressSolver` — 1D Euler-Bernoulli via FEM primitives (*parallel with above; depends on Phase 1 FEM*)
  - [ ] Support: cantilever, simply supported, fixed-fixed
  - [ ] Loads: point loads, distributed loads, moments
  - [ ] Output: deflection curve, bending moment, shear force, stress distribution
  - [ ] E, I from `EngineeringMaterial` + cross-section dimensions

## Phase 3 — Export & Visualization Support

- [ ] `MultiphysicsBinaryExporter` — Unity-compatible binary format (same MAGIC+header pattern as GIS `UnityBinaryExporter`) (*depends on Phase 2*)
- [ ] `MultiphysicsJsonExporter` — JSON format for web visualization
- [ ] Snapshot system (*parallel with exporters*)
  - [ ] `FieldSnapshot` — time-stamped 2D scalar field (like GIS `GridSnapshot`)
  - [ ] `BeamSnapshot` — time-stamped 1D arrays (deflection, moment, shear, stress)
  - [ ] `SimulationTimeline` — ordered collection of snapshots with interpolation

## Phase 4 — ML & Monte Carlo Integration

- [ ] `MultiphysicsMonteCarloModel : IMonteCarloModel` — runs solver N times with sampled parameters (*depends on Phase 2; pattern from GIS `PlumeMonteCarloModel`*)
  - [ ] `ParameterVariation` — ranges for material properties, BCs, loads
  - [ ] Output: scenario matrix, percentile maps
- [ ] `SurrogateTrainer` — trains regression model (Linear/SVR/MLP) on MC scenario data (*depends on MC model*)
  - [ ] Predict simulation output without re-running solver
- [ ] `MultiphysicsClusterAnalyzer` — wraps `ClusteringExperiment` on scenario matrix (*parallel with surrogate; pattern from GIS `ScenarioClusterAnalyzer`*)
  - [ ] Identifies representative scenarios from MC ensemble

## Phase 5 — Tests & Documentation

- [ ] Unit tests with analytical validation
  - [ ] Heat plate: steady-state convergence vs known Fourier solutions
  - [ ] Pipe flow: velocity profile vs analytical Poiseuille
  - [ ] Electric field: capacitor/line charge vs analytical E-field
  - [ ] Beam stress: cantilever deflection vs PL³/3EI
  - [ ] FEM primitives: stiffness matrix assembly, known beam cases
  - [ ] MC integration: scenario matrix dimensions, clustering produces valid results
  - [ ] Export: binary format header parsing, dimension checks
- [ ] Update `Engines/README.md` — add Multiphysics engine section
- [ ] Update `Numerics/README.md` — add FiniteElement section & Poisson solver
- [ ] Update `Physics/README.md` — add EngineeringMaterial section

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
-  `Numerics/Numerics/Physics/SolidExtensions.cs` — Beam stress

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