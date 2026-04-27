# Cylinder Flow (Kármán Vortex Street) Roadmap

## Overview

Add a 2D incompressible Navier-Stokes solver to `Engines/Multiphysics/` using Chorin's projection method on `Grid2D`. The headline scenario is **vortex shedding around a circular cylinder** (Kármán vortex street) — a classic CFD benchmark that complements the existing 1D Hagen-Poiseuille pipe flow.

All required FD infrastructure already exists in `GridOperators`: `Laplacian2D`, `Gradient2D`, `Divergence2D`, `Advection2D`, `SolvePoisson2D`. This feature connects them into a time-stepping loop with immersed-boundary (mask) handling for the cylinder.

**Target complexity:** ~300–450 LOC solver + ~100 LOC builder/result extensions + ~150 LOC tests.

### Algorithm — Chorin Projection (Fractional Step)

Each time step:

1. **Predict** — compute intermediate velocity v* from advection + diffusion:
   v* = vⁿ + Δt · (−(v·∇)v + ν∇²v)
2. **Pressure Poisson** — enforce incompressibility:
   ∇²p = (∇·v*) / Δt
3. **Correct** — project onto divergence-free field:
   vⁿ⁺¹ = v* − Δt · ∇p
4. **Enforce BCs** — no-slip on cylinder mask, inlet, outlet, walls.

### Boundary Conditions

| Boundary | Type |
|----------|------|
| Inlet (left) | Uniform velocity U∞ (Dirichlet) |
| Outlet (right) | Zero-gradient ∂v/∂x = 0 (Neumann) |
| Top / Bottom | Free-slip or periodic |
| Cylinder | No-slip v = 0 (immersed boundary mask) |

### Validation Targets

| Quantity | Expected | Source |
|----------|----------|--------|
| Steady symmetric wake | Re < 47 | Analytical / Stokes |
| Onset of shedding | Re ≈ 47 | Empirical |
| Strouhal number | St ≈ 0.2 at Re = 100–200 | Roshko (1954) |
| Drag coefficient | Cd ≈ 1.3–1.5 at Re = 100 | Literature |

---

## Phase 1 — Solver Core

- [x] Add `CylinderFlow` to `MultiphysicsType` enum
- [x] Create `CylinderFlowSolver : IMultiphysicsSolver` in `Engines/Multiphysics/Solvers/`
  - [x] Build `Grid2D` from builder geometry (width × height, Nx × Ny)
  - [x] Generate cylinder mask: `bool[]` — true for cells inside cylinder (centre + radius from builder)
  - [x] Implement projection method time loop:
    - [x] Step 1 — advection via `GridOperators.Advection2D` (upwind, already exists)
    - [x] Step 2 — diffusion via `GridOperators.Laplacian2D`
    - [x] Step 3 — divergence of v* via `GridOperators.Divergence2D`
    - [x] Step 4 — pressure Poisson via `GridOperators.SolvePoisson2D`
    - [x] Step 5 — velocity correction via `GridOperators.Gradient2D`
  - [x] Enforce boundary conditions each step: inlet (Dirichlet), outlet (Neumann copy), top/bottom (free-slip), cylinder mask (no-slip)
  - [x] CFL-based adaptive Δt clamping: Δt ≤ min(Δx/|u|, Δx²/(4ν))
  - [x] Store vx, vy, p snapshots for timeline

## Phase 2 — Builder & Result Extensions

- [x] Extend `SimulationBuilder` with cylinder-flow configuration:
  - [x] `.WithCylinder(double centerX, double centerY, double radius)` — cylinder position and size
  - [x] `.WithInletVelocity(double u)` — freestream velocity
  - [x] Reuse existing `.WithMaterial()` (ν from `EngineeringMaterial.KinematicViscosity`) and `.WithGeometry(width, height, nx, ny)`
- [x] Extend `SimulationResult` with 2D velocity fields:
  - [x] `Vx` — double[,] x-velocity field
  - [x] `Vy` — double[,] y-velocity field
  - [x] `Pressure` — double[,] pressure field
  - [x] `CylinderMask` — bool[,] for visualization
- [x] Wire `CylinderFlow` case into `SimulationBuilder.Solve()` dispatch
- [ ] Add `IViscousFlowModel` delegation for physics calculations (follows existing PipeFlowSolver pattern)

## Phase 3 — Post-Processing & Diagnostics

- [x] Compute vorticity field: ω = ∂vy/∂x − ∂vx/∂y (via `GridOperators.Gradient2D` on vx, vy)
- [x] Add `Vorticity` double[,] to `SimulationResult`
- [x] Compute drag and lift coefficients from pressure + viscous stress integration around cylinder
  - [x] `DragCoefficient` and `LiftCoefficient` on `SimulationResult`
- [x] Strouhal number estimation: detect shedding frequency from lift coefficient time series
  - [ ] Use existing `PeakFitting.FindPeaks` or `LombScarglePeriodogram` from Statistics
  - [x] `StrouhalNumber` on `SimulationResult`

## Phase 4 — Export & Snapshot Support

- [ ] Extend `FieldSnapshot.FromResult()` to handle CylinderFlow (vx, vy, vorticity layers)
- [x] Extend `MultiphysicsBinaryExporter` — write vx, vy, p, vorticity as separate layers per timestep
- [x] Extend `MultiphysicsJsonExporter` — include velocity components and vorticity

## Phase 5 — Tests & Documentation

- [x] Unit tests in `NumericTest/MultiphysicsTests.cs` (or new `CylinderFlowTests.cs`):
  - [ ] Poiseuille channel flow (no cylinder): verify parabolic profile u(y) = (ΔP/(2μL))·y·(H−y)
  - [ ] Low-Re cylinder (Re ≈ 20): symmetric wake, no shedding
  - [ ] Moderate-Re cylinder (Re ≈ 100): verify Strouhal ≈ 0.16–0.22
  - [x] Mass conservation: verify ∇·v ≈ 0 everywhere after projection
  - [x] CFL stability: verify solver rejects Δt violating CFL condition
- [x] Update `Engines/Multiphysics/README.md` with CylinderFlow section and usage examples
- [ ] Update `Engines/README.md` if simulation type table needs a new row

---

## Key Files

### Existing — to modify
- `Engines/Multiphysics/Enums/MultiphysicsType.cs` — add `CylinderFlow`
- `Engines/Multiphysics/SimulationBuilder.cs` — add cylinder geometry + inlet velocity methods
- `Engines/Multiphysics/SimulationResult.cs` — add Vx, Vy, Pressure, Vorticity, CylinderMask, drag/lift/Strouhal
- `Engines/Multiphysics/SimulationType.cs` — add dispatch for CylinderFlow
- `Engines/Multiphysics/Export/MultiphysicsBinaryExporter.cs` — handle new layers
- `Engines/Multiphysics/Export/MultiphysicsJsonExporter.cs` — handle new fields
- `Engines/Multiphysics/Snapshots/FieldSnapshot.cs` — handle CylinderFlow

### Existing — to reuse (no changes needed)
- `Numerics/FiniteDifference/GridOperators.cs` — Laplacian2D, Gradient2D, Divergence2D, Advection2D, SolvePoisson2D
- `Numerics/FiniteDifference/Grid2D.cs` — computational grid
- `Physics/FluidDynamics/Interfaces/IViscousFlowModel.cs` — physics delegation
- `Physics/Materials/Engineering/EngineeringMaterial.cs` — KinematicViscosity
- `Statistics/TimeSeriesAnalysis/LombScarglePeriodogram.cs` — Strouhal estimation (optional)
- `Statistics/TimeSeriesAnalysis/PeakFitting.cs` — shedding frequency detection (optional)

### New files
- `Engines/Multiphysics/Solvers/CylinderFlowSolver.cs` — projection method solver

### Usage Example

```csharp
using CSharpNumerics.Engines.Multiphysics;
using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Physics.Materials.Engineering;

var result = SimulationType.Create(MultiphysicsType.CylinderFlow)
    .WithMaterial(EngineeringLibrary.Water)
    .WithGeometry(width: 2.0, height: 0.8, nx: 200, ny: 80)
    .WithCylinder(centerX: 0.4, centerY: 0.4, radius: 0.05)
    .WithInletVelocity(0.1)                 // U∞ = 0.1 m/s
    .WithTimeStep(dt: 0.001, steps: 5000)
    .Solve();

// Inspect results
double[,] vorticity = result.Vorticity;     // ω field for visualization
double cd = result.DragCoefficient;          // ≈ 1.3–1.5 at Re ≈ 100
double st = result.StrouhalNumber;           // ≈ 0.2
bool[,] mask = result.CylinderMask;         // for rendering the obstacle
```
