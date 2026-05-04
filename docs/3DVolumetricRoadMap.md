# 3D Volumetric Computing Roadmap

## Overview

Extend CSharpNumerics with **3D finite difference** primitives, **3D multiphysics solvers**, and **3D volumetric contamination spread** in the GIS engine. The work builds on existing 2D infrastructure:

- `Grid2D` + `GridOperators` (Laplacian2D, Advection2D, Gradient2D, Divergence2D, SolvePoisson2D)
- `Engines/Multiphysics/` solvers (HeatPlate, FluidFlow2D, CylinderFlow)
- `GeoGrid` — already supports 3D indexing (Nx × Ny × Nz) but lacks FD operators
- `WaterContaminationSimulator` — currently 1D along river network, no volumetric depth

### Motivation

| Question | Current answer | After 3D |
|----------|----------------|----------|
| "How does heat diffuse in a 3D block?" | N/A (2D only) | `HeatBlock3D` solver |
| "How does a contaminant spread in volume over time?" | 1D along river | Full 3D advection-diffusion in lake/ocean |
| "Which depths are affected?" | N/A | Per-layer concentration + depth profiles |
| "What does flow look like around a 3D cylinder?" | 2D cross-section | 3D velocity field with axial variation |

### Architecture

```
Numerics/FiniteDifference/
├── Grid2D.cs               ← existing
├── Grid3D.cs               ← NEW
├── GridOperators.cs         ← existing (1D + 2D)
├── GridOperators3D.cs       ← NEW (3D stencils)
└── BoundaryCondition.cs     ← reuse (Dirichlet, Neumann, Periodic)

Engines/Multiphysics/
├── Enums/MultiphysicsType.cs   ← extend with 3D types
├── Solvers/
│   ├── HeatPlateSolver.cs      ← existing 2D
│   ├── HeatBlock3DSolver.cs    ← NEW
│   ├── FluidFlow2DSolver.cs    ← existing 2D
│   ├── FluidDiffusion3DSolver.cs  ← NEW
│   └── CylinderFlow3DSolver.cs    ← NEW (stretch goal)
└── SimulationBuilder.cs        ← extend with Nz/depth

Engines/GIS/
├── Grid/GeoGrid.cs             ← already 3D
├── Spread/
│   ├── WaterContamination/     ← existing 1D river
│   └── VolumetricContamination/  ← NEW 3D lake/ocean
│       ├── VolumetricContaminationSimulator.cs
│       ├── VolumetricContaminationParameters.cs
│       ├── VolumetricContaminationResult.cs
│       ├── DepthProfile.cs
│       └── Enums/
│           └── ContaminationCellState3D.cs
```

### Cross-section dependencies

```
VolumetricContamination  →  GridOperators3D / Grid3D  +  GeoGrid (already 3D)
Multiphysics 3D solvers  →  GridOperators3D / Grid3D
Grid3D + GridOperators3D →  VectorN, BoundaryCondition (existing)
```

---

## Phase 1 — 3D Finite Difference Primitives

Foundation: `Grid3D` and all 3D FD operators. Everything else depends on this.

### Grid3D

- [x] `Grid3D` class in `Numerics/FiniteDifference/Grid3D.cs`
  - [x] Constructor: `Grid3D(nx, ny, nz, dx, dy, dz)`
  - [x] `Nx`, `Ny`, `Nz` — cell counts
  - [x] `Dx`, `Dy`, `Dz` — uniform spacing
  - [x] `Length` → `Nx × Ny × Nz`
  - [x] `Index(ix, iy, iz)` → int (row-major: ix fastest, then iy, then iz)
  - [x] `Index3D(flatIndex)` → `(int ix, int iy, int iz)`
  - [x] `ToVector(double[,,])` → `VectorN` — pack 3D array to flat
  - [x] `ToArray(VectorN)` → `double[,,]` — unpack flat to 3D array
  - [x] `Initialize(Func<double, double, double, double>)` → `VectorN` — evaluate f(x,y,z) on cell centres
  - [x] `Zeros()` → `VectorN`

### GridOperators3D

- [x] `GridOperators3D` static class in `Numerics/FiniteDifference/GridOperators3D.cs`
  - [x] `Laplacian3D(u, grid, bc)` → `VectorN` — 7-point stencil: ∂²u/∂x² + ∂²u/∂y² + ∂²u/∂z²
  - [x] `Gradient3D(u, grid, bc)` → `(VectorN dux, VectorN duy, VectorN duz)` — central differences
  - [x] `Divergence3D(fx, fy, fz, grid, bc)` → `VectorN` — ∂Fx/∂x + ∂Fy/∂y + ∂Fz/∂z
  - [x] `Advection3D(u, vx, vy, vz, grid, bc)` → `VectorN` — first-order upwind v·∇u
  - [x] `SolvePoisson3D(rhs, grid, boundaryMask, boundaryValues, maxIter, tol)` → `(VectorN, int)` — Gauss-Seidel iterative

### Tests

- [x] `Grid3D_Index_RoundTrips` — verify Index(ix,iy,iz) ↔ Index3D(flat) for all cells
- [x] `Grid3D_ToVector_ToArray_RoundTrips` — pack/unpack identity
- [x] `Grid3D_Initialize_SetsCorrectPositions`
- [x] `Laplacian3D_Quadratic_IsConstant` — ∇²(x²+y²+z²) = 6
- [x] `Laplacian3D_Linear_IsZero` — ∇²(ax+by+cz) = 0
- [x] `Gradient3D_Linear_ReturnsCoefficients` — ∇(ax+by+cz) = (a,b,c)
- [x] `Divergence3D_IdentityField_IsThree` — ∇·(x,y,z) = 3
- [x] `Advection3D_UniformFlow_ShiftsProfile`
- [x] `SolvePoisson3D_KnownSolution_Converges`
- [x] `Laplacian3D_Periodic_Wraps`
- [x] `Laplacian3D_Neumann_MirrorsInterior`

---

## Phase 2 — 3D Heat Equation (Multiphysics)

First 3D solver: heat diffusion in a solid block. Validates the FD primitives in a physics context.

### Equation

∂T/∂t = α ∇²T + S

where α = k/(ρ·cp) is thermal diffusivity, S is a volumetric source term.

### Implementation

- [x] Add `HeatBlock3D` to `MultiphysicsType` enum
- [x] Extend `SimulationBuilder` with 3D geometry:
  - [x] `.WithGeometry3D(width, height, depth, nx, ny, nz)` sets `GeomWidth`, `GeomHeight`, `GeomDepth`, `Nx`, `Ny`, `Nz`
  - [x] `.WithBoundary3D(top, bottom, left, right, front, back)` — six Dirichlet faces
  - [x] `.AddSource3D(ix, iy, iz, value)` — point source in volume
- [x] Extend `SimulationResult` with 3D fields:
  - [x] `Field3D` → `double[,,]` (Nx × Ny × Nz)
  - [x] `Timeline3D` → `List<double[,,]>` — temporal snapshots
  - [x] `SliceXY(iz)` → `double[,]` — extract horizontal slice at given z
  - [x] `SliceXZ(iy)` → `double[,]` — extract vertical slice at given y
  - [x] `SliceYZ(ix)` → `double[,]` — extract vertical slice at given x
- [x] Create `HeatBlock3DSolver : IMultiphysicsSolver`
  - [x] Build `Grid3D` from builder geometry
  - [x] Forward Euler: `T(n+1) = T(n) + dt * (α * Laplacian3D(T) + source)`
  - [x] Apply 6-face Dirichlet BCs each step
  - [x] Store timeline snapshots at configurable interval
  - [x] ~150–200 LOC

### Tests

- [x] `HeatBlock3D_UniformIC_BoundaryDriven` — symmetric cooling from all faces
- [x] `HeatBlock3D_PointSource_DiffusesRadially` — source at centre, check spherical symmetry
- [x] `HeatBlock3D_SteadyState_MatchesAnalytical` — 1D analytical for Dirichlet T=100 top, T=0 bottom (linear gradient in z)
- [x] `HeatBlock3D_SliceXY_MatchesExpected` — slice accessor returns correct plane
- [x] `HeatBlock3D_Timeline_GrowsWithSteps`

---

## Phase 3 — 3D Diffusion + Advection Fluid (Multiphysics)

Scalar transport in a 3D velocity field: advection-diffusion equation. Simpler than full 3D Navier-Stokes but physically meaningful for contaminant/thermal plume scenarios.

### Equation

∂c/∂t = D ∇²c − v·∇c + S

where D is diffusion coefficient, v = (vx, vy, vz) is a prescribed (or solved) velocity field.

### Implementation

- [x] Add `FluidDiffusion3D` to `MultiphysicsType` enum
- [x] Extend `SimulationBuilder`:
  - [x] `.WithVelocityField3D(Func<double,double,double,(double,double,double)>)` — prescribed velocity (x,y,z) → (vx,vy,vz)
  - [x] `.WithDiffusionCoefficient(double D)` — scalar diffusion
- [x] Create `FluidDiffusion3DSolver : IMultiphysicsSolver`
  - [x] Forward Euler time stepping
  - [x] Each step: `c(n+1) = c(n) + dt * (D * Laplacian3D(c) - Advection3D(c, vx, vy, vz) + source)`
  - [x] CFL check: `max(|vx|,|vy|,|vz|) * dt / min(dx,dy,dz) < 1`
  - [x] 6-face boundary conditions
  - [x] ~150–200 LOC

### Tests

- [x] `FluidDiffusion3D_PureDiffusion_SphericalSpread` — zero velocity, point source, Gaussian-like spread
- [x] `FluidDiffusion3D_UniformAdvection_ShiftsPlume` — uniform v, check plume translates
- [x] `FluidDiffusion3D_MassConservation` — total concentration conserved (Neumann BCs)
- [x] `FluidDiffusion3D_CFL_Warning` — large velocity triggers CFL flag

---

## Phase 4 — 3D Cylinder Flow (Multiphysics, stretch goal)

3D incompressible Navier-Stokes around a cylinder. Extension of the existing 2D `CylinderFlowSolver` with axial (z) variation.

### Equation

∂v/∂t + (v·∇)v = −(1/ρ)∇p + ν∇²v  
∇·v = 0

### Implementation

- [x] Add `CylinderFlow3D` to `MultiphysicsType` enum
- [x] Create `CylinderFlow3DSolver : IMultiphysicsSolver`
  - [x] 3D Chorin projection on `Grid3D`
  - [x] Cylinder mask: `bool[]` for cells inside cylinder (extended along z-axis)
  - [x] Predict: `v* = v + dt * (-Advection3D(v,v) + ν * Laplacian3D(v))`
  - [x] Pressure Poisson: `SolvePoisson3D(∇·v* / dt)`
  - [x] Correct: `v(n+1) = v* - dt * Gradient3D(p)`
  - [x] BCs: inlet (Dirichlet), outlet (Neumann), top/bottom (free-slip or periodic), cylinder (no-slip)
  - [x] Periodic or free-slip in z for initial version
  - [x] Drag/lift coefficients + Strouhal number
  - [x] ~400–500 LOC (most complex solver)
- [x] Extend `SimulationResult`:
  - [x] `Vx3D, Vy3D, Vz3D` — 3D velocity components
  - [x] `Pressure3D` — 3D pressure field

### Tests

- [x] `CylinderFlow3D_LowRe_SymmetricSteadyWake` — Re < 40, check symmetric flow
- [x] `CylinderFlow3D_Drag_ReasonableRange` — Cd ≈ 1.0–1.5 at moderate Re
- [x] `CylinderFlow3D_AxialSymmetry_WithPeriodicZ` — check z-slices are similar for periodic z BCs
- [x] `CylinderFlow3D_VelocityMagnitude_PositiveEverywhere`

### Validation Targets

| Quantity | Expected | Source |
|----------|----------|--------|
| Steady symmetric wake | Re < 47 | Same as 2D |
| Drag coefficient | Cd ≈ 1.0–1.5 at Re ≈ 100 | Literature |
| 3D instability onset | Re ≈ 190 (Mode A) | Williamson (1996) |

---

## Phase 4.5 — 2D Water Contamination (Intermediate Step)

Stepping stone from the existing 1D river-network `WaterContaminationSimulator` to the full 3D volumetric contamination in Phase 5. Implements **2D advection-diffusion** on the full `GeoGrid` for open water bodies (lakes, estuaries, coastal areas).

### Equation

∂c/∂t = Dx·∂²c/∂x² + Dy·∂²c/∂y² − vx·∂c/∂x − vy·∂c/∂y − λ·c + S

### Implementation

- [x] `Engines/GIS/Spread/WaterContamination2D/` folder
- [x] `WaterContamination2DParameters` — sources, anisotropic diffusion (Dx, Dy), velocity field, land mask
- [x] `WaterContamination2DSimulator : ISpreadSimulator` — forward Euler, upwind advection, central-difference diffusion, exponential decay, Neumann BCs, land mask barriers, CFL check
- [x] `WaterContamination2DResult` — snapshots, `MaxConcentration`, `PeakArrivalTimeSeconds`, `AffectedCellCount`, `AffectedAreaM2`, `ExceedanceDurationSeconds()`, `GenerateContaminationExtent()`
- [x] `WaterContamination2DScenarioBuilder` — fluent API: `WithSource`, `WithContaminant`, `WithDiffusion`, `WithCurrentField`, `WithCurrent`, `WithLandMask`, `OverGrid`, `OverTime`, `RunSingle`
- [x] Wire into `RiskScenario.ForWaterContamination2D()`

### Tests (22 tests)

- [x] `PureDiffusion_SpreadsRadially`
- [x] `PureDiffusion_IsSymmetric`
- [x] `Advection_ShiftsPlumeDownstream`
- [x] `Decay_ReducesConcentration`
- [x] `LandMask_BlocksSpread`
- [x] `CflViolation_DetectedForHighVelocity`
- [x] `CflViolation_DetectedForHighDiffusion`
- [x] `NoCflViolation_ForSmallTimeStep`
- [x] `Source_InjectsFixedConcentration`
- [x] `Source_StopsAfterDuration`
- [x] `Result_MaxConcentration_IsPositive`
- [x] `Result_AffectedAreaM2_Grows`
- [x] `Result_ExceedanceDuration_IsPositive`
- [x] `Result_PeakArrivalTime_IsNonNegative`
- [x] `Result_GenerateContaminationExtent_ReturnsPolygon`
- [x] `ScenarioBuilder_RunSingle_ProducesResult`
- [x] `ScenarioBuilder_WithCurrent_AppliesAdvection`
- [x] `ScenarioBuilder_WithLandMask_Blocks`
- [x] `AnisotropicDiffusion_SpreadsMoreInXThanY`
- [x] `MultipleSources_BothContribute`
- [x] `Snapshots_ContainExpectedLayers`
- [x] `ConservativeTracer_MassConserved`

---

## Phase 5 — 3D Volumetric Contamination (GIS Engine)

"How does a contaminant spread in a lake or ocean volume over time?" and "Which depths are affected?"

Extends the GIS engine with a **3D advection-diffusion simulator** on `GeoGrid` (which already supports Nx × Ny × Nz). Unlike the existing `WaterContaminationSimulator` (1D along river), this models volumetric spread in open water bodies.

### Equation

∂c/∂t = Dh(∂²c/∂x² + ∂²c/∂y²) + Dv·∂²c/∂z² − (vx·∂c/∂x + vy·∂c/∂y + vz·∂c/∂z) − λc + S

where:
- Dh = horizontal diffusivity (m²/s), Dv = vertical diffusivity (typically Dv ≪ Dh in stratified water)
- v = (vx, vy, vz) = current velocity field (prescribed or depth-varying)
- λ = first-order decay rate (s⁻¹)
- S = source injection (kg/m³/s)

### Implementation

- [x] Create `Engines/GIS/Spread/VolumetricContamination/` folder
- [x] `ContaminationCellState3D` enum: `Clean = 0`, `Contaminated = 1`, `Source = 2`, `Land = 3`
- [x] `VolumetricContaminationParameters`
  - [x] `Sources` — list of `(ix, iy, iz, concentrationMgL, durationSeconds)`
  - [x] `HorizontalDiffusivity` (m²/s, default ~1.0 for lakes)
  - [x] `VerticalDiffusivity` (m²/s, default ~0.01, much smaller due to stratification)
  - [x] `CurrentVelocity` — `Func<int,int,int,(double,double,double)>` or uniform `(vx, vy, vz)`
  - [x] `DecayRate` (s⁻¹, default 0 = conservative tracer)
  - [x] `LandMask` — `bool[]` optional array marking land cells (non-water)
- [x] `VolumetricContaminationSimulator : ISpreadSimulator`
  - [x] Uses `GeoGrid` (already 3D) directly — no `Grid3D` needed since GeoGrid has `Index(ix,iy,iz)` and uniform `Step`
  - [x] But internally uses `GridOperators3D` math with grid dimensions mapped from `GeoGrid.Nx/Ny/Nz` and `GeoGrid.Step`
  - [x] Forward Euler time stepping
  - [x] Each step:
    - Horizontal diffusion: Dh × (∂²c/∂x² + ∂²c/∂y²)
    - Vertical diffusion: Dv × ∂²c/∂z²
    - Advection: upwind scheme using current velocity
    - Decay: c *= exp(−λ·dt)
    - Source injection: add concentration at source cells
  - [x] Land mask enforcement: zero concentration in land cells
  - [x] Neumann (zero-flux) at domain boundaries
  - [x] CFL check on velocity
  - [x] Store `SpreadSnapshot` per time step with layers: `"concentration"`, `"contaminationState"`, `"exposureTime"`
- [x] `VolumetricContaminationResult`
  - [x] `Snapshots` — time series of `SpreadSnapshot`
  - [x] `Grid` — the `GeoGrid`
  - [x] `PeakConcentration3D()` → `double[]` — per-cell peak over time (flat array, Nx×Ny×Nz)
  - [x] `DepthProfile(ix, iy)` → `DepthProfile` — concentration vs depth at a horizontal point
  - [x] `SurfaceSlice(timeIndex)` → `double[]` — iz=Nz-1 (surface) concentration at time t
  - [x] `BottomSlice(timeIndex)` → `double[]` — iz=0 (bottom) concentration at time t
  - [x] `HorizontalSlice(timeIndex, iz)` → `double[]` — arbitrary depth slice
  - [x] `MaxAffectedDepth(threshold)` → `double` — deepest z where peak concentration > threshold
- [x] `DepthProfile` data class
  - [x] `Depths` — `double[]` (z values, metres)
  - [x] `Concentrations` — `double[]` (mg/L at each depth)
  - [x] `MaxConcentrationDepth` → depth at which concentration is highest

### Fluent Builder

- [x] `VolumetricContaminationScenarioBuilder` — fluent API following same pattern as `WildfireScenarioBuilder`
  - [x] `.OverGrid(GeoGrid)` — must have Nz > 1
  - [x] `.OverTime(start, end, step)`
  - [x] `.WithDiffusivity(horizontal, vertical)`
  - [x] `.WithCurrent(vx, vy, vz)` or `.WithCurrentField(Func<...>)`
  - [x] `.WithSource(ix, iy, iz, concentrationMgL, durationS)`
  - [x] `.WithDecay(rateSec)`
  - [x] `.WithLandMask(bool[])` or `.WithLandMaskFromTerrain(TerrainGrid, seaLevel)`
  - [x] `.RunSingle()` → `VolumetricContaminationResult`
  - [ ] `.RunMonteCarlo(iterations, seed)` — vary current + diffusivity

### GeoJSON Export

- [x] `GeoJsonExporter.VolumetricContaminationToGeoJson(result, iz)` — export a horizontal slice
- [x] `GeoJsonExporter.DepthProfileToGeoJson(result, ix, iy)` — export depth profile at point
- [x] Include `"depth"`, `"isLand"` properties per feature

### Tests

- [x] `VolumetricContamination_PureDiffusion_SphericalSpread` — no current, point source, check radial symmetry
- [x] `VolumetricContamination_VerticalStratification_LimitedDepthSpread` — Dv ≪ Dh, plume stays near source depth
- [x] `VolumetricContamination_HorizontalCurrent_TranslatesPlume` — uniform current shifts plume
- [x] `VolumetricContamination_Decay_ReducesConcentration` — λ > 0 reduces peak over time
- [x] `VolumetricContamination_LandMask_BlocksSpread` — contamination does not enter land cells
- [x] `VolumetricContamination_DepthProfile_ReturnsCorrectShape` — depth profile at source has peak
- [x] `VolumetricContamination_SurfaceSlice_MatchesTopLayer` — surface slice consistent with Nz-1
- [x] `VolumetricContamination_MaxAffectedDepth_IncreaseOverTime` — depth penetration grows
- [x] `VolumetricContamination_MassConservation_NeumannBCs` — total mass conserved (no decay case)
- [x] `VolumetricContamination_CFL_Warning` — large velocity triggers warning

---

## Phase 6 — Integration, Export & Documentation

### Integration

- [x] Wire `VolumetricContaminationScenarioBuilder` into `RiskScenario.ForVolumetricContamination()`
- [x] Ensure `ExposurePolygonGenerator` works with 3D snapshots (pick surface slice by default)
- [ ] Monte Carlo clustering pipeline compatibility for volumetric results

### Export

- [x] `CesiumExporter.ToVolumetricContaminationCzml()` — 3D point cloud or layer-by-layer
- [x] `UnityBinaryExporter.SaveVolumetric()` — binary format for 3D data
- [x] Depth-profile CSV export

### Documentation

- [x] Update `Numerics/Numerics/Numerics/README.md` — document Grid3D + GridOperators3D
- [x] Update `Numerics/Numerics/Engines/README.md` — document 3D Multiphysics solvers
- [x] Update `Numerics/Numerics/Engines/GIS/README.md` — document volumetric contamination with examples

---

## Complexity Estimates

| Component | Approx. LOC | Depends on |
|-----------|-------------|------------|
| `Grid3D` | ~80 | — |
| `GridOperators3D` | ~250 | Grid3D |
| `HeatBlock3DSolver` | ~180 | GridOperators3D |
| `FluidDiffusion3DSolver` | ~180 | GridOperators3D |
| `CylinderFlow3DSolver` | ~450 | GridOperators3D |
| `VolumetricContaminationSimulator` | ~350 | GeoGrid + GridOperators3D |
| `VolumetricContaminationResult` | ~150 | GeoGrid |
| Builder + parameters | ~200 | — |
| Export extensions | ~200 | — |
| Tests (all phases) | ~600 | — |
| **Total** | **~2,640** | |

## Phase Order

```
Phase 1 (Grid3D + Operators) ──┬── Phase 2 (3D Heat)
                                ├── Phase 3 (3D Diffusion-Advection)
                                ├── Phase 4 (3D Cylinder Flow) ← stretch goal
                                └── Phase 5 (Volumetric Contamination GIS)
                                        │
                                        └── Phase 6 (Integration + Export + Docs)
```

Phase 1 is the critical path. Phases 2–5 can be done in any order once Phase 1 is complete, though Phase 2 (heat) is the simplest validation of the 3D FD primitives.
