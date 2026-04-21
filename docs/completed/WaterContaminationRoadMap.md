# Water Contamination Spread — Waterway Transport Roadmap

## Overview

A **waterway contamination spread simulation** inside `Engines/GIS/`, modelling point-source pollutant transport in river and stream networks. The physics are based on the **1-D advection–diffusion equation** with optional first-order decay, adsorption (retardation), and tributary mixing — the standard used by EPA QUAL2E, HEC-RAS WQ, and MIKE 11.

The framework reuses the existing GIS engine infrastructure: `GeoGrid` for the spatial domain, `GridSnapshot` + named layers for per-cell state (concentration, velocity, exposure), the Monte Carlo / clustering pipeline for stochastic discharge ensembles, `ExposurePolygonGenerator` for contamination-zone extraction, and the export pipeline (GeoJSON, Cesium, Unity) for visualization.

### Architecture

```
Physics/
├── Materials/
│   ├── Chemical/              ← existing (ChemicalSubstance, ChemicalLibrary)
│   ├── Nuclear/               ← existing (Isotope, IsotopeLibrary, Decay, …)
│   └── Water/                 ← NEW — aquatic contaminant descriptors
│       ├── AquaticContaminant.cs    Name, half-life, Kd, toxicity thresholds
│       ├── ContaminantLibrary.cs    Static registry (Cs-137, Benzene, E. coli, …)
│       └── Enums/
│           └── ContaminantType.cs   Enum: Radioactive, Chemical, Biological, Thermal
├── Environmental/
│   ├── DiffusionExtensions.cs            ← existing (Fick's laws)
│   ├── TransportExtensions.cs            ← existing (advection–diffusion)
│   └── Water/                            ← NEW — river hydraulics & dispersion
│       ├── ManningEquation.cs                 Open-channel flow (discharge, velocity)
│       ├── LongitudinalDispersion.cs          Fischer et al. (1979) dispersion coefficient
│       └── MixingZoneModel.cs                 Tributary dilution / mixing length


Engines/GIS/
├── Grid/                  ← existing (GeoGrid, GridSnapshot, GeoCell)
├── Terrain/               ← existing (TerrainGrid, FuelMap)
│   ├── RiverNetwork.cs    ← NEW — directed flow graph over GeoGrid cells
│   └── ChannelMap.cs      ← NEW — per-cell channel properties (depth, width, Manning n)
├── Spread/
│   ├── ISpreadSimulator.cs           ← existing interface
│   ├── SpreadSnapshot.cs             ← existing (extend with contamination layers)
│   ├── Wildfire/                     ← existing
│   └── WaterContamination/           ← NEW — contaminant transport simulator
│       ├── WaterContaminationSimulator.cs   Implements ISpreadSimulator
│       ├── WaterContaminationParameters.cs  Runtime config (sources, discharge, decay)
│       ├── WaterContaminationResult.cs      Contamination-specific result
│       ├── WaterContaminationMonteCarloResult.cs
│       ├── WaterContaminationVariation.cs   Stochastic ranges
│       └── Enums/
│           └── CellContaminationState.cs    Clean, Contaminated, Decayed, Source
├── Scenario/
│   ├── RiskScenario.cs                      + ForWaterContamination() entry
│   └── WaterContaminationScenarioBuilder.cs ← NEW — fluent configuration
├── Analysis/              ← reuse ExposurePolygonGenerator for contamination extent
└── Export/                ← extend GeoJSON/Cesium with contamination features
```

### Cross-section Dependencies

```
Physics.Environmental.Water.ManningEquation       →  (no internal deps, pure math)
Physics.Environmental.Water.LongitudinalDispersion →  (no internal deps, pure math)
Physics.Materials.Water.AquaticContaminant         →  (no internal deps, data)
Physics.Materials.Water.ContaminantLibrary         →  AquaticContaminant

Engines.GIS.WaterContaminationSimulator  →  Physics.Environmental.Water.*   (hydraulics)
                                         →  Physics.Environmental.TransportExtensions (advection–diffusion)
                                         →  Physics.Environmental.DiffusionExtensions  (Fick's laws)
                                         →  Physics.Materials.Water.AquaticContaminant (contaminant data)
                                         →  Engines.GIS.RiverNetwork                  (flow connectivity)
                                         →  Engines.GIS.ChannelMap                    (per-cell hydraulics)
                                         →  Engines.GIS.GeoGrid                       (spatial domain)
                                         →  Engines.GIS.GridSnapshot                  (layers)
```

### Governing Equation

1-D advection–diffusion with first-order decay and retardation:

$$R_f \frac{\partial C}{\partial t} = E_L \frac{\partial^2 C}{\partial x^2} - u \frac{\partial C}{\partial x} - \lambda C + S$$

| Symbol | Meaning |
|--------|---------|
| $C$ | Contaminant concentration (mg/L or Bq/m³) |
| $u$ | Cross-sectional mean velocity (m/s) from Manning's equation |
| $E_L$ | Longitudinal dispersion coefficient (m²/s) — Fischer (1979) |
| $\lambda$ | First-order decay constant (1/s): $\lambda = \ln 2 / t_{1/2}$ |
| $R_f$ | Retardation factor: $R_f = 1 + \rho_s K_d / n$ (adsorption) |
| $S$ | Source term (mass injection rate per unit volume) |

### Manning's Equation (Open-Channel Velocity)

$$u = \frac{1}{n} R_h^{2/3} S_0^{1/2}$$

| Symbol | Meaning |
|--------|---------|
| $n$ | Manning's roughness coefficient |
| $R_h$ | Hydraulic radius (m) = A / P |
| $S_0$ | Channel bed slope |

### Fischer Dispersion Coefficient

$$E_L = 0.011 \frac{u^2 W^2}{H \cdot u_*}$$

| Symbol | Meaning |
|--------|---------|
| $W$ | Channel width (m) |
| $H$ | Mean depth (m) |
| $u_*$ | Shear velocity = $\sqrt{g R_h S_0}$ |

---

## Phase 1 — River Hydraulics & Contaminant Data (Physics section)

Pure math and data — no grid or simulation. Two namespaces following existing conventions:
- **Materials:** `CSharpNumerics.Physics.Materials.Water` (alongside Chemical, Nuclear, Fire)
- **Physics:** `CSharpNumerics.Physics.Environmental.Water` (alongside Fire)

### Contaminant models — `Physics/Materials/Water/`

- [x] `ContaminantType` enum: `Radioactive`, `Chemical`, `Biological`, `Thermal`
- [x] `AquaticContaminant` immutable record: `string Name`, `ContaminantType Type`, `double HalfLifeSeconds` (0 = conservative tracer), `double PartitionCoefficient` Kd (L/kg, 0 = no adsorption), `double ToxicityThresholdMgL` (drinking water limit), `double LethalThresholdMgL`
- [x] `ContaminantLibrary` static class: `Get(string name)`, `TryGet(name, out c)`, `All`, `Register(AquaticContaminant)` — pre-loaded with:
  - Radioactive: Cs-137, Sr-90, I-131
  - Chemical: Benzene, Toluene, Cyanide, Mercury, Lead, Arsenic
  - Biological: E. coli, Enterococcus
  - Thermal: GenericHeat (no decay, no Kd, threshold in °C delta)

### River hydraulics — `Physics/Environmental/Water/`

- [x] `ManningEquation` static class:
  - `Velocity(double manningN, double hydraulicRadius, double bedSlope)` → `double` (m/s)
  - `Discharge(double manningN, double hydraulicRadius, double bedSlope, double area)` → `double` (m³/s)
  - `HydraulicRadius(double area, double wettedPerimeter)` → `double` (m)
  - `RectangularHydraulicRadius(double width, double depth)` → `double` shorthand
  - `TrapezoidalHydraulicRadius(double bottomWidth, double depth, double sideSlope)` → `double`
- [x] `LongitudinalDispersion` static class:
  - `FischerCoefficient(double velocity, double width, double depth, double shearVelocity)` → `double` EL (m²/s)
  - `ShearVelocity(double hydraulicRadius, double bedSlope)` → `double` u* (m/s)
  - `DecayConstant(double halfLifeSeconds)` → `double` λ (1/s)
  - `RetardationFactor(double bulkDensity, double Kd, double porosity)` → `double` Rf
- [x] `MixingZoneModel` static class:
  - `TributaryMixing(double mainConc, double mainQ, double tribConc, double tribQ)` → `double` downstream concentration (mass balance)
  - `MixingLength(double velocity, double width, double transverseDispersion)` → `double` (m) distance to complete mixing

### Tests — Phase 1

- [x] Manning velocity: rectangular channel (W=10 m, H=2 m, n=0.035, S=0.001) → verify against handbook values
- [x] Manning discharge = velocity × area
- [x] Fischer EL: wide shallow river → higher EL than narrow deep
- [x] Decay constant: Cs-137 half-life ≈ 30.17 yr → λ ≈ 7.28 × 10⁻¹⁰ s⁻¹
- [x] Retardation factor Rf = 1 when Kd = 0 (conservative tracer)
- [x] Tributary mixing: equal flows, one clean + one at 100 mg/L → 50 mg/L downstream
- [x] ContaminantLibrary round-trip: Get/Register/All
- [x] All pre-loaded contaminants have valid thresholds

---

## Phase 2 — River Network & Channel Map (GIS engine)

Build the directed flow graph and per-cell channel properties. These are grid-aware wrappers in `Engines/GIS/Terrain/` that consume `Physics.Environmental.Water`.

### River network — `Engines/GIS/Terrain/RiverNetwork.cs`

- [x] `RiverNetwork` class over a `GeoGrid` (Nz=1):
  - Internal representation: per-cell `bool IsRiver` flag + per-cell `List<(int ix, int iy)> DownstreamNeighbours` (directed graph)
  - `FromElevation(TerrainGrid terrain, double flowThreshold)` — D8 flow-direction algorithm: each cell drains to its steepest downhill neighbour, accumulate flow, mark cells as river where accumulation ≥ threshold
  - `FromManual(GeoGrid grid)` — builder for hand-drawn river paths: `AddSegment(List<(int,int)> cells)`, `SetConfluence(ix, iy, List<upstream>)`, `Build()`
  - `IsRiverCell(ix, iy)` → `bool`
  - `GetDownstream(ix, iy)` → `IReadOnlyList<(int ix, int iy)>` (usually 1, >1 at bifurcation)
  - `GetUpstream(ix, iy)` → `IReadOnlyList<(int ix, int iy)>` (>1 at confluence)
  - `GetReachCells()` → topologically sorted list of all river cells (upstream → downstream)
  - `RiverCellCount` → `int`

### Channel map — `Engines/GIS/Terrain/ChannelMap.cs`

- [x] `ChannelMap` class: assigns hydraulic properties per river cell:
  - `SetChannel(ix, iy, double width, double depth, double manningN)`
  - `SetUniformChannel(double width, double depth, double manningN)` — all river cells
  - `SetChannelByStreamOrder(RiverNetwork net, TerrainGrid terrain, ...)` — wider/deeper downstream
  - `GetWidth(ix, iy)` → `double` (m)
  - `GetDepth(ix, iy)` → `double` (m)
  - `GetManningN(ix, iy)` → `double`
  - `GetBedSlope(ix, iy, RiverNetwork net, TerrainGrid terrain)` → `double` — elevation drop per reach length
  - `GetVelocity(ix, iy, RiverNetwork net, TerrainGrid terrain)` → `double` — Manning's equation using local properties

### Tests — Phase 2

- [x] `RiverNetwork.FromElevation` on a V-shaped valley → single stream line at valley bottom
- [x] Confluence: two upstream branches merge → downstream cell has 2 upstream neighbours
- [x] `FromManual` builder creates specified topology
- [x] `GetReachCells()` returns topologically sorted (upstream before downstream)
- [x] `ChannelMap` set/get round-trip
- [x] `GetVelocity` matches `ManningEquation.Velocity` with same input parameters
- [x] `SetChannelByStreamOrder` produces wider channels downstream

---

## Phase 3 — Advection–Diffusion Contaminant Transport Simulator

Wire `Physics.Environmental.Water` and `Physics.Environmental.TransportExtensions` to the `GeoGrid` via an Eulerian finite-difference scheme that transports contaminant along the river network.

### Spread engine

- [x] `CellContaminationState` enum: `Clean`, `Contaminated`, `Decayed`, `Source`
- [x] `WaterContaminationParameters`:
  - `Sources` — list of `(int ix, int iy, double concentrationMgL, double durationSeconds)` (point-source injections)
  - `Contaminant` — `AquaticContaminant` descriptor (decay, adsorption, thresholds)
  - `BaseDischargeM3s` — reference river discharge (m³/s)
  - `BedPorosity` — for retardation factor (default 0.4)
  - `BedBulkDensity` — kg/m³ sediment (default 1600)
- [x] `WaterContaminationSimulator : ISpreadSimulator`
  - **Initialization**: set concentration = 0 everywhere, mark source cells
  - **Time-step loop** (upstream → downstream topological order):
    1. **Source injection**: at active source cells, add mass: $\Delta C = S \cdot \Delta t / R_f$
    2. **Advection**: forward-in-space finite-difference along flow path: $C_i^{n+1} = C_i^n - u \cdot \Delta t / \Delta x \cdot (C_i^n - C_{i-1}^n)$ (upwind scheme)
    3. **Dispersion**: explicit central difference: $\Delta C = E_L \cdot \Delta t / \Delta x^2 \cdot (C_{i+1} - 2C_i + C_{i-1})$
    4. **Decay**: $C_i^{n+1} = C_i^{n+1} \cdot e^{-\lambda \Delta t}$ (first-order)
    5. **Confluence mixing**: mass-balance at junction cells using `MixingZoneModel.TributaryMixing`
    6. **Record**: produce `SpreadSnapshot` with layers
  - **CFL stability check**: warn if $u \cdot \Delta t / \Delta x > 1$ (Courant condition)
  - **Output layers**:
    - `"concentration"` (mg/L)
    - `"contaminationState"` (Clean=0, Contaminated=1, Decayed=2, Source=3)
    - `"velocity"` (m/s per cell)
    - `"exposureTime"` (seconds above toxicity threshold)

### SpreadSnapshot extensions

- [x] Extend `SpreadSnapshot` (or subclass) with contamination-specific convenience properties:
  - `ContaminatedCellCount` — cells where concentration > toxicity threshold
  - `MaxConcentration` — peak across all cells
  - `AffectedReachLengthKm` — sum of river-cell spacing for contaminated cells
  - `GenerateContaminationExtent(double threshold)` → `ExposurePolygon`

### Tests — Phase 3

- [x] Conservative tracer (no decay, no adsorption): total mass conserved across all time steps
- [x] Single point source in straight uniform channel → concentration plume advects downstream at Manning velocity
- [x] Dispersion: plume spreads (standard deviation increases with time) — verify $\sigma \propto \sqrt{2 E_L t}$
- [x] First-order decay: peak concentration drops exponentially with time
- [x] Retardation: Rf > 1 slows plume front velocity by factor 1/Rf
- [x] Confluence mixing: two branches at known concentrations → correct downstream value
- [x] CFL warning when Δt too large for grid spacing
- [x] `ContaminatedCellCount` increases then decreases as plume passes
- [x] `Clean` cells have concentration = 0, `Source` cells maintain injection concentration

---

## Phase 4 — Fluent API & Scenario Integration

Expose water contamination through the same fluent builder pattern as plume and wildfire, including Monte Carlo for discharge and weather uncertainty.

### Builder

- [x] `RiskScenario.ForWaterContamination()` → returns `WaterContaminationScenarioBuilder`
- [x] `WaterContaminationScenarioBuilder` fluent chain:
  ```csharp
  RiskScenario
      .ForWaterContamination()
      .WithRiverNetwork(riverNetwork)
      .WithChannels(channelMap)
      .WithTerrain(terrainGrid)
      .WithSource(ix, iy, concentrationMgL, durationSeconds)
      .WithContaminant(ContaminantLibrary.Get("Benzene"))
      .WithDischarge(baseDischargeM3s)
      .OverGrid(grid)
      .OverTime(0, 86400, 300)             // 24 hours, 5-min steps
      .RunSingle();                         // → WaterContaminationResult
  ```
- [x] `WaterContaminationResult`:
  - `Snapshots` — `List<SpreadSnapshot>`
  - `MaxConcentration` (mg/L)
  - `PeakArrivalTimeSeconds` — time when max concentration reaches grid outlet
  - `TotalAffectedReachKm`
  - `ExceedanceDurationSeconds(double threshold)` — total time any cell exceeds threshold
  - `GenerateContaminationExtent(timeIndex, threshold)` → `ExposurePolygon`

### Monte Carlo

- [x] `WaterContaminationVariation` — stochastic parameter ranges:
  - Discharge range (low/high flow)
  - Source concentration uncertainty
  - Manning's n uncertainty range
- [x] `.WithVariation(v => v.Discharge(5, 50).SourceConcentration(80, 120))`
- [x] `.RunMonteCarlo(iterations)` → `WaterContaminationMonteCarloResult`
  - Per-cell exceedance probability across all iterations
  - Mean / P95 peak concentration downstream
  - Worst-case plume arrival time
- [x] `.AnalyzeWith(clustering)` — reuse `ScenarioClusterAnalyzer` on exceedance probability matrix
- [x] `.Build(threshold)` → `WaterContaminationResult` with probability-weighted outputs

### Tests — Phase 4

- [x] Fluent API deterministic round-trip: build → run → inspect affected reach
- [x] Monte Carlo 20 iterations: exceedance probability ∈ [0, 1] for all cells
- [x] Higher discharge → faster plume arrival, lower peak concentration (dilution)
- [x] Monte Carlo with variation produces valid iterations and statistics
- [x] `GenerateContaminationExtent()` returns valid polygon enclosing contaminated cells

---

## Phase 5 — Export & Visualization

Extend the existing export pipeline with contamination-specific outputs.

### GeoJSON

- [x] Point features with `concentration`, `velocity`, `contaminationState` properties
- [x] Contamination extent as `Polygon` geometry per time step (via `ExposurePolygonGenerator`)
- [x] Exceedance probability heatmap export (from Monte Carlo)
- [x] River centreline as `LineString` with time-varying concentration attribute

### Cesium (CZML)

- [x] Time-dynamic contamination extent polygons (animated transport)
- [x] Colour ramp: blue (clean) → yellow (low concentration) → red (above toxicity) → black (source)

### Unity binary

- [x] Extend `UnityBinaryExporter` to write contamination layers (concentration, velocity)

### Tests — Phase 5

- [x] GeoJSON export produces valid FeatureCollection with contamination polygons
- [x] CZML export has time intervals matching simulation time steps
- [x] Binary round-trip: write → read → verify layer values

---

## Phase 6 — Documentation & Polish

- [x] Update `Physics/README.md` with Water Hydraulics section:
  - ManningEquation overview and formula
  - LongitudinalDispersion (Fischer) reference
  - AquaticContaminant / ContaminantLibrary / ContaminantType reference
  - Standalone usage examples (velocity, dispersion calculation without grid)
- [x] Update `Engines/GIS/README.md` with Water Contamination section:
  - RiverNetwork, ChannelMap overview
  - WaterContaminationSimulator usage (consuming Physics.Environmental.Water)
  - Fluent API examples
  - Monte Carlo discharge uncertainty workflow
  - Export examples
- [x] Add code samples for common scenarios:
  - Industrial spill in straight river
  - Tributary dilution at confluence
  - Radioactive release (Cs-137) with decay and adsorption
  - MC ensemble with discharge uncertainty
- [x] Validate all public types have XML doc summaries
- [x] Final test pass — all water contamination tests green

---

## Future Extensions (Out of Scope for MVP)

These are **not** part of the MVP but inform the architecture to keep extensibility open:

- **2-D lateral dispersion** — transverse mixing in wide rivers or estuaries
- **Sediment transport** — suspended load, bed load, erosion/deposition coupling
- **Bioaccumulation** — contaminant transfer through aquatic food chain
- **Stormwater runoff** — non-point-source loading from urban/agricultural areas
- **Lake/reservoir stratification** — vertical mixing in still-water bodies
- **Estuary tidal mixing** — bi-directional flow with salinity wedge
- **Real-time data ingestion** — USGS/SMHI gauge data feeding discharge time series
- **DEM-based network extraction** — automated watershed delineation from GeoTIFF
- **Water quality standards** — automated compliance checking against EU WFD / US CWA limits
- **Multi-species reactions** — coupled reaction networks (e.g., nitrification, BOD–DO)
- **RL suppression** — optimal boom placement / treatment timing via reinforcement learning
