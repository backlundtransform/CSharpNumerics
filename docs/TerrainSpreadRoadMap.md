# Terrain Spread Framework — Wildfire (Rothermel) MVP Roadmap

## Overview

A **terrain-aware spread simulation framework** inside `Engines/GIS/`, starting with **wildfire** as the first scenario. The fire physics are based on the **Rothermel (1972) surface fire spread model** — the standard used by FARSITE, FlamMap, and BehavePlus.

The framework reuses the existing GIS engine infrastructure: `GeoGrid` for the spatial domain, `GridSnapshot` + named layers for per-cell state (fuel, moisture, flame length, burn time), the Monte Carlo / clustering pipeline for stochastic weather ensembles, `ExposurePolygonGenerator` for fire-perimeter extraction, and the export pipeline (GeoJSON, Cesium, Unity) for visualization.

### Architecture

```
Physics/
├── Materials/
│   ├── Chemical/              ← existing (ChemicalSubstance, ChemicalLibrary)
│   ├── Nuclear/               ← existing (Isotope, IsotopeLibrary, Decay, …)
│   ├── Engineering/           ← existing (EngineeringMaterial, EngineeringLibrary)
│   └── Fire/                  ← NEW — fuel model data (same pattern as Chemical/Nuclear)
│       ├── FuelModel.cs            Rothermel fuel parameters (Anderson 13)
│       ├── FuelLibrary.cs          Static registry of standard fuel models
│       └── Enums/
│           └── FuelModelType.cs    Enum: ShortGrass, TimberLitter, …
├── Environmental/
│   └── Fire/                  ← NEW — fire spread physics (pure math, no grid)
│       └── RothermelModel.cs       Core Rothermel equations (R, IR, φw, φs)
├── EnvironmentalExtensions.cs ← existing (GaussianPlume, GaussianPuff)

Engines/GIS/
├── Grid/                  ← existing (GeoGrid, GridSnapshot, GeoCell)
├── Terrain/               ← NEW — terrain model + fuel map (grid-aware)
│   ├── TerrainGrid.cs          Elevation surface, slope/aspect
│   └── FuelMap.cs              Per-cell fuel assignment on grid
├── Spread/                ← NEW — generic spread engine
│   ├── ISpreadSimulator.cs     Interface for any terrain spread model
│   ├── SpreadResult.cs         Timeline of SpreadSnapshot (layers per step)
│   ├── SpreadSnapshot.cs       Per-step cell state (burning, burned, unburned, flame length, ROS)
│   └── Wildfire/
│       ├── WildfireSimulator.cs     Cell-automaton spread using Physics.Environmental.Fire.RothermelModel
│       ├── WildfireParameters.cs    Runtime config (ignition, wind, moisture)
│       └── Enums/
│           └── CellBurnState.cs     Unburned, Burning, Burned, Firebreak
├── Scenario/              ← extend existing fluent API
│   ├── RiskScenario.cs          + ForWildfire() entry point
│   ├── WildfireScenarioBuilder.cs  Fluent configuration for fire scenarios
│   └── WildfireScenarioResult.cs   Fire-specific result (perimeters, area burned, …)
├── Simulation/            ← existing plume (unchanged)
├── Analysis/              ← reuse ExposurePolygonGenerator for fire perimeters
├── Export/                ← extend GeoJSON/Cesium with fire-specific features
└── RL/                    ← future: fire suppression RL environment
```

### Cross-section Dependencies

```
Physics.Environmental.Fire.RothermelModel  →  Physics.Materials.Fire.FuelModel  (pure math, no grid deps)
Physics.Materials.Fire.FuelLibrary         →  Physics.Materials.Fire.FuelModel

Engines.GIS.WildfireSimulator  →  Physics.Environmental.Fire.RothermelModel  (fire spread physics)
                               →  Physics.Materials.Fire.FuelModel           (fuel data)
                               →  Engines.GIS.TerrainGrid                   (slope/aspect)
                               →  Engines.GIS.FuelMap                       (fuel per cell, wraps FuelModel)
                               →  Engines.GIS.GeoGrid                      (spatial domain, existing)
                               →  Engines.GIS.GridSnapshot                 (layers, existing)
```

This follows the same separation as the plume pipeline:
- **Materials:** `Physics.Materials.Fire.FuelModel` sits alongside `Physics.Materials.Chemical` and `Physics.Materials.Nuclear`
- **Physics:** `Physics.Environmental.Fire.RothermelModel` extends the environmental physics alongside `EnvironmentalExtensions.GaussianPlume()`
- **Engine:** `Engines.GIS.WildfireSimulator` wires everything to the grid, just like `Engines.GIS.PlumeSimulator`

### Rothermel Model Summary

Rate of spread (m/min):

$$R = \frac{I_R \cdot \xi \cdot (1 + \phi_w + \phi_s)}{\rho_b \cdot \varepsilon \cdot Q_{ig}}$$

| Symbol | Meaning |
|--------|---------|
| $I_R$ | Reaction intensity (kJ/m²·min) |
| $\xi$ | Propagating flux ratio |
| $\phi_w$ | Wind correction factor |
| $\phi_s$ | Slope correction factor |
| $\rho_b$ | Ovendry bulk density (kg/m³) |
| $\varepsilon$ | Effective heating number |
| $Q_{ig}$ | Heat of pre-ignition (kJ/kg) |

Fuel parameters per model: surface-area-to-volume ratio σ, fuel bed depth δ, ovendry fuel load $w_0$, dead fuel moisture of extinction $M_x$, low heat content $h$.

---

## Phase 1 — Fire Physics & Fuel Library (Physics section)

Pure math and data — no grid or simulation. Two namespaces following existing conventions:
- **Materials:** `CSharpNumerics.Physics.Materials.Fire` (alongside Chemical, Nuclear, Engineering)
- **Physics:** `CSharpNumerics.Physics.Environmental.Fire` (alongside GaussianPlume)

### Fuel models (Anderson 13) — `Physics/Materials/Fire/`

- [ ] `FuelModel` immutable record: `FuelModelType Type`, `string Name`, `double SurfaceAreaToVolumeRatio` (1/m), `double FuelBedDepth` (m), `double OvendryFuelLoad` (kg/m²), `double MoistureOfExtinction` (fraction), `double LowHeatContent` (kJ/kg), `double ParticleDensity` (kg/m³ — default 513 for wood)
- [ ] `FuelModelType` enum for Anderson 13: `ShortGrass (1)`, `TimberGrassUnderstory (2)`, `TallGrass (3)`, `Chaparral (4)`, `Brush (5)`, `DormantBrush (6)`, `SouthernRough (7)`, `ClosedTimberLitter (8)`, `HardwoodLitter (9)`, `TimberLitterUnderstory (10)`, `LightLoggingSlash (11)`, `MediumLoggingSlash (12)`, `HeavyLoggingSlash (13)`, `NoFuel (0)`
- [ ] `FuelLibrary` static class: `Get(FuelModelType)`, `TryGet(type, out fuel)`, `All`, `Register(FuelModel)` — pre-loaded with all 13 Anderson models

### Rothermel core equations — `Physics/Environmental/Fire/RothermelModel.cs`

- [ ] `RothermelModel` static class with:
  - `RateOfSpread(FuelModel fuel, double moistureContent, double windSpeed, double slopeRadians)` → `double` (m/min)
  - `ReactionIntensity(FuelModel fuel, double moistureContent)` → `double` IR (kJ/m²·min)
  - `WindFactor(FuelModel fuel, double midflameWindSpeed)` → `double` φw
  - `SlopeFactor(double packingRatio, double slopeRadians)` → `double` φs
  - `PropagatingFluxRatio(FuelModel fuel)` → `double` ξ
  - `HeatOfPreignition(double moistureContent)` → `double` Qig (kJ/kg)
  - `EffectiveHeatingNumber(FuelModel fuel)` → `double` ε
  - `PackingRatio(FuelModel fuel)` → `double` β = ρb/ρp
  - `OptimalPackingRatio(FuelModel fuel)` → `double` β_op
  - `FlameLength(double reactionIntensity, double rateOfSpread)` → `double` (m) — Byram's fireline intensity → flame length correlation

### Tests — Phase 1

- [ ] `FuelLibrary` returns all 13 Anderson models, verify Short Grass (σ=3500 1/ft ≈ 11483 1/m, δ=1 ft ≈ 0.305 m)
- [ ] Short Grass (model 1), moisture 0.05, wind 5 mph, flat → R ≈ 23–25 m/min (BehavePlus reference)
- [ ] Chaparral (model 4), moisture 0.10, wind 10 mph, flat → verify against BehavePlus
- [ ] Slope factor: flat terrain → φs = 0
- [ ] Slope factor: 30° slope → φs > 0, increasing with slope
- [ ] Wind factor: zero wind → φw = 0
- [ ] Moisture at extinction → R = 0 (no spread)
- [ ] Moisture above extinction → R = 0 (clamp)
- [ ] Flame length proportional to fireline intensity

---

## Phase 2 — Terrain Model & Fuel Map (GIS engine)

Build the elevation surface and per-cell fuel assignment. These are grid-aware wrappers that live in `Engines/GIS/Terrain/` and consume `Physics.Fire.FuelModel`.

### Terrain grid — `Engines/GIS/Terrain/TerrainGrid.cs`

- [ ] `TerrainGrid` class wrapping a `GeoGrid` (ground plane, Nz=1) with a `double[] Elevation` array (one height per (ix,iy) cell)
- [ ] `FromFunction(grid, Func<double,double,double> elevationFn)` — procedural elevation from f(x,y)
- [ ] `FromArray(grid, double[,] elevation)` — load from 2D array (row = iy, col = ix)
- [ ] `Slope(ix, iy)` — terrain slope in radians using central-difference gradient: $\tan(\theta) = \sqrt{(\partial z/\partial x)^2 + (\partial z/\partial y)^2}$
- [ ] `Aspect(ix, iy)` — downslope direction in radians (0=N, π/2=E, π=S, 3π/2=W)
- [ ] `SlopeInDirection(ix, iy, Vector2 direction)` — slope component along a given heading (needed for directional Rothermel $\phi_s$)

### Fuel map — `Engines/GIS/Terrain/FuelMap.cs`

- [ ] `FuelMap` class: assigns a `Physics.Materials.Fire.FuelModel` per (ix,iy) cell on a `GeoGrid`
  - `SetFuel(ix, iy, FuelModelType)`
  - `SetUniformFuel(FuelModelType)` — fill entire grid
  - `SetFuelByElevation(ranges)` — assign fuel models to elevation bands
  - `GetFuel(ix, iy)` → `FuelModel`
  - `GetMoisture(ix, iy)` → `double` — per-cell dead fuel moisture content (fraction)
  - `SetMoisture(ix, iy, double)` / `SetUniformMoisture(double)`

### Tests — Phase 2

- [ ] `TerrainGrid` slope/aspect on flat surface → slope ≈ 0
- [ ] `TerrainGrid` slope/aspect on known tilted plane → verify against analytical result
- [ ] `SlopeInDirection` matches full slope when direction = aspect, zero when perpendicular
- [ ] `FuelMap` set/get round-trip, uniform fill, moisture defaults

---

## Phase 3 — Cell-Automaton Fire Spread Simulator

Wire `Physics.Environmental.Fire.RothermelModel` to the `GeoGrid` via a cellular automaton that propagates fire across the terrain surface.

### Spread engine

- [ ] `CellBurnState` enum: `Unburned`, `Burning`, `Burned`, `Firebreak`
- [ ] `ISpreadSimulator` interface:
  ```csharp
  IReadOnlyList<SpreadSnapshot> Run(GeoGrid grid, TerrainGrid terrain, FuelMap fuelMap,
      WildfireParameters parameters, TimeFrame timeFrame);
  ```
- [ ] `SpreadSnapshot` class — extends or wraps `GridSnapshot` with:
  - Layer `"burnState"` (0=Unburned, 1=Burning, 2=Burned, 3=Firebreak)
  - Layer `"flameLength"` (metres)
  - Layer `"rateOfSpread"` (m/min)
  - Layer `"burnTime"` (seconds since ignition, 0 if unburned)
  - `BurningCellCount`, `BurnedCellCount`, `BurnedAreaHectares`
- [ ] `WildfireParameters`:
  - `IgnitionPoints` — list of (ix,iy) or world positions
  - `MidflameWindSpeed` (m/s)
  - `WindDirection` (Vector — same convention as PlumeSimulator)
  - `BurnDuration` (seconds — how long a cell stays in Burning state before transitioning to Burned)
  - `SpotFireEnabled` (bool) — future, default false
- [ ] `WildfireSimulator : ISpreadSimulator`
  - 8-neighbour spread (N/S/E/W + diagonals, diagonal distance = step√2)
  - For each Burning cell, compute ROS towards each unburned neighbour:
    1. Direction from burning cell to neighbour
    2. Slope in that direction from `TerrainGrid.SlopeInDirection()`
    3. Wind component along that direction
    4. Neighbour fuel + moisture → `RothermelModel.RateOfSpread()`
    5. Travel time = distance / ROS
    6. If accumulated time ≥ travel time → ignite neighbour
  - Time step loop: advance time, update burning→burned transitions, ignite reachable neighbours
  - Produce one `SpreadSnapshot` per time step

### Tests — Phase 3

- [ ] Single ignition on flat uniform Short Grass with no wind → near-circular spread pattern
- [ ] Ignition with steady wind → elliptical spread (elongated downwind)
- [ ] Uphill slope accelerates spread, downhill decelerates
- [ ] `NoFuel` cells block fire (act as firebreaks)
- [ ] `Firebreak` cell state blocks spread
- [ ] `BurnedAreaHectares` increases monotonically
- [ ] Zero ROS at high moisture (saturated fuel) → fire does not spread

---

## Phase 4 — Fluent API & Scenario Integration

Expose wildfire through the same fluent builder pattern as the plume scenario, including Monte Carlo for weather uncertainty.

### Builder

- [ ] `RiskScenario.ForWildfire()` → returns `WildfireScenarioBuilder`
- [ ] `WildfireScenarioBuilder` fluent chain:
  ```csharp
  RiskScenario
      .ForWildfire()
      .WithTerrain(terrainGrid)
      .WithFuel(fuelMap)
      .WithIgnition(ix, iy)              // or .WithIgnition(position)
      .WithWind(speed, direction)
      .WithMoisture(0.08)                // global dead fuel moisture
      .OverGrid(grid)
      .OverTime(0, 7200, 60)            // 2 hours, 1-min steps
      .RunSingle();                      // → WildfireScenarioResult
  ```
- [ ] `WildfireScenarioResult` — like `ScenarioResult` but fire-specific:
  - `Snapshots` — `List<SpreadSnapshot>`
  - `FinalBurnedArea` (hectares)
  - `MaxFlameLength` (metres)
  - `FirePerimeters` — list of `ExposurePolygon` (one per time step, from burn boundary)
  - `GenerateFirePerimeter(timeIndex)` → `ExposurePolygon` via `ExposurePolygonGenerator` on `burnState ≥ 1`

### Monte Carlo

- [ ] `WildfireVariation` — stochastic parameter ranges:
  - Wind speed range
  - Wind direction jitter
  - Moisture content range
  - Ignition location offset radius
- [ ] `.WithVariation(v => v.WindSpeed(3, 8).Moisture(0.04, 0.12))`
- [ ] `.RunMonteCarlo(iterations)` → `WildfireMonteCarloResult`
  - Per-cell burn probability across all iterations
  - Mean / max burned area statistics
- [ ] `.AnalyzeWith(clustering)` — reuse existing `ScenarioClusterAnalyzer` on the burn-probability matrix
- [ ] `.Build()` → `WildfireScenarioResult` with probability-weighted outputs

### Tests — Phase 4

- [ ] Fluent API deterministic round-trip: build → run → inspect burned area
- [ ] Monte Carlo 50 iterations: burn probability ∈ [0, 1] for all cells
- [ ] Clustering identifies distinct fire spread regimes (high-wind vs. low-wind clusters)
- [ ] `GenerateFirePerimeter()` returns valid polygon enclosing burned cells

---

## Phase 5 — Export & Visualization

Extend the existing export pipeline with fire-specific outputs.

### GeoJSON

- [ ] Point features with `burnState`, `flameLength`, `rateOfSpread` properties
- [ ] Fire perimeter as `Polygon` geometry per time step (via `ExposurePolygonGenerator`)
- [ ] Burn probability heatmap export (from Monte Carlo)

### Cesium (CZML)

- [ ] Time-dynamic fire perimeter polygons (animated spread)
- [ ] Colour ramp: red (burning) → grey (burned) → green (unburned)

### Unity binary

- [ ] Extend `UnityBinaryExporter` to write fire layers (burnState, flameLength)

### Tests — Phase 5

- [ ] GeoJSON export produces valid FeatureCollection with fire perimeter polygons
- [ ] CZML export has time intervals matching simulation time steps
- [ ] Binary round-trip: write → read → verify layer values

---

## Phase 6 — Documentation & Polish

- [ ] Update `Physics/README.md` with Fire section:
  - RothermelModel overview and equations
  - FuelModel / FuelLibrary / FuelModelType reference (Physics.Materials.Fire)
  - Standalone usage examples (ROS calculation without grid)
- [ ] Update `Engines/GIS/README.md` with Wildfire / Terrain Spread section:
  - TerrainGrid, FuelMap overview
  - WildfireSimulator usage (consuming Physics.Environmental.Fire + Physics.Materials.Fire)
  - Fluent API examples
  - Monte Carlo fire probability workflow
  - Export examples
- [ ] Add code samples for common scenarios:
  - Flat grassland fire
  - Mountainous terrain with mixed fuel
  - MC ensemble with wind uncertainty
- [ ] Validate all public types have XML doc summaries
- [ ] Final test pass — all wildfire tests green

---

## Future Extensions (Out of Scope for MVP)

These are **not** part of the MVP but inform the architecture to keep extensibility open:

- **Spot fires** — ember transport via lofting model (wind + convection column)
- **Crown fire** — Van Wagner (1977) crown fire initiation + spread
- **Scott & Burgan 40 fuel models** — expanded fuel library
- **Weather timeline** — time-varying wind speed/direction/moisture during simulation
- **Suppression modelling** — firebreaks, water drops, crew lines (RL environment candidate)
- **Additional spread scenarios** — flood/inundation, landslide, disease/pest, pollutant propagation (all reuse `ISpreadSimulator` + `TerrainGrid`)
- **DEM import** — read GeoTIFF / ASCII grid elevation data
- **Fuel map import** — read LANDFIRE / Corine raster land cover → fuel model mapping
