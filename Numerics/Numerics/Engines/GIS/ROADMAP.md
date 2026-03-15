# GIS / Geo Engine — Roadmap

## Vision

A **fluent simulation engine** for modelling atmospheric gas dispersion over geographic areas.  
Pipeline: **Physics model → Monte Carlo (N scenarios) → ML clustering → probability map → export (GeoJSON / Unity)**.

```
Fluent API (implemented):
─────────────────────────
var result = RiskScenario
    .ForGaussianPlume(emissionRate: 5.0)
    .FromSource(new Vector(0, 0, 50))
    .WithWind(speed: 10, direction: new Vector(1, 0, 0))
    .WithStability(StabilityClass.D)
    .WithVariation(v => v.WindSpeed(8, 12).WindDirectionJitter(15))
    .OverGrid(new GeoGrid(-500, 500, -500, 500, 0, 100, step: 10))
    .OverTime(startSeconds: 0, endSeconds: 3600, stepSeconds: 60)
    .RunMonteCarlo(iterations: 1000, seed: 42)
    .AnalyzeWith(new KMeans { Seed = 42 }, new SilhouetteEvaluator(), 2, 6)
    .Build(threshold: 1e-6);
```

---

## Architecture

```
Engines/
├── Common/                          ← shared (exists)
│   ├── ISimulationEngine.cs         ← re-used
│   ├── SimulationClock.cs           ← re-used
│   ├── FieldSerializer.cs           ← extended
│   └── Extensions/FileExtensions.cs ← extended
│
├── GIS/                             ← NEW — the Geo Engine
│   ├── Scenario/
│   │   ├── RiskScenario.cs               — fluent entry point
│   │   ├── RiskScenarioBuilder.cs        — builder (wind, source, stability, …)
│   │   ├── ScenarioResult.cs             — holds MC results + clustering + maps
│   │   └── TimeFrame.cs                  — time-step config (start, end, step)
│   │
│   ├── Grid/
│   │   ├── GeoGrid.cs                    — 3-D spatial grid (x/y/z bounds + step)
│   │   ├── GeoCell.cs                    — single cell: position, concentration, probability
│   │   └── GridSnapshot.cs               — field values at one time step
│   │
│   ├── Simulation/
│   │   ├── PlumeSimulator.cs             — runs one physics scenario on a grid
│   │   ├── PlumeMonteCarloModel.cs       — IMonteCarloModel adapter (varies wind, stability)
│   │   └── ScenarioVariation.cs          — parameter variation ranges for MC
│   │
│   ├── Analysis/
│   │   ├── ScenarioClusterAnalyzer.cs    — wraps ClusteringGrid for scenario vectors
│   │   ├── ProbabilityMap.cs             — cell → P(concentration > threshold)
│   │   └── TimeAnimator.cs              — assembles time-stepped probability maps
│   │
│   └── Export/
│       ├── GeoJsonExporter.cs            — GeoJSON FeatureCollection output
│       ├── UnityBinaryExporter.cs        — compact binary grid for Unity3D
│       └── CesiumExporter.cs             — GeoJSON + Cesium-compatible metadata
│
Physics/
│   EnvironmentalExtensions.cs       ← exists — Gaussian plume, diffusion, advection
│   (possible additions to EnvironmentalExtensions or new files)
│
Statistics/
│   MonteCarlo/                      ← exists
│   Distributions/                   ← extend with wind/stability distributions
│
ML/
│   Clustering/                      ← exists — KMeans, DBSCAN, Silhouette, …
│   (ClusteringGrid already supports the fluent pattern)
```

---

## Phases

### Phase 0 — Foundation & Grid (`Engines/GIS/Grid/`) ✅

> Goal: a reusable 3-D spatial grid that can hold scalar values per cell per time step.

| # | Item | Location | New / Extend | Depends on | Status |
|---|------|----------|--------------|------------|--------|
| 0.1 | **`GeoGrid`** — 3-D uniform grid (xMin/xMax, yMin/yMax, zMin/zMax, step). Enumerates cell centres. Index ↔ (ix, iy, iz) mapping. | `Engines/GIS/Grid/GeoGrid.cs` | New | — | ✅ |
| 0.2 | **`GeoCell`** — lightweight struct: `Vector Position`, `double Value`, `int TimeIndex`. | `Engines/GIS/Grid/GeoCell.cs` | New | 0.1 | ✅ |
| 0.3 | **`GridSnapshot`** — stores all cell values for one time step; provides lookup by index and by nearest cell. | `Engines/GIS/Grid/GridSnapshot.cs` | New | 0.1, 0.2 | ✅ |
| 0.4 | **`TimeFrame`** — value-object (start, end, stepSeconds); generates `double[] Times`. | `Engines/GIS/Scenario/TimeFrame.cs` | New | — | ✅ |

**Deliverable:** `GeoGrid` can be constructed with the syntax `new GeoGrid(-500,500,-500,500,0,100,10)` and iterates cells. ✅ *25 tests passing.*

---

### Phase 1 — Single-scenario physics simulation (`Engines/GIS/Simulation/`) ✅

> Goal: evaluate the Gaussian plume (from `EnvironmentalExtensions`) on a `GeoGrid` for every time step.

| # | Item | Location | New / Extend | Depends on | Status |
|---|------|----------|--------------|------------|--------|
| 1.1 | **`PlumeSimulator`** — takes emission rate, wind, stability, grid, time frame. For each time step produces a `GridSnapshot` by evaluating `GaussianPlume` or `AdvectionDiffusionRate` at each cell. Supports `SteadyState` and `Transient` (puff) modes. | `Engines/GIS/Simulation/PlumeSimulator.cs` | New | Phase 0, `EnvironmentalExtensions` | ✅ |
| 1.2 | **Time-dependent plume** — added `GaussianPuff()` extension to `EnvironmentalExtensions`: transient puff model with axial + lateral + vertical Gaussian terms, ground reflection, advection at wind speed. | `Physics/EnvironmentalExtensions.cs` | Extend | — | ✅ |
| 1.3 | **`ScenarioVariation`** — fluent value-object defining parameter ranges: wind speed (min/max), wind direction jitter (std-dev degrees), emission rate (min/max), stability class categorical weights. | `Engines/GIS/Simulation/ScenarioVariation.cs` | New | — | ✅ |

**Deliverable:** produce a `List<GridSnapshot>` (one per time step) for a single deterministic scenario. ✅ *14 tests passing.*

---

### Phase 2 — Monte Carlo scenario generation (`Engines/GIS/Simulation/` + `Statistics/`) ✅

> Goal: generate N stochastic scenarios by varying physics parameters; collect per-cell distributions.

| # | Item | Location | New / Extend | Depends on | Status |
|---|------|----------|--------------|------------|--------|
| 2.1 | **`PlumeMonteCarloModel`** — implements `IMonteCarloModel`. `Evaluate` returns peak concentration scalar. `RunBatch(iterations, seed)` returns full `MonteCarloScenarioResult` with scenario matrix + snapshots. Samples wind speed (uniform), emission rate (uniform), wind direction jitter (Gaussian), stability class (categorical weights). | `Engines/GIS/Simulation/PlumeMonteCarloModel.cs` | New | Phase 1, `MonteCarloSimulator` | ✅ |
| 2.2 | **Wind / stability distributions** — sampling built into `PlumeMonteCarloModel`: uniform for wind/emission, Gaussian rotation for direction jitter, categorical `SampleStability()` for stability class. | In `PlumeMonteCarloModel.cs` | In-place | — | ✅ |
| 2.3 | **Batch runner** — `RunBatch(iterations, seed)` collects `ScenarioMatrix` (Matrix of shape N × features), `Snapshots` (List<List<GridSnapshot>>), and `MonteCarloScenarioResult` with `GetCellDistribution` / `GetScenarioVector` helpers. Compatible with `MonteCarloSimulator.Run()`. | `Engines/GIS/Simulation/PlumeMonteCarloModel.cs` | New | 2.1 | ✅ |

**Deliverable:** a `Matrix` of shape (N, cells×timesteps) where each row is one scenario's concentration field. ✅ *10 tests passing.*

---

### Phase 3 — ML clustering & analysis (`Engines/GIS/Analysis/`) ✅

> Goal: cluster the N scenarios to find representative outcome groups; identify the most-probable scenario.

| # | Item | Location | New / Extend | Depends on | Status |
|---|------|----------|--------------|------------|--------|
| 3.1 | **`ScenarioClusterAnalyzer`** — fluent wrapper around `ClusteringExperiment`. Accepts `MonteCarloScenarioResult`, returns `ClusterAnalysisResult` with labels, dominant cluster, `GetClusterMeanSnapshots()`, `GetClusterIterations()`. | `Engines/GIS/Analysis/ScenarioClusterAnalyzer.cs` | New | Phase 2, `ClusteringExperiment` | ✅ |
| 3.2 | **`ProbabilityMap`** — per-cell P(concentration > threshold). Supports cluster-filtered iteration subsets via `iterationFilter`. `CellsAbove(p)`, `At(position)`, 3-D indexer. | `Engines/GIS/Analysis/ProbabilityMap.cs` | New | Phase 2 | ✅ |
| 3.3 | **`TimeAnimator`** — time-stepped sequence of `ProbabilityMap` objects. `ProbabilityAt(position, time)`, `CumulativeProbabilityAt(position, time)`, `ToArray()`. | `Engines/GIS/Analysis/TimeAnimator.cs` | New | 3.2 | ✅ |
| 3.4 | (Optional) **PCA reduction** — before clustering, reduce high-dimensional scenario vectors with existing `PCA` to improve K-Means quality. Already supported via `ClusteringGrid.AddReducer<PCA>(...)`. | Use existing `ML/DimensionalityReduction/` | Re-use | — | ⏳ |

**Deliverable:** `ProbabilityMap[]` array (one per time step) with per-cell risk probabilities. ✅ *19 tests passing.*

---

### Phase 4 — Fluent API & `RiskScenario` (`Engines/GIS/Scenario/`) ✅

> Goal: wire everything behind a clean fluent builder.

| # | Item | Location | New / Extend | Depends on | Status |
|---|------|----------|--------------|------------|--------|
| 4.1 | **`RiskScenario`** — static entry point: `RiskScenario.ForGaussianPlume(Q)`. | `Engines/GIS/Scenario/RiskScenario.cs` | New | — | ✅ |
| 4.2 | **`RiskScenarioBuilder`** — fluent builder collecting `.FromSource()`, `.WithWind()`, `.WithStability()`, `.WithMode()`, `.WithVariation()`, `.OverGrid()`, `.OverTime()`. Terminal `.RunMonteCarlo(N)` → `MonteCarloStageResult`. `.AnalyzeWith(algo, eval)` → `AnalysisStageResult`. `.Build(threshold)` → `ScenarioResult`. Also `.RunSingle()` for deterministic runs. | `Engines/GIS/Scenario/RiskScenarioBuilder.cs` | New | Phases 1-3 | ✅ |
| 4.3 | **`ScenarioResult`** — terminal result holding snapshots, MC data, clustering analysis, `TimeAnimator`, probability maps. Exposes `.ProbabilityAt()`, `.CumulativeProbabilityAt()`, `.ProbabilityMapAt()`. | `Engines/GIS/Scenario/ScenarioResult.cs` | New | Phases 1-3 | ✅ |

**Implemented fluent chain:**

```csharp
var result = RiskScenario
    .ForGaussianPlume(5.0)                           // entry point
    .FromSource(new Vector(0, 0, 50))                // source position + stack height
    .WithWind(10, new Vector(1, 0, 0))               // speed + direction
    .WithStability(StabilityClass.D)                  // Pasquill–Gifford class
    .WithVariation(v => v                             // MC parameter ranges
        .WindSpeed(8, 12)                             //   uniform sampling
        .WindDirectionJitter(15)                      //   degrees std-dev (Gaussian)
        .SetStabilityWeights(d: 0.6, c: 0.2, e: 0.2))  // categorical weights
    .OverGrid(new GeoGrid(-500, 500, -500, 500, 0, 100, 10))
    .OverTime(0, 3600, 60)
    .RunMonteCarlo(1000, seed: 42)                    // → MonteCarloStageResult
    .AnalyzeWith(                                     // → AnalysisStageResult
        new KMeans { Seed = 42 },
        new SilhouetteEvaluator(), minK: 2, maxK: 6)
    .Build(threshold: 1e-6);                          // → ScenarioResult

// Interrogate — probability maps filtered to dominant cluster
double pAt30min = result.ProbabilityAt(
    new Vector(200, 50, 0), timeSeconds: 1800);
double cumulative = result.CumulativeProbabilityAt(
    new Vector(200, 50, 0), timeSeconds: 1800);

// Access underlying data
var clusterLabels = result.ClusterAnalysis.Labels;
var dominantSnaps = result.ClusterAnalysis.GetClusterMeanSnapshots(
    result.ClusterAnalysis.DominantCluster);
var map = result.ProbabilityMapAt(timeIndex: 30);

// Deterministic single-scenario run (no MC)
var single = RiskScenario
    .ForGaussianPlume(5.0)
    .FromSource(new Vector(0, 0, 50))
    .WithWind(10, new Vector(1, 0, 0))
    .WithMode(PlumeMode.Transient, releaseSeconds: 30)
    .OverGrid(new GeoGrid(-500, 500, -500, 500, 0, 100, 10))
    .OverTime(0, 3600, 60)
    .RunSingle();

// Export
result.ExportGeoJson("output/plume.geojson");
result.ExportUnity("output/plume.bin");
result.ExportCesium("output/plume.czml");
```

**Pipeline stages:**

```
RiskScenario.ForGaussianPlume(Q)         → RiskScenarioBuilder
  .FromSource / .WithWind / ...          → RiskScenarioBuilder (config)
  .RunMonteCarlo(N, seed)                → MonteCarloStageResult
    .AnalyzeWith(algo, eval, minK, maxK) → AnalysisStageResult
      .Build(threshold)                  → ScenarioResult
    .Build(threshold)                    → ScenarioResult (skip clustering)
  .RunSingle()                           → ScenarioResult (deterministic)
```

---

### Phase 5 — Export (`Engines/GIS/Export/`) ✅

> Goal: serialize probability maps & time animations to formats consumable by external tools.

| # | Item | Location | New / Extend | Depends on | Status |
|---|------|----------|--------------|------------|--------|
| 5.1 | **`GeoJsonExporter`** — exports `GridSnapshot`, `ProbabilityMap`, combined snapshot+probability, and `TimeAnimator` to GeoJSON FeatureCollections. `ToGeoJson()` string + `Save()` file. `SavePerTimeStep()` for one file per time step. Supports `ExportMetadata`. | `Engines/GIS/Export/GeoJsonExporter.cs` | New | Phase 3 | ✅ |
| 5.2 | **`UnityBinaryExporter`** — compact binary format (GPLM magic, header with grid dims + time range, float[] concentration + probability layers). `Save()` + `Read()` round-trip. Three overloads: concentration-only, probability-only, both. | `Engines/GIS/Export/UnityBinaryExporter.cs` | New | Phase 3 | ✅ |
| 5.3 | **`CesiumExporter`** — CZML output with time-dynamic rgba colour per cell (probability → colour ramp). Also `ToGeoJsonCesium()` with Cesium-specific metadata (3D Tiles hints, `cesium:color`). `ProbabilityToColor()` maps [0..1] to RGBA. | `Engines/GIS/Export/CesiumExporter.cs` | New | 5.1 | ✅ |
| 5.4 | **`ScenarioResult` export methods** — `.ExportGeoJson()`, `.ExportUnity()`, `.ExportCesium()` convenience methods wired into the fluent terminal result. Works for both MC and deterministic runs. | `Engines/GIS/Scenario/ScenarioResult.cs` | Extend | 5.1 | ✅ |

### GeoJSON output structure (per time step)

```json
{
  "type": "FeatureCollection",
  "metadata": {
    "simulation": "GaussianPlume",
    "emissionRate": 5.0,
    "timeStep": 1800,
    "unit": "kg/m³",
    "crs": "local"
  },
  "features": [
    {
      "type": "Feature",
      "geometry": {
        "type": "Point",
        "coordinates": [200, 50, 0]
      },
      "properties": {
        "concentration": 1.23e-5,
        "probability": 0.87,
        "cluster": 1,
        "exceedance_p": 0.72
      }
    }
  ]
}
```

### Unity binary format

```
Header:
  4 bytes  magic "GPLM"
  4 bytes  version (1)
  4 bytes  nx, ny, nz
  4 bytes  timeStepCount
  8 bytes  xMin, xMax, yMin, yMax, zMin, zMax (6 × double)
  8 bytes  tStart, tEnd, tStep (3 × double)

Body (per time step):
  float[nx * ny * nz]  concentration
  float[nx * ny * nz]  probability
```

---

### Phase 6 — GIS coordinate support ✅

> Goal: allow real-world coordinates (lat/lon/altitude) and project to/from the local Cartesian grid.

| # | Item | Location | New / Extend | Status |
|---|------|----------|--------------|--------|
| 6.1 | **`GeoCoordinate`** — (latitude, longitude, altitude) readonly struct with Haversine distance, equality, validation. | `Engines/GIS/Coordinates/GeoCoordinate.cs` | New | ✅ |
| 6.2 | **`Projection`** — Local Tangent Plane (ENU) and UTM projection. WGS-84 ellipsoid. Forward + inverse. Auto-detects UTM zone. | `Engines/GIS/Coordinates/Projection.cs` | New | ✅ |
| 6.3 | **`GeoGrid.FromLatLon()`** — factory method: lat/lon bounding box → auto-projected local-metre grid with attached `Projection`. `CellCentreGeo()` inverse lookup. | `Engines/GIS/Grid/GeoGrid.cs` | Extend | ✅ |
| 6.4 | **GeoJSON/Cesium with CRS** — writes WGS-84 [lon,lat,alt] when grid has projection; `"crs":"WGS84"` metadata. Both `GeoJsonExporter` and `CesiumExporter` updated. | `Engines/GIS/Export/` | Extend | ✅ |

**Tests:** 28 tests in `NumericTest/CoordinateTests.cs` — GeoCoordinate, Projection (LTP + UTM round-trips), GeoGrid.FromLatLon, GeoJSON CRS output.

---

## Where new code lives — decision matrix

| Concept | Location | Rationale |
|---------|----------|-----------|
| Scenario builder, RiskScenario | `Engines/GIS/Scenario/` | Domain-specific orchestration |
| GeoGrid, GeoCell | `Engines/GIS/Grid/` | GIS-specific spatial grid |
| PlumeSimulator, MC adapter | `Engines/GIS/Simulation/` | Geo simulation logic |
| ScenarioClusterAnalyzer, ProbabilityMap, TimeAnimator | `Engines/GIS/Analysis/` | Geo-specific analysis |
| GeoJsonExporter, UnityBinaryExporter, CesiumExporter | `Engines/GIS/Export/` | Geo-specific I/O |
| Gaussian plume physics (existing + transient) | `Physics/EnvironmentalExtensions.cs` | Core physics, not engine-specific |
| New distributions (Weibull wind, wrapped-normal) | `Statistics/Distributions/` | Reusable statistical building blocks |
| Clustering (KMeans, Silhouette etc.) | `ML/Clustering/` (existing) | Already implemented — re-use |
| `MonteCarloSimulator`, `IMonteCarloModel` | `Statistics/MonteCarlo/` (existing) | Already implemented — re-use |
| FieldSerializer extensions | `Engines/Common/` | Shared across engines |

---

## Implementation order (suggested)

```
Phase 0  ████████████████████  Grid foundation              ✅ Done
Phase 1  ████████████████████  Single-scenario physics       ✅ Done
Phase 2  ████████████████████  Monte Carlo generation        ✅ Done
Phase 3  ████████████████████  Clustering & probability maps ✅ Done
Phase 4  ████████████████████  Fluent API (RiskScenario)     ✅ Done
Phase 5  ████████████████████  Export (GeoJSON / Unity / Cesium) ✅ Done
Phase 6  ████████████████████  GIS coordinates (WGS84 / UTM)  ✅ Done
```

Each phase is independently testable. Phase 4 (fluent API) can be started in parallel with Phase 3 by wiring the builder steps that don't yet depend on clustering.

---

## Existing assets to leverage

| Asset | Where | Re-use in Geo Engine |
|-------|-------|---------------------|
| `EnvironmentalExtensions.GaussianPlume()` | `Physics/` | Core physics model in `PlumeSimulator` |
| `EnvironmentalExtensions.AdvectionDiffusionRate()` | `Physics/` | Time-stepping transport |
| `EnvironmentalExtensions.DiffusionPointSource()` | `Physics/` | Transient analytical solution |
| `MonteCarloSimulator` + `IMonteCarloModel` | `Statistics/MonteCarlo/` | Scenario generation |
| `MonteCarloResult` | `Statistics/MonteCarlo/` | Per-cell statistics |
| `ClusteringExperiment` + `ClusteringGrid` | `ML/Clustering/` | Scenario classification |
| `KMeans`, `DBSCAN`, `AgglomerativeClustering` | `ML/Clustering/Algorithms/` | Clustering algorithms |
| `SilhouetteEvaluator` etc. | `ML/Clustering/Evaluators/` | Cluster quality metrics |
| `PCA` | `ML/DimensionalityReduction/` | Optional pre-clustering reduction |
| `ScalarField`, `VectorField`, `Vector` | `Numerics/Objects/` | Field evaluation |
| `FieldSerializer` | `Engines/Common/` | Binary/CSV serialization basis |
| `ISimulationEngine`, `SimulationClock` | `Engines/Common/` | Engine interface |
| `StabilityClass` enum | `Physics/Enums/` | Pasquill–Gifford classes |
