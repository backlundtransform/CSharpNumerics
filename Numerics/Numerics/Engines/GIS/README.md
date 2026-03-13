## 🌍 GIS / Geo Engine

A fluent simulation engine for modelling atmospheric gas dispersion over geographic areas. Combines Gaussian plume physics, Monte Carlo uncertainty quantification, ML clustering, and probability mapping into a single pipeline — with export to GeoJSON (Cesium) and Unity3D binary.

**Namespace:** `CSharpNumerics.Engines.GIS`

---

### Pipeline Overview

```
Physics model  →  Monte Carlo (N scenarios)  →  ML clustering  →  Probability map  →  Export
```

The engine builds on existing library capabilities:

| Step | Uses |
|------|------|
| Physics | `EnvironmentalExtensions.GaussianPlume` / `GaussianPuff` |
| Monte Carlo | `MonteCarloSimulator`, `IMonteCarloModel` |
| Clustering | `ClusteringGrid`, `KMeans`, `SilhouetteEvaluator` |
| Numerics | `ScalarField`, `VectorField`, `Vector`, `Matrix` |

---

### GeoGrid — 3-D spatial grid

A uniform 3-D grid defined by axis-aligned bounds and a cell step size:

```csharp
using CSharpNumerics.Engines.GIS.Grid;

// 1 km × 1 km area, ground to 100 m altitude, 10 m resolution
var grid = new GeoGrid(
    xMin: -500, xMax: 500,
    yMin: -500, yMax: 500,
    zMin: 0,    zMax: 100,
    step: 10);

int totalCells = grid.CellCount;   // Nx × Ny × Nz
Vector centre = grid.CellCentre(5, 10, 3);  // world-space position

// Index round-trip
int flat = grid.Index(5, 10, 3);
var (ix, iy, iz) = grid.Index3D(flat);

// Find nearest cell to a world position
var (nx, ny, nz) = grid.NearestIndex(new Vector(123, -45, 8));

// Enumerate all cell centres
foreach (Vector pos in grid)
    Console.WriteLine(pos);
```

### GeoCell & GridSnapshot

```csharp
// GridSnapshot holds scalar values for one time step
double[] values = new double[grid.CellCount];
var snapshot = new GridSnapshot(grid, values, time: 60.0, timeIndex: 1);

// Access by index
double val = snapshot[5, 10, 0];            // by (ix, iy, iz)
double val2 = snapshot.ValueAt(new Vector(200, 50, 0));  // by nearest cell

// Query
double max = snapshot.Max();
var hotCells = snapshot.CellsAbove(threshold: 1e-6);

// Iterate all cells with position + value
foreach (GeoCell cell in snapshot.AllCells())
    Console.WriteLine($"{cell.Position} = {cell.Value:E3}");
```

### TimeFrame

```csharp
using CSharpNumerics.Engines.GIS.Scenario;

// 1 hour simulation, 1-minute steps
var tf = new TimeFrame(start: 0, end: 3600, stepSeconds: 60);

double[] times = tf.ToArray();             // [0, 60, 120, ..., 3600]
double t30min = tf.TimeAt(30);             // 1800.0
int idx = tf.NearestIndex(1800);           // 30
```

---

### PlumeSimulator — single-scenario physics

Evaluates the Gaussian plume (steady-state or transient puff) on a `GeoGrid` for every time step:

```csharp
using CSharpNumerics.Engines.GIS.Simulation;
using CSharpNumerics.Physics.Enums;

// Steady-state: same concentration field at every time step
var sim = new PlumeSimulator(
    emissionRate: 5.0,           // kg/s
    windSpeed: 10,               // m/s
    windDirection: new Vector(1, 0, 0),
    stackHeight: 50,             // metres
    sourcePosition: new Vector(0, 0, 50),
    stability: StabilityClass.D,
    mode: PlumeMode.SteadyState);

List<GridSnapshot> snapshots = sim.Run(grid, tf);

// Transient: puff advects downwind over time
var simTransient = new PlumeSimulator(
    emissionRate: 5.0,
    windSpeed: 10,
    windDirection: new Vector(1, 0, 0),
    stackHeight: 50,
    sourcePosition: new Vector(0, 0, 50),
    stability: StabilityClass.D,
    mode: PlumeMode.Transient);
simTransient.ReleaseSeconds = 10;

List<GridSnapshot> transientSnaps = simTransient.Run(grid, tf);

// Single time step
GridSnapshot snap = sim.RunSingle(grid, time: 300);
```

---

### ScenarioVariation — parameter ranges for Monte Carlo

Defines how physics parameters are varied across stochastic iterations:

```csharp
var variation = new ScenarioVariation()
    .WindSpeed(8, 12)                    // uniform [8, 12] m/s
    .WindDirectionJitter(15)             // Gaussian jitter, σ = 15°
    .EmissionRate(3, 7)                  // uniform [3, 7] kg/s
    .SetStabilityWeights(                // categorical sampling
        d: 0.6, c: 0.2, e: 0.2);

bool hasVar = variation.HasVariation;    // true
```

---

### PlumeMonteCarloModel — N stochastic scenarios

Runs N scenarios by sampling from `ScenarioVariation` and collecting per-cell concentration distributions:

```csharp
var mcModel = new PlumeMonteCarloModel(
    emissionRate: 5.0,
    windSpeed: 10,
    windDirection: new Vector(1, 0, 0),
    stackHeight: 50,
    sourcePosition: new Vector(0, 0, 50),
    grid: grid,
    timeFrame: tf,
    variation: variation,
    stability: StabilityClass.D,
    mode: PlumeMode.SteadyState);

// Full batch: returns scenario matrix + all snapshots
MonteCarloScenarioResult result = mcModel.RunBatch(
    iterations: 1000,
    seed: 42);

// Scenario matrix: (1000 × cells·timeSteps)
Matrix scenarioMatrix = result.ScenarioMatrix;

// Per-cell distribution across all scenarios
double[] cellDist = result.GetCellDistribution(
    cellIndex: grid.NearestFlatIndex(new Vector(200, 0, 0)),
    timeIndex: 30);                    // at t = 30 min

// Single scenario vector
double[] row = result.GetScenarioVector(0);
```

**Compatible with the existing `MonteCarloSimulator`** for scalar summary statistics (peak concentration):

```csharp
using CSharpNumerics.Statistics.MonteCarlo;

var sim = new MonteCarloSimulator(seed: 42);
MonteCarloResult peakStats = sim.Run(mcModel, iterations: 1000);

double meanPeak = peakStats.Mean;
double p95 = peakStats.Percentile(95);
```

---

### Fluent API (planned)

The full fluent chain will be available once all phases are complete:

```csharp
var result = RiskScenario
    .ForGaussianPlume(5.0)
    .FromSource(new Vector(0, 0, 50))
    .WithWind(10, new Vector(1, 0, 0))
    .WithStability(StabilityClass.D)
    .WithVariation(v => v
        .WindSpeed(8, 12)
        .WindDirectionJitter(15)
        .SetStabilityWeights(d: 0.6, c: 0.2, e: 0.2))
    .OverGrid(new GeoGrid(-500, 500, -500, 500, 0, 100, 10))
    .OverTime(0, 3600, 60)
    .RunMonteCarlo(1000)
    .AnalyzeWith(new ClusteringGrid()
        .AddModel<KMeans>(g => g.Add("K", 3, 5))
        .WithEvaluator(new SilhouetteEvaluator()))
    .Build();

// Query probability at a point and time
double p = result.ProbabilityAt(
    new Vector(200, 50, 0), timeSeconds: 1800, thresholdKgM3: 1e-6);

// Export
result.ExportGeoJson("output/plume.geojson");
result.ExportUnity("output/plume.bin");
```

---

### Architecture

```
Engines/GIS/
├── Grid/
│   ├── GeoGrid.cs              — 3-D spatial grid with index mapping
│   ├── GeoCell.cs              — position + value + time struct
│   └── GridSnapshot.cs         — cell values at one time step
├── Scenario/
│   ├── TimeFrame.cs            — time range value-object
│   ├── RiskScenario.cs         — fluent entry point (planned)
│   ├── RiskScenarioBuilder.cs  — builder (planned)
│   └── ScenarioResult.cs       — result container (planned)
├── Simulation/
│   ├── PlumeSimulator.cs       — single-scenario physics evaluation
│   ├── PlumeMonteCarloModel.cs — MC batch runner + IMonteCarloModel
│   └── ScenarioVariation.cs    — parameter variation ranges
├── Analysis/                    — (planned)
│   ├── ScenarioClusterAnalyzer.cs
│   ├── ProbabilityMap.cs
│   └── TimeAnimator.cs
└── Export/                      — (planned)
    ├── GeoJsonExporter.cs
    ├── UnityBinaryExporter.cs
    └── CesiumExporter.cs
```

### Status

| Phase | Description | Status |
|-------|-------------|--------|
| 0 | Grid foundation (`GeoGrid`, `GeoCell`, `GridSnapshot`, `TimeFrame`) | ✅ Done |
| 1 | Single-scenario physics (`PlumeSimulator`, `GaussianPuff`) | ✅ Done |
| 2 | Monte Carlo generation (`PlumeMonteCarloModel`, `ScenarioVariation`) | ✅ Done |
| 3 | ML clustering & probability maps | 🔲 Planned |
| 4 | Fluent API (`RiskScenario` builder) | 🔲 Planned |
| 5 | Export (GeoJSON / Unity / Cesium) | 🔲 Planned |
| 6 | Real-world coordinates (lat/lon) | 🔲 Future |

See [ROADMAP.md](ROADMAP.md) for the full design and phase details.
