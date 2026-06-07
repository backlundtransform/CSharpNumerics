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

### ScenarioClusterAnalyzer — ML clustering of scenarios

Finds the most probable emission scenario by clustering the Monte Carlo scenario matrix:

```csharp
using CSharpNumerics.Engines.GIS.Analysis;

var analysis = ScenarioClusterAnalyzer
    .For(mcResult)
    .WithAlgorithm(new ClusteringGrid().AddModel<KMeans>(g => g.Add("K", 3, 5)))
    .WithEvaluator(new SilhouetteEvaluator())
    .Run();

int dominant = analysis.DominantCluster;      // largest cluster
int[] members = analysis.GetClusterIterations(dominant);
List<GridSnapshot> mean = analysis.GetClusterMeanSnapshots(dominant);
```

### ProbabilityMap & TimeAnimator — probability mapping

```csharp
// Exceedance probability per cell at one time step
var pMap = ProbabilityMap.Build(
    mcResult.Snapshots,
    timeIndex: 30,
    threshold: 1e-6,
    iterationFilter: members);   // optional: filter to dominant cluster

double pAt = pMap.At(new Vector(200, 50, 0));
var hotCells = pMap.CellsAbove(0.5);

// Time-animated probability
var animation = TimeAnimator.Build(mcResult.Snapshots, threshold: 1e-6);
double p = animation.ProbabilityAt(new Vector(200, 50, 0), timeSeconds: 1800);
double cumP = animation.CumulativeProbabilityAt(new Vector(200, 50, 0), 1800);
```

---

### Fluent API — RiskScenario

End-to-end pipeline in a single fluent chain:

```csharp
using CSharpNumerics.Engines.GIS.Scenario;

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
    .AnalyzeWith(
        new ClusteringGrid().AddModel<KMeans>(g => g.Add("K", 3, 5)),
        new SilhouetteEvaluator())
    .Build(threshold: 1e-6);

// Query probability at a point and time
double p = result.ProbabilityAt(new Vector(200, 50, 0), timeSeconds: 1800);
double cumP = result.CumulativeProbabilityAt(new Vector(200, 50, 0), 1800);
ProbabilityMap map = result.ProbabilityMapAt(timeIndex: 30);

// Export
result.ExportGeoJson("output/plume.geojson");
result.ExportUnity("output/plume.bin");
result.ExportCesium("output/plume.czml");
```

---

### Export — GeoJSON, Unity binary, Cesium CZML

```csharp
using CSharpNumerics.Engines.GIS.Export;

// GeoJSON — single snapshot or full animation
GeoJsonExporter.Save(snapshot, "snap.geojson");
GeoJsonExporter.Save(animation, "anim.geojson");
var paths = GeoJsonExporter.SavePerTimeStep(animation, "out/step");

// Unity binary — compact float[] format with GPLM header
UnityBinaryExporter.Save(snapshots, animation, "plume.bin");
var data = UnityBinaryExporter.Read("plume.bin");  // round-trip

// Cesium CZML — time-animated colour-coded entities
CesiumExporter.Save(animation, "plume.czml", name: "Plume Risk");
string json = CesiumExporter.ToGeoJsonCesium(probMap);
```

---

### GIS Coordinates — WGS84, UTM, local tangent plane

Create geo-referenced grids from latitude/longitude bounding boxes with automatic projection:

```csharp
using CSharpNumerics.Engines.GIS.Coordinates;

// Geographic coordinate
var stockholm = new GeoCoordinate(59.3293, 18.0686, altitude: 25);
var gothenburg = new GeoCoordinate(57.7089, 11.9746);
double dist = stockholm.DistanceTo(gothenburg);  // ~398 km

// Local Tangent Plane projection
var proj = new Projection(stockholm, ProjectionType.LocalTangentPlane);
Vector local = proj.ToLocal(59.34, 18.08);        // → metres (east, north, up)
GeoCoordinate back = proj.ToGeo(local);            // → lat/lon round-trip

// UTM projection (auto-detects zone)
var utmProj = new Projection(stockholm, ProjectionType.UTM);
int zone = utmProj.UtmZone;  // 34

// Geo-referenced grid from lat/lon bounding box
var sw = new GeoCoordinate(59.32, 18.06);
var ne = new GeoCoordinate(59.34, 18.10);
var geoGrid = GeoGrid.FromLatLon(sw, ne, altMin: 0, altMax: 100, step: 50);

// Cell centres in WGS84
GeoCoordinate geo = geoGrid.CellCentreGeo(0);  // lat, lon, altitude

// GeoJSON export automatically uses [lon, lat, alt] with "crs":"WGS84"
var snapshot = new GridSnapshot(geoGrid, values, 0, 0);
string json = GeoJsonExporter.ToGeoJson(snapshot);  // WGS84 coordinates
```

---


#### Materials & Fluent Pipeline Integration

```csharp
using CSharpNumerics.Physics.Materials;

// Create a material descriptor (auto-attaches known decay chain)
var material = Materials.Radioisotope("Cs137");

// Use in the fluent pipeline — adds "activity" and "dose" layers to snapshots
var result = RiskScenario
    .ForGaussianPlume(5.0)
    .FromSource(new Vector(0, 0, 50))
    .WithWind(10, new Vector(1, 0, 0))
    .WithStability(StabilityClass.D)
    .WithMaterial(Materials.Radioisotope("Cs137"))    // ← attach isotope
    .WithVariation(v => v.WindSpeed(8, 12))
    .OverGrid(new GeoGrid(-500, 500, -500, 500, 0, 100, 10))
    .OverTime(0, 3600, 60)
    .RunMonteCarlo(1000)
    .Build(threshold: 1e-6);

// Query named layers on any snapshot
GridSnapshot snap = result.Snapshots[0];
double[] activityLayer = snap.GetLayer("activity");   // Bq/m³
double[] doseLayer = snap.GetLayer("dose");           // Sv
bool hasIt = snap.HasLayer("activity");               // true

// GeoJSON export automatically includes activity & dose properties
result.ExportGeoJson("fallout.geojson");
```

Output GeoJSON includes per-feature nuclear properties:

```json
{
  "type": "Feature",
  "geometry": { "type": "Point", "coordinates": [200, 50, 0] },
  "properties": {
    "concentration": 0.0012,
    "activity": 120.5,
    "dose": 0.004
  }
}
```

---

### Exposure Polygons — Peak & Time-Integrated Contours

Generate closed polygons delineating areas where an exposure metric (peak or time-integrated)
exceeds a threshold. Works with concentration or any named layer (activity, dose).

```csharp
using CSharpNumerics.Engines.GIS.Analysis;

// Run a scenario with radioactive material
var result = RiskScenario
    .ForGaussianPlume(5.0)
    .FromSource(new Vector(0, 0, 50))
    .WithWind(10, new Vector(1, 0, 0))
    .WithStability(StabilityClass.D)
    .WithMaterial(Materials.Radioisotope("Cs137"))
    .OverGrid(new GeoGrid(-500, 500, -500, 500, 0, 100, 50))
    .OverTime(0, 3600, 60)
    .RunSingle();

// Peak exceedance polygon — where the maximum dose at any time step exceeds threshold
ExposurePolygon peakDose = result.GeneratePeakExposurePolygon(
    threshold: 1e-6, layerName: "dose");

// Time-integrated exposure polygon — where cumulative dose exceeds threshold
ExposurePolygon integratedDose = result.GenerateIntegratedExposurePolygon(
    threshold: 1e-4, layerName: "dose");

// Polygon properties
List<Vector> boundary = peakDose.Boundary;      // closed ring of vertices
int cells = peakDose.ExceedanceCellCount;       // cells above threshold
double area = peakDose.AreaSquareMetres;         // approximate area (m²)
ExposureType type = peakDose.Type;               // Peak or Integrated

// Also works with plain concentration (no material needed)
var concPoly = result.GeneratePeakExposurePolygon(threshold: 1e-6);

// Export to GeoJSON Polygon
string json = GeoJsonExporter.ToGeoJson(peakDose);
GeoJsonExporter.Save(peakDose, "peak_dose_zone.geojson");

// Export multiple polygons as FeatureCollection
string fc = GeoJsonExporter.ToGeoJson(
    new List<ExposurePolygon> { peakDose, integratedDose });
```

Output GeoJSON polygon:

```json
{
  "type": "Feature",
  "geometry": {
    "type": "Polygon",
    "coordinates": [[[100, 50, 0], [150, 75, 0], [200, 50, 0], [100, 50, 0]]]
  },
  "properties": {
    "threshold": 1e-6,
    "layerName": "dose",
    "exposureType": "peak",
    "exceedanceCellCount": 42,
    "areaSquareMetres": 105000
  }
}
```

---

### Terrain & Fuel Maps

The terrain subsystem provides elevation surfaces and per-cell fuel assignment for terrain-aware spread simulations.

**TerrainGrid** — elevation surface with slope and aspect:

```csharp
using CSharpNumerics.Engines.GIS.Terrain;

var grid = new GeoGrid(0, 1000, 0, 1000, 0, 0, 25);

// Procedural terrain from a function
var terrain = TerrainGrid.FromFunction(grid, (x, y) => 100 + 0.1 * x);

// Or from a 2D array
var terrain2 = TerrainGrid.FromArray(grid, elevationArray);

// Slope and aspect at a cell
double slope  = terrain.Slope(10, 10);         // radians
double aspect = terrain.Aspect(10, 10);        // radians (0=N, π/2=E)
double sDir   = terrain.SlopeInDirection(10, 10, windDir);  // slope along heading
```

**FuelMap** — per-cell fuel model and moisture:

```csharp
using CSharpNumerics.Physics.Materials.Fire.Enums;

var fuelMap = new FuelMap(grid);
fuelMap.SetUniformFuel(FuelModelType.ShortGrass);
fuelMap.SetUniformMoisture(0.08);  // 8% dead fuel moisture

// Mixed fuels by elevation band
fuelMap.SetFuelByElevation(terrain, new[]
{
    (0.0,   300.0, FuelModelType.ShortGrass),
    (300.0, 600.0, FuelModelType.Brush),
    (600.0, 999.0, FuelModelType.TimberLitterUnderstory)
});

// Per-cell moisture variation
fuelMap.SetMoisture(5, 5, 0.15);   // wet spot
```

---

### Wildfire Simulator

A cellular automaton that propagates fire across terrain using the Rothermel (1972) rate-of-spread model. Each time step, burning cells attempt to ignite their 8 neighbours based on directional slope, wind, and fuel parameters.

```csharp
using CSharpNumerics.Engines.GIS.Spread.Wildfire;
using CSharpNumerics.Engines.GIS.Scenario;

var parameters = new WildfireParameters(
    ignitionPoints: new List<(int, int)> { (20, 20) }.AsReadOnly(),
    midflameWindSpeed: 3.0,              // m/s
    windDirection: new Vector(1, 0, 0),  // east
    burnDuration: 600);                  // seconds

var simulator = new WildfireSimulator(parameters);
var snapshots = simulator.Run(grid, terrain, fuelMap,
    new TimeFrame(0, 3600, 60));  // 1 hour, 1-min steps

// Inspect results
var last = snapshots[snapshots.Count - 1];
Console.WriteLine($"Burned: {last.BurnedAreaHectares:F1} ha");
Console.WriteLine($"Burning cells: {last.BurningCellCount}");
Console.WriteLine($"Flame length: {last.Snapshot.GetLayer("flameLength").Max()} m");
```

**SpreadSnapshot layers:**

| Layer | Description |
|-------|-------------|
| `burnState` | 0 = Unburned, 1 = Burning, 2 = Burned, 3 = Firebreak |
| `flameLength` | Flame length in metres |
| `rateOfSpread` | Rate of spread in m/min |
| `burnTime` | Seconds since cell ignition |

---

### Wildfire Fluent API

Wildfire scenarios use the same fluent pattern as the plume pipeline, accessed via `RiskScenario.ForWildfire()`:

**Deterministic run:**

```csharp
var result = RiskScenario
    .ForWildfire()
    .WithTerrain(terrain)
    .WithFuel(fuelMap)
    .WithIgnition(20, 20)
    .WithWind(5.0, new Vector(1, 0, 0))
    .WithMoisture(0.08)
    .OverGrid(grid)
    .OverTime(0, 7200, 60)        // 2 hours, 1-min steps
    .RunSingle();

double area = result.FinalBurnedArea;       // hectares
double flame = result.MaxFlameLength;       // metres
var perimeter = result.GenerateFirePerimeter(result.Snapshots.Count - 1);
```

**Monte Carlo ensemble with weather uncertainty:**

```csharp
var mcResult = RiskScenario
    .ForWildfire()
    .WithTerrain(terrain)
    .WithFuel(fuelMap)
    .WithIgnition(20, 20)
    .WithWind(5.0, new Vector(1, 0, 0))
    .WithVariation(v => v
        .WindSpeed(2, 8)              // uniform [2, 8] m/s
        .WindDirectionJitter(30)      // ±30° Gaussian jitter
        .Moisture(0.04, 0.12))        // uniform [4%, 12%]
    .OverGrid(grid)
    .OverTime(0, 3600, 60)
    .RunMonteCarlo(100, seed: 42);

// Per-cell burn probability across all iterations
double[] burnProb = mcResult.BurnProbability;
double meanArea   = mcResult.MeanBurnedArea;   // hectares
double maxArea    = mcResult.MaxBurnedArea;
```

**Clustering + probability map:**

```csharp
using CSharpNumerics.ML.Clustering.Algorithms;
using CSharpNumerics.ML.Clustering.Evaluators;

var probMap = mcResult
    .AnalyzeWith(
        new KMeans { Seed = 42 },
        new SilhouetteEvaluator(),
        minK: 2, maxK: 5)
    .Build();   // → WildfireScenarioResult from dominant cluster

// probMap.Snapshots contain probability-weighted burn state (0–1)
// Use for risk assessment visualizations
```

---

### Wildfire Export

All three export formats support fire-specific data:

**GeoJSON — fire snapshots with burn properties:**

```csharp
string json = GeoJsonExporter.ToGeoJson(snapshots[snapshots.Count - 1]);
// Features with: burnState, flameLength, rateOfSpread, timeStep
```

**GeoJSON — fire perimeters as polygons:**

```csharp
var perimeters = result.FirePerimeters;
string json = GeoJsonExporter.FirePerimetersToGeoJson(perimeters, result.Snapshots);
// Polygon features with timeIndex, timeStep, areaSquareMetres
```

**GeoJSON — burn probability heatmap:**

```csharp
string json = GeoJsonExporter.BurnProbabilityToGeoJson(mcResult);
// Point features with burnProbability, iterations
```

**CZML — animated fire spread:**

```csharp
CesiumExporter.SaveFireCzml(result.Snapshots, timeFrame, "fire.czml");
// Per-cell time-sampled colour: red (burning) → grey (burned)
```

**Unity binary — fire layers:**

```csharp
UnityBinaryExporter.SaveFire(result.Snapshots, grid, timeFrame, "fire.bin");
var data = UnityBinaryExporter.ReadFire("fire.bin");
// data.Concentration = burnState[][], data.Probability = flameLength[][]
```

---

### River Network & Channel Geometry

Build hydrological connectivity from elevation data (D8 flow routing) or define rivers manually with a fluent builder. `ChannelMap` attaches hydraulic geometry (width, depth, Manning's n) to each river cell.

**Automatic extraction from DEM:**

```csharp
using CSharpNumerics.Engines.GIS.Terrain;

var grid = new GeoGrid(0, 1000, 0, 1000, 0, 0, 100);
var terrain = TerrainGrid.FromFunction(grid, (x, y) => 100 - y * 0.001);

// D8 flow direction → river cells above accumulation threshold
var river = RiverNetwork.FromElevation(terrain, accumulationThreshold: 4);
int count = river.RiverCellCount;
(int dr, int dc)? ds = river.GetDownstream(5, 5);
```

**Manual builder:**

```csharp
var river = RiverNetwork.FromManual(grid)
    .AddReach(new[] { (0,5), (1,5), (2,5), (3,5), (4,5) })
    .AddTributary(parentRow: 2, parentCol: 5,
        cells: new[] { (2,3), (2,4), (2,5) })
    .Build();

var sorted = river.TopologicalOrder();   // upstream → downstream (Kahn's algorithm)
var reach  = river.GetReachCells(0, 5);  // all cells downstream from (0,5)
```

**ChannelMap — hydraulic geometry:**

```csharp
var channels = new ChannelMap(grid);
channels.SetUniformChannel(width: 10, depth: 2, manningN: 0.035);

// Vary geometry by Strahler stream order
channels.SetChannelByStreamOrder(river, new Dictionary<int, (double w, double d, double n)>
{
    [1] = (3, 0.5, 0.04),
    [2] = (8, 1.5, 0.035),
    [3] = (15, 3.0, 0.03),
});

double slope = channels.GetBedSlope(terrain, 3, 5);  // from elevation difference
double u = channels.GetVelocity(terrain, 3, 5);       // Manning velocity at cell
```

---

### Water Contamination Simulator

A 1-D advection-diffusion simulator that transports dissolved or adsorbed contaminants along a river network. Implements `ISpreadSimulator` so it integrates directly with the GIS pipeline.

**Governing equation:**

$$R_f \frac{\partial C}{\partial t} = E_L \frac{\partial^2 C}{\partial x^2} - u \frac{\partial C}{\partial x} - \lambda C + S$$

Where $R_f$ = retardation factor, $E_L$ = Fischer dispersion, $u$ = Manning velocity, $\lambda$ = first-order decay, $S$ = source term.

```csharp
using CSharpNumerics.Engines.GIS.Spread.WaterContamination;
using CSharpNumerics.Physics.Materials.Water;

var parameters = new WaterContaminationParameters(
    sources: new List<(int row, int col, double concentrationMgL)>
        { (0, 5, 100.0) }.AsReadOnly(),
    contaminant: AquaticContaminant.Cs137,
    baseDischargeM3s: 20.0,
    bedPorosity: 0.4,
    bedBulkDensity: 1600);

var simulator = new WaterContaminationSimulator(parameters);
var snapshots = simulator.Run(grid, terrain, river, channels,
    new TimeFrame(0, 3600, 60));  // 1 hour, 1-min steps

var last = snapshots[^1];
double peak = last.MaxConcentration();              // mg/L
int affected = last.ContaminatedCellCount(1e-6);    // cells > threshold
double reachKm = last.AffectedReachLengthKm(1e-6); // km of river contaminated
```

**SpreadSnapshot layers:**

| Layer | Description |
|-------|-------------|
| `concentration` | Dissolved concentration (mg/L) |
| `contaminationState` | 0 = Clean, 1 = Contaminated, 2 = Decayed, 3 = Source |
| `velocity` | Manning velocity at cell (m/s) |
| `exposureTime` | Seconds the cell has been above threshold |

---

### Water Contamination Fluent API

Accessed via `RiskScenario.ForWaterContamination()`:

**Deterministic run:**

```csharp
using CSharpNumerics.Engines.GIS.Scenario;

var result = RiskScenario
    .ForWaterContamination()
    .WithRiverNetwork(river)
    .WithChannels(channels)
    .WithTerrain(terrain)
    .WithSource(0, 5, 100.0)                      // row, col, mg/L
    .WithContaminant(AquaticContaminant.Benzene)
    .WithDischarge(20.0)                           // m³/s
    .WithBedProperties(porosity: 0.4, bulkDensity: 1600)
    .OverGrid(grid)
    .OverTime(0, 3600, 60)
    .RunSingle();

double peak = result.MaxConcentration;             // mg/L
double arrival = result.PeakArrivalTimeSeconds;    // seconds to peak
double reachKm = result.TotalAffectedReachKm;      // contaminated river length
var extent = result.GenerateContaminationExtent(1e-4); // polygon of affected area
```

**Monte Carlo ensemble with parameter uncertainty:**

```csharp
var mcResult = RiskScenario
    .ForWaterContamination()
    .WithRiverNetwork(river)
    .WithChannels(channels)
    .WithTerrain(terrain)
    .WithSource(0, 5, 100.0)
    .WithContaminant(AquaticContaminant.Cs137)
    .WithDischarge(20.0)
    .WithBedProperties(0.4, 1600)
    .WithVariation(v => v
        .Discharge(10, 40)                         // uniform [10, 40] m³/s
        .SourceConcentration(50, 200)              // uniform [50, 200] mg/L
        .ManningN(0.025, 0.045))                   // uniform [0.025, 0.045]
    .OverGrid(grid)
    .OverTime(0, 3600, 60)
    .RunMonteCarlo(200, seed: 42);

double meanPeak = mcResult.MeanPeakConcentration;
double p95Peak  = mcResult.P95PeakConcentration;
double[] exceedProb = mcResult.ExceedanceProbability;  // per-cell [0,1]
```

**Clustering + analysis:**

```csharp
using CSharpNumerics.ML.Clustering.Algorithms;
using CSharpNumerics.ML.Clustering.Evaluators;

var analysisResult = mcResult
    .AnalyzeWith(
        new KMeans { Seed = 42 },
        new SilhouetteEvaluator(),
        minK: 2, maxK: 5)
    .Build();   // → WaterContaminationResult from dominant cluster
```

---

### Water Contamination Export

All three export formats support contamination data:

**GeoJSON — contamination snapshots:**

```csharp
string json = GeoJsonExporter.ContaminationToGeoJson(snapshots[^1]);
// Point features with: concentration, contaminationState, velocity, exposureTime

// Exceedance probability heatmap
string heatmap = GeoJsonExporter.ExceedanceProbabilityToGeoJson(
    mcResult.ExceedanceProbability, grid);

// Contamination extent as polygon
string extent = GeoJsonExporter.ContaminationExtentToGeoJson(
    result.GenerateContaminationExtent(1e-4), grid);

// River centreline as LineString
string centreline = GeoJsonExporter.RiverCentrelineToGeoJson(river, grid);
```

**CZML — animated contamination plume:**

```csharp
CesiumExporter.SaveContaminationCzml(snapshots, timeFrame, "contamination.czml");
// Per-cell time-sampled colour: blue → yellow → red → purple (by concentration)
```

**Unity binary — contamination layers:**

```csharp
UnityBinaryExporter.SaveContamination(snapshots, grid, timeFrame, "water.bin");
var data = UnityBinaryExporter.ReadContamination("water.bin");
// data.Concentration[][], data.Velocity[][]
```

---

### 2D Water Contamination Simulator

Extends the 1D river model to 2D open water (lakes, estuaries, coastal areas). Implements full-grid advection-diffusion with anisotropic diffusion and a prescribed velocity field.

```csharp
using CSharpNumerics.Engines.GIS.Scenario;

var result = RiskScenario
    .ForWaterContamination2D()
    .WithSource(50, 50, 100.0, double.MaxValue)
    .WithContaminant(ContaminantLibrary.Get("Benzene"))
    .WithDiffusion(1.0, 0.5)                        // Dx, Dy (anisotropic)
    .WithCurrent(0.1, 0)                             // uniform flow
    .WithLandMask(coastlineMask)
    .OverGrid(grid)
    .OverTime(0, 86400, 60)
    .RunSingle();

double peak = result.MaxConcentration;
double area = result.AffectedAreaM2;
var polygon = result.GenerateContaminationExtent(result.Snapshots.Count - 1, 0.01);
```

---

### Volumetric (3D) Water Contamination

Full 3D advection-diffusion for contaminant transport with depth variation. Models vertical stratification (Dv ≪ Dh), 3D current fields, and depth profiles.

**Governing equation:**

$$\frac{\partial c}{\partial t} = D_h\left(\frac{\partial^2 c}{\partial x^2} + \frac{\partial^2 c}{\partial y^2}\right) + D_v\frac{\partial^2 c}{\partial z^2} - \mathbf{v} \cdot \nabla c - \lambda c + S$$

**Fluent API:**

```csharp
var grid = new GeoGrid(0, 1000, 0, 1000, 0, 50, 10); // 100×100×5 cells

var result = RiskScenario
    .ForVolumetricContamination()
    .WithSource(50, 50, 3, 100.0, double.MaxValue)   // ix, iy, iz, mg/L, duration
    .WithDiffusivity(1.0, 0.01)                        // Dh=1, Dv=0.01 (stratified)
    .WithCurrent(0.1, 0, 0)                            // surface current
    .WithLandMask(seabedMask)
    .OverGrid(grid)
    .OverTime(0, 86400, 60)
    .RunSingle();

// 3D analysis
double peak = result.MaxConcentration;
double[] peakField = result.PeakConcentration3D();     // per-cell peak over time
double maxDepth = result.MaxAffectedDepth(0.01);       // deepest contaminated z

// Depth profile at a point
var profile = result.GetDepthProfile(50, 50);
double peakDepth = profile.MaxConcentrationDepth;

// Horizontal slices
double[] surface = result.SurfaceSlice(result.Snapshots.Count - 1);
double[] bottom  = result.BottomSlice(result.Snapshots.Count - 1);
double[] atDepth = result.HorizontalSlice(result.Snapshots.Count - 1, iz: 2);

// Surface contamination extent polygon
var polygon = result.GenerateContaminationExtent(
    result.Snapshots.Count - 1, 0.01);
```

**Export:**

```csharp
// GeoJSON — horizontal slice at depth layer
string json = GeoJsonExporter.VolumetricContaminationToGeoJson(result, iz: 3);
// Includes: concentration, contaminationState, depth, isLand properties

// GeoJSON — depth profile at a point
string profile = GeoJsonExporter.DepthProfileToGeoJson(result, ix: 50, iy: 50);

// CZML — full 3D time-animated (uses cell altitude for true volumetric vis)
CesiumExporter.SaveVolumetricContaminationCzml(result, timeFrame, "vol.czml");

// Unity binary
UnityBinaryExporter.SaveVolumetric(result, timeFrame, "vol3d.bin");

// CSV depth profile
DepthProfileCsvExporter.Save(result, 50, 50, "depth_profile.csv");
```

**SpreadSnapshot layers (same as 2D):**

| Layer | Description |
|-------|-------------|
| `concentration` | Dissolved concentration (mg/L) |
| `contaminationState` | 0 = Clean, 1 = Contaminated, 2 = Source, 3 = Land |
| `velocity` | Current speed magnitude (m/s) |
| `exposureTime` | Seconds above toxicity threshold |

---

### Architecture

```
Engines/GIS/
├── Coordinates/
│   ├── GeoCoordinate.cs        — WGS-84 lat/lon/alt value-type
│   └── Projection.cs           — LTP & UTM projection (forward + inverse)
├── Grid/
│   ├── GeoGrid.cs              — 3-D spatial grid + FromLatLon() factory
│   ├── GeoCell.cs              — position + value + time struct
│   └── GridSnapshot.cs         — cell values at one time step + named layers
├── Terrain/
│   ├── TerrainGrid.cs          — elevation surface, slope/aspect
│   ├── FuelMap.cs              — per-cell fuel model + moisture
│   ├── RiverNetwork.cs         — D8 flow routing + manual builder + topo sort
│   └── ChannelMap.cs           — per-cell width/depth/Manning's n
├── Spread/
│   ├── ISpreadSimulator.cs     — generic spread interface
│   ├── SpreadSnapshot.cs       — per-step state (+ contamination helpers)
│   ├── Wildfire/
│   │   ├── WildfireSimulator.cs          — 8-neighbour Rothermel CA
│   │   ├── WildfireParameters.cs         — ignition, wind, burn duration
│   │   ├── WildfireScenarioBuilder.cs    — fluent builder
│   │   ├── WildfireScenarioResult.cs     — deterministic result
│   │   ├── WildfireMonteCarloResult.cs   — MC result + AnalyzeWith()
│   │   ├── WildfireAnalysisResult.cs     — clustering + Build()
│   │   ├── WildfireVariation.cs          — stochastic parameter ranges
│   │   └── Enums/
│   │       └── CellBurnState.cs          — Unburned, Burning, Burned, Firebreak
│   └── WaterContamination/
│       ├── WaterContaminationSimulator.cs      — 1-D advection-diffusion along river
│       ├── WaterContaminationParameters.cs     — sources, contaminant, discharge
│       ├── WaterContaminationScenarioBuilder.cs — fluent builder
│       ├── WaterContaminationResult.cs         — deterministic result
│       ├── WaterContaminationMonteCarloResult.cs — MC result + AnalyzeWith()
│       ├── WaterContaminationAnalysisResult.cs — clustering + Build()
│       ├── WaterContaminationVariation.cs      — discharge/conc/Manning ranges
│       └── Enums/
│           └── CellContaminationState.cs       — Clean, Contaminated, Decayed, Source
│   ├── WaterContamination2D/
│   │   ├── WaterContamination2DSimulator.cs    — 2-D advection-diffusion on grid
│   │   ├── WaterContamination2DParameters.cs   — anisotropic diffusion, velocity field
│   │   ├── WaterContamination2DScenarioBuilder.cs — fluent builder
│   │   └── WaterContamination2DResult.cs       — deterministic result
│   └── VolumetricContamination/
│       ├── VolumetricContaminationSimulator.cs  — 3-D advection-diffusion (Nx×Ny×Nz)
│       ├── VolumetricContaminationParameters.cs — Dh/Dv, 3D velocity, land mask, decay
│       ├── VolumetricContaminationScenarioBuilder.cs — fluent builder
│       ├── VolumetricContaminationResult.cs    — depth profiles, slices, 3D peak
│       ├── DepthProfile.cs                     — concentration-vs-depth at (ix,iy)
│       └── Enums/
│           └── ContaminationCellState3D.cs     — Clean, Contaminated, Source, Land
├── Scenario/
│   ├── TimeFrame.cs            — time range value-object
│   ├── RiskScenario.cs         — fluent entry point (+ ForWildfire())
│   ├── RiskScenarioBuilder.cs  — pipeline builder + stage results
│   └── ScenarioResult.cs       — terminal result with export methods
├── Simulation/
│   ├── PlumeSimulator.cs       — single-scenario physics evaluation (+ material)
│   ├── PlumeMonteCarloModel.cs — MC batch runner + IMonteCarloModel
│   └── ScenarioVariation.cs    — parameter variation ranges
├── Analysis/
│   ├── ScenarioClusterAnalyzer.cs    — ML clustering of MC scenarios
│   ├── ProbabilityMap.cs             — per-cell exceedance probability
│   ├── TimeAnimator.cs               — time-animated probability maps
│   ├── ExposurePolygon.cs            — polygon result type
│   └── ExposurePolygonGenerator.cs   — marching squares contour extraction
└── Export/
    ├── GeoJsonExporter.cs      — GeoJSON (+ fire + contamination snapshots/extent/heatmap/centreline)
    ├── UnityBinaryExporter.cs  — compact binary for Unity3D (+ GFIR fire + GWCN contamination)
    └── CesiumExporter.cs       — CZML (+ fire + contamination time-dynamic animation)

Physics/Materials/Fire/
├── FuelModel.cs                — Rothermel fuel parameters (Anderson 13)
├── FuelLibrary.cs              — static registry of standard fuel models
└── Enums/
    └── FuelModelType.cs        — ShortGrass, Chaparral, …

Physics/Materials/Water/
├── AquaticContaminant.cs       — immutable contaminant descriptor (decay, Kd, toxicity)
├── ContaminantLibrary.cs       — static registry of 12 built-in contaminants
└── Enums/
    └── ContaminantType.cs      — Radioactive, Chemical, Biological, Thermal

Physics/Environmental/Fire/
└── RothermelModel.cs           — Rothermel (1972) equations (R, I_R, φw, φs)

Physics/Environmental/Water/
├── ManningEquation.cs          — open-channel velocity & discharge
├── LongitudinalDispersion.cs   — Fischer coefficient, shear velocity, decay/retardation
└── MixingZoneModel.cs          — tributary confluence mass-balance mixing
```

### Status

| Phase | Description | Tests | Status |
|-------|-------------|-------|--------|
| 0 | Grid foundation (`GeoGrid`, `GeoCell`, `GridSnapshot`, `TimeFrame`) | 25 | ✅ Done |
| 1 | Single-scenario physics (`PlumeSimulator`, `GaussianPuff`) | 14 | ✅ Done |
| 2 | Monte Carlo generation (`PlumeMonteCarloModel`, `ScenarioVariation`) | 10 | ✅ Done |
| 3 | ML clustering & probability maps (`ScenarioClusterAnalyzer`, `ProbabilityMap`, `TimeAnimator`) | 19 | ✅ Done |
| 4 | Fluent API (`RiskScenario` → `ScenarioResult`) | 14 | ✅ Done |
| 5 | Export (GeoJSON / Unity binary / Cesium CZML) | 24 | ✅ Done |
| 6 | GIS coordinates (WGS84, UTM, `GeoCoordinate`, `Projection`, `FromLatLon`) | 28 | ✅ Done |
| 7 | Radioactive fallout & nuclear materials (`Isotope`, `DecayChain`, `RadiationDose`, `.WithMaterial()`) | 58 | ✅ Done |
| 7.5 | Exposure polygons (`ExposurePolygonGenerator`, peak & time-integrated contours) | 17 | ✅ Done |
| **W1** | **Fire physics & fuel library** (`RothermelModel`, `FuelModel`, `FuelLibrary`) | **33** | ✅ Done |
| **W2** | **Terrain model & fuel map** (`TerrainGrid`, `FuelMap`) | **18** | ✅ Done |
| **W3** | **Cell-automaton fire spread** (`WildfireSimulator`, `SpreadSnapshot`) | **9** | ✅ Done |
| **W4** | **Fluent API & Monte Carlo** (`WildfireScenarioBuilder`, `AnalyzeWith`, clustering) | **13** | ✅ Done |
| **W5** | **Fire export** (GeoJSON fire features/perimeters/heatmap, CZML animation, Unity binary) | **14** | ✅ Done |
| **C1** | **Water physics & contaminant library** (`ManningEquation`, `LongitudinalDispersion`, `MixingZoneModel`, `AquaticContaminant`) | **31** | ✅ Done |
| **C2** | **River network & channel geometry** (`RiverNetwork`, `ChannelMap`) | **12** | ✅ Done |
| **C3** | **Advection-diffusion simulator** (`WaterContaminationSimulator`, `SpreadSnapshot` extensions) | **13** | ✅ Done |
| **C4** | **Fluent API & Monte Carlo** (`WaterContaminationScenarioBuilder`, `AnalyzeWith`, clustering) | **11** | ✅ Done |
| **C5** | **Contamination export** (GeoJSON extent/heatmap/centreline, CZML animation, Unity binary) | **7** | ✅ Done |

**Total: 370 GIS + nuclear + wildfire + water contamination tests**

See [ROADMAP.md](ROADMAP.md) for the full design and phase details.
