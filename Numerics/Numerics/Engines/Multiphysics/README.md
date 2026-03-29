## ⚡ Multiphysics Engine

A fluent simulation engine for four student-friendly PDE simulations on simple geometries, designed for visualization in Unity or web platforms. Built on existing finite difference infrastructure (`Grid2D`, `GridOperators`, `ITimeStepper`) with engineering materials, binary/JSON export, and optional ML/Monte Carlo integration.

**Namespace:** `CSharpNumerics.Engines.Multiphysics`

---

### Pipeline Overview

```
Material → Geometry → Boundary → Source → Solve → Result → Export / ML
```

The engine builds on existing library capabilities:

| Step | Uses |
|------|------|
| Materials | `EngineeringMaterial`, `EngineeringLibrary` (Physics) |
| FD Spatial | `Grid2D`, `GridOperators.Laplacian2D`, `GridOperators.SolvePoisson2D` (Numerics) |
| Time Integration | Forward Euler (internal), `ITimeStepper` compatible |
| Beam Theory | `SolidExtensions` — Euler-Bernoulli, Hooke's law (Physics) |
| Monte Carlo | `IMonteCarloModel`, `MonteCarloSimulator` (Statistics) |
| ML Clustering | `ClusteringExperiment`, `KMeans`, evaluators (ML) |
| ML Surrogate | `IModel`, `Ridge`, `MLPRegressor` (ML) |

### Simulation Types

| Type | PDE | Method | Dimension |
|------|-----|--------|-----------|
| `HeatPlate` | $\partial T/\partial t = \alpha \nabla^2 T + S$ | FD (Grid2D + Laplacian2D + forward Euler) | 2D |
| `PipeFlow` | Hagen-Poiseuille (cylindrical Laplacian) | FD (1D transient) | 1D |
| `ElectricField` | $\nabla^2 \varphi = -\rho/\varepsilon$ | Poisson solver (Gauss-Seidel) | 2D |
| `BeamStress` | $EI u'''' = q$ | Analytical (SolidExtensions) | 1D |

---

### Fluent API — SimulationType → SimulationBuilder → SimulationResult

All simulations use the same fluent entry point:

```csharp
using CSharpNumerics.Engines.Multiphysics;
using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Physics.Materials.Engineering;

var result = SimulationType.Create(MultiphysicsType.HeatPlate)
    .WithMaterial(EngineeringLibrary.Aluminum)
    .WithGeometry(width: 0.1, height: 0.1, nx: 20, ny: 20)
    .WithBoundary(top: 100.0, bottom: 0.0, left: 0.0, right: 0.0)
    .WithInitialCondition(0.0)
    .AddSource(10, 10, 1e6)
    .Solve(dt: 0.0001, steps: 500);

double[,] finalField = result.Field;           // final temperature field
List<double[,]> timeline = result.Timeline;    // per-step snapshots
double peak = result.MaxValue;                 // 100.0
```

### HeatPlate — 2D Transient Heat Equation

$$\frac{\partial T}{\partial t} = \alpha \nabla^2 T + \frac{S}{\rho c_p}$$

Thermal diffusivity $\alpha = k / (\rho c_p)$ is automatically derived from the `EngineeringMaterial`.

```csharp
var result = SimulationType.Create(MultiphysicsType.HeatPlate)
    .WithMaterial(EngineeringLibrary.Copper)       // α = k/(ρ·cp) ≈ 1.17e-4
    .WithGeometry(width: 0.1, height: 0.1, nx: 50, ny: 50)
    .WithBoundary(top: 200.0, bottom: 0.0, left: 0.0, right: 0.0)
    .WithInitialCondition((x, y) => 20.0)          // function-based IC
    .AddSource(25, 25, 5e5)                         // point heat source
    .Solve(dt: 0.00005, steps: 1000);

// Access timeline for animation
for (int step = 0; step < result.Timeline.Count; step++)
    double[,] frame = result.Timeline[step];
```

### PipeFlow — 1D Hagen-Poiseuille

Transient cylindrical channel flow converging to the parabolic velocity profile:

$$v(r) = \frac{R^2 - r^2}{4\mu} \left(-\frac{dP}{dx}\right)$$

```csharp
var result = SimulationType.Create(MultiphysicsType.PipeFlow)
    .WithMaterial(EngineeringLibrary.Water)        // uses KinematicViscosity
    .WithGeometry(length: 1.0, radius: 0.005, nodes: 21)
    .WithBoundary(pressureGradient: -100.0)
    .Solve(dt: 0.02, steps: 3000);

double[] velocity = result.Values;     // radial velocity profile
double[] positions = result.Positions; // radial positions [0..R]
double vCenter = result.Values[0];     // peak centreline velocity
```

### ElectricField — 2D Poisson/Laplace

Solves $\nabla^2 \varphi = -\rho/\varepsilon$ with Dirichlet BCs, then computes $\mathbf{E} = -\nabla \varphi$:

```csharp
var result = SimulationType.Create(MultiphysicsType.ElectricField)
    .WithMaterial(EngineeringLibrary.Air)          // uses ElectricPermittivity
    .WithGeometry(width: 1.0, height: 1.0, nx: 40, ny: 40)
    .WithBoundary(top: 0.0, bottom: 0.0, left: 0.0, right: 100.0)
    .AddSource(20, 20, 1e-6)                        // point charge
    .Solve(maxIterations: 20000, tolerance: 1e-8);

double[,] potential = result.Field;   // φ(x, y)
double[,] Ex = result.Ex;            // -∂φ/∂x
double[,] Ey = result.Ey;            // -∂φ/∂y
```

### BeamStress — 1D Euler-Bernoulli

Analytical beam solutions with multiple support types and load combinations:

```csharp
var result = SimulationType.Create(MultiphysicsType.BeamStress)
    .WithMaterial(EngineeringLibrary.Steel)         // uses YoungsModulus
    .WithGeometry(length: 2.0, nodes: 201)
    .WithCrossSection(width: 0.05, height: 0.1)    // rectangular
    // .WithCrossSection(radius: 0.02)              // circular
    // .WithSecondMoment(1.5e-6)                    // custom I
    .WithBoundary(BeamSupport.Cantilever)           // or SimplySupported, FixedFixed
    .AddSource(position: 2.0, value: 1000)          // point load at tip
    .WithSource(500)                                 // uniform distributed load
    .Solve();

double[] deflection = result.Values;        // δ(x)
double[] moment = result.BendingMoment;     // M(x)
double[] shear = result.ShearForce;         // V(x)
double[] stress = result.Stress;            // σ(x) = M·y_max/I
double[] x = result.Positions;             // node positions
```

#### Analytical Reference

| Support | Point Load P at tip | Uniform Load q |
|---------|--------------------|--------------  |
| Cantilever | $\delta_{max} = \frac{PL^3}{3EI}$ | $\delta_{max} = \frac{qL^4}{8EI}$ |
| Simply Supported | — | $\delta_{mid} = \frac{5qL^4}{384EI}$ |
| Fixed-Fixed | — | $\delta_{mid} = \frac{qL^4}{384EI}$ |

---

### Snapshots & Timeline

Time-stamped field data for animation and export:

```csharp
using CSharpNumerics.Engines.Multiphysics.Snapshots;

// Build timeline from simulation result
var timeline = SimulationTimeline.FromResult(result, dt: 0.0001, dx: 0.005, dy: 0.005);

int frames = timeline.Count;                     // initial + N steps
double t0 = timeline.StartTime;
double tEnd = timeline.EndTime;
FieldSnapshot first = timeline[0];               // access by index

// Interpolate between frames (linear blend)
double[,] interp = timeline.InterpolateAt(0.00015);

// Individual snapshots
FieldSnapshot snap = timeline[5];
double val = snap[10, 15];                       // indexing
double max = snap.Max();
double min = snap.Min();

// Beam snapshot
var beamSnap = BeamSnapshot.FromResult(result, BeamSupport.Cantilever, length: 2.0);
double maxDef = beamSnap.MaxDeflection;
double maxStress = beamSnap.MaxStress;
```

---

### Binary Export (MPHY Format)

Unity-compatible binary format with magic bytes, header, and float layers:

```csharp
using CSharpNumerics.Engines.Multiphysics.Export;

// Save 2D timeline
MultiphysicsBinaryExporter.Save(timeline, "heat_sim.mphy");

// Save beam snapshot
MultiphysicsBinaryExporter.Save(beamSnap, "beam.mphy");

// Read back
var header = MultiphysicsBinaryExporter.ReadHeader("heat_sim.mphy");
// header.Type, header.Nx, header.Ny, header.TimeStepCount, header.LayerCount

var data = MultiphysicsBinaryExporter.Read("heat_sim.mphy");
float[] frame0 = data.Layers[0];    // first time step, flattened
```

### JSON Export

Zero-dependency JSON for web visualization:

```csharp
// Direct result export
string json = MultiphysicsJsonExporter.ToJson(result);
MultiphysicsJsonExporter.Save(result, "simulation.json");

// Timeline export (all frames)
string timelineJson = MultiphysicsJsonExporter.ToJson(timeline);

// Single snapshot
string snapJson = MultiphysicsJsonExporter.ToJson(snapshot);

// Beam snapshot with metadata
string beamJson = MultiphysicsJsonExporter.ToJson(beamSnap);
```

---

### Monte Carlo Integration

Stochastic parameter variation for uncertainty quantification:

```csharp
using CSharpNumerics.Engines.Multiphysics.MonteCarlo;

// Define parameter variation ranges
var variation = new ParameterVariation()
    .ThermalConductivity(40, 60)        // k ∈ [40, 60] W/(m·K)
    .BoundaryTemperature(90, 110)       // T_top ∈ [90, 110] °C
    .SourceIntensity(8e5, 1.2e6);       // S ∈ [8e5, 1.2e6] W/m³

// Build MC model from a baseline simulation
var baseline = SimulationType.Create(MultiphysicsType.HeatPlate)
    .WithMaterial(EngineeringLibrary.Steel)
    .WithGeometry(width: 0.1, height: 0.1, nx: 10, ny: 10)
    .WithBoundary(top: 100.0, bottom: 0, left: 0, right: 0)
    .WithInitialCondition(0.0);

var mc = new MultiphysicsMonteCarloModel(baseline, variation);

// Full batch: scenario matrix + per-iteration results
MultiphysicsScenarioResult mcResult = mc.RunBatch(iterations: 100, seed: 42);

Matrix scenarioMatrix = mcResult.ScenarioMatrix;  // 100 × (Nx·Ny)
int features = mcResult.FeatureCount;              // Nx·Ny

// Per-cell distribution across scenarios
double[] dist = mcResult.GetFeatureDistribution(featureIndex: 50);

// Percentile maps
double[] p95 = mcResult.ComputePercentile(95);
double[] p50 = mcResult.ComputePercentile(50);
```

**Compatible with `MonteCarloSimulator`** for scalar summary:

```csharp
using CSharpNumerics.Statistics.MonteCarlo;

var sim = new MonteCarloSimulator(seed: 42);
MonteCarloResult stats = sim.Run(mc, iterations: 1000);
double meanPeak = stats.Mean;
double p95scalar = stats.Percentile(95);
```

### Surrogate Model

Train a regression model on MC data for fast prediction without re-running the solver:

```csharp
var surrogate = new SurrogateTrainer(mcResult, new Ridge());
surrogate.Train(r => r.MaxValue);       // target = peak temperature

// Fast prediction without running the solver
VectorN predictions = surrogate.Predict(mcResult.ScenarioMatrix);
double single = surrogate.PredictOne(mcResult.GetScenarioVector(0));
```

### Cluster Analysis

Identify representative scenario groups from Monte Carlo ensembles:

```csharp
using CSharpNumerics.ML.Clustering.Algorithms;
using CSharpNumerics.ML.Clustering.Evaluators;

var analysis = MultiphysicsClusterAnalyzer
    .For(mcResult)
    .WithAlgorithm(new KMeans { K = 3 })
    .TryClusterCounts(2, 5)
    .WithEvaluator(new SilhouetteEvaluator())
    .Run();

int bestK = analysis.BestClusterCount;
int dominant = analysis.DominantCluster;             // most frequent group
int[] iterations = analysis.GetClusterIterations(dominant);
double[] meanOutput = analysis.GetClusterMeanOutput(dominant);
```

---

### Architecture

```
Engines/Multiphysics/
├── SimulationType.cs             ← static entry point
├── SimulationBuilder.cs          ← fluent builder chain
├── SimulationResult.cs           ← unified result type
├── Enums/
│   ├── MultiphysicsType.cs       ← HeatPlate | PipeFlow | ElectricField | BeamStress
│   └── BeamSupport.cs            ← Cantilever | SimplySupported | FixedFixed
├── Solvers/
│   ├── IMultiphysicsSolver.cs    ← internal solver interface
│   ├── HeatPlateSolver.cs        ← 2D FD heat equation
│   ├── PipeFlowSolver.cs         ← 1D cylindrical Laplacian
│   ├── ElectricFieldSolver.cs    ← 2D Poisson + gradient
│   └── BeamStressSolver.cs       ← analytical Euler-Bernoulli
├── Snapshots/
│   ├── FieldSnapshot.cs          ← time-stamped 2D scalar field
│   ├── BeamSnapshot.cs           ← immutable beam curves
│   └── SimulationTimeline.cs     ← ordered snapshot collection
├── Export/
│   ├── MultiphysicsBinaryExporter.cs  ← MPHY binary format
│   └── MultiphysicsJsonExporter.cs    ← zero-dependency JSON
└── MonteCarlo/
    ├── ParameterVariation.cs          ← stochastic parameter ranges
    ├── MultiphysicsMonteCarloModel.cs ← IMonteCarloModel, RunBatch
    ├── MultiphysicsScenarioResult.cs  ← scenario matrix + results
    ├── SurrogateTrainer.cs            ← regression surrogate
    └── MultiphysicsClusterAnalyzer.cs ← clustering wrapper
```
