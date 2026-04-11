## 🪐 Exoplanet Engine

An analysis engine for exoplanet transit detection and ML-assisted classification. Combines transit physics, time-series analysis, feature extraction, and sequence model classification into a single pipeline — from raw light curve to transit prediction.

**Namespace:** `CSharpNumerics.Engines.Exoplanet`

---

### Pipeline Overview

```
LightCurve  →  Preprocess  →  BLS Period Search  →  Transit Fit  →  Validate  →  Feature Extract  →  ML Classify
```

The engine builds on existing library capabilities:

| Step | Uses |
|------|------|
| Preprocessing | `TimeSeriesDetrending`, `LightCurveSanitizer` |
| Period search | `BoxFittingLeastSquares`, `LombScarglePeriodogram` |
| Transit model | `Physics.Astro.TransitModel`, `LimbDarkening`, `TransitGeometry` |
| Fitting | `NonlinearLeastSquaresFitter`, `GoodnessOfFit` |
| Feature extraction | `TransitFeatureExtractor`, `WindowedFeatureExtractor` |
| ML classification | `CNN1DClassifier`, `LSTMClassifier`, `BiLSTMClassifier` via `SupervisedExperiment` |
| Statistics | `SigmaClipping`, `SlidingWindowStatistics`, `FalseAlarmProbability` |
| Engine infra | `ISimulationEngine`, `EventBus` from `Engines.Common` |

---

### Data Model

Light curves use a format compatible with Kepler/TESS FITS columns:

```csharp
using CSharpNumerics.Engines.Exoplanet.Data;

// Create from arrays
var lc = new LightCurve(time, flux, fluxError, qualityFlags,
    new LightCurveMetadata
    {
        TargetId = "TIC-12345678",
        Mission = "TESS",
        TimeOffset = 2457000.0,
        Cadence = CadenceType.Short
    });

// Or from a TimeSeries
var lc2 = LightCurve.FromTimeSeries(ts, "TIME", "PDCSAP_FLUX", "PDCSAP_FLUX_ERR");

// Sanitize
var clean = LightCurveSanitizer.RemoveBadQuality(lc, qualityMask: 0);
clean = LightCurveSanitizer.RemoveOutliers(clean, sigmaThreshold: 5.0);
clean = LightCurveSanitizer.NormalizeFlux(clean);
```

### Classical Detection

Detect transit candidates without ML — uses BLS period search, transit model fitting, and validation:

```csharp
using CSharpNumerics.Engines.Exoplanet.Pipeline;

var config = new TransitDetectionConfig
{
    MinPeriodDays = 0.5,
    MaxPeriodDays = 100.0,
    SnrThreshold = 7.0,
    MinTransitDepthPpm = 100,
    MaxPlanets = 3,
    DetrendingMethod = DetrendingMethod.MedianFilter,
    PeriodSearchMethod = PeriodSearchMethod.BLS
};

TransitCandidate[] candidates = TransitDetectionPipeline.Detect(lc, config);

foreach (var c in candidates)
{
    Console.WriteLine($"Period: {c.Parameters.Period:F4} d");
    Console.WriteLine($"Depth:  {c.Parameters.Depth * 1e6:F0} ppm");
    Console.WriteLine($"Score:  {c.Score:F3}");
}
```

### Feature Extraction

Extract transit-specific features for ML classification:

```csharp
using CSharpNumerics.Engines.Exoplanet.Features;

// Per-candidate features
TransitFeatureSet features = TransitFeatureExtractor.Extract(candidate, lc);
double blsSnr = features["SnrBls"];
double vShape = features["VShapeMetric"];

// Windowed training data (phase-folded + feature columns)
var (X, y) = WindowedFeatureExtractor.CreateTrainingData(
    curves, candidates, labels, windowSize: 50);
```

Features extracted: `Depth`, `Duration`, `Period`, `SnrBls`, `OddEvenRatio`, `VShapeMetric`, `SecondaryDepth`, `ScatterInTransit`, `ScatterOutTransit`, `LimbDarkeningU1`, `LimbDarkeningU2`, `IngressEgressRatio`.

### ML Training & Inference

Train a sequence classifier (CNN1D / LSTM / BiLSTM) on labeled transit data:

```csharp
using CSharpNumerics.Engines.Exoplanet.Pipeline;
using CSharpNumerics.ML.Experiment;

var trainerConfig = new TrainerConfig
{
    WindowSize = 50,
    CandidateModels = new[] { TransitModelType.CNN1D, TransitModelType.LSTM },
    Epochs = 150,
    LearningRate = 0.02,
    Filters = 8,
    KernelSize = 5,
    HiddenUnits = 8,
    Cv = CrossValidatorConfig.KFold(folds: 3)
};

// Train with grid search + cross-validation
TrainedTransitModel model = TransitClassifierTrainer.Train(
    curves, labels, trainerConfig, detectionConfig);

Console.WriteLine($"Best model: {model.ModelName}");
Console.WriteLine($"Accuracy:   {model.Metrics.Accuracy:P1}");
Console.WriteLine($"F1 Score:   {model.Metrics.F1Score:F3}");

// Serialize for deployment
byte[] bytes = ModelSerializer.Serialize(model);
ModelSerializer.SaveToFile(model, "transit_model.bin");

// Load and predict
var restored = ModelSerializer.LoadFromFile("transit_model.bin", model);
TransitPrediction[] predictions = restored.Predict(newLightCurve);
```

### Engine Integration

`ExoplanetEngine` implements `ISimulationEngine` for batch processing with event-driven notifications:

```csharp
using CSharpNumerics.Engines.Exoplanet;

var engineConfig = new ExoplanetEngineConfig
{
    Detection = new TransitDetectionConfig
    {
        MinPeriodDays = 0.5,
        MaxPeriodDays = 50.0,
        SnrThreshold = 7.0
    }
};

var engine = new ExoplanetEngine(engineConfig);

// Subscribe to detection events
engine.Bus.Subscribe<TransitDetectedEvent>(evt =>
{
    Console.WriteLine($"[t={evt.Timestamp:F1}s] Transit found: " +
        $"P={evt.Candidate.Parameters.Period:F3} d, " +
        $"depth={evt.Candidate.Parameters.Depth * 1e6:F0} ppm");
});

// Classical detection (no ML)
engine.Init();

// Or with a pre-trained model for ML-assisted classification:
// engine.Init(trainedModel);

// Process light curves
engine.Enqueue(lightCurve1);
engine.Enqueue(lightCurve2);
engine.Step(dt: 1.0);

Console.WriteLine($"Processed: {engine.ProcessedCount}");
Console.WriteLine($"Detections: {engine.Detections.Count}");
```

---

### Architecture

```
Engines/Exoplanet/
├── Data/
│   ├── LightCurve.cs              # Time-series photometry container
│   ├── LightCurveMetadata.cs      # Mission, cadence, target ID
│   ├── LightCurveSanitizer.cs     # Quality filtering, outlier removal
│   ├── TransitParameters.cs       # Period, epoch, depth, duration, etc.
│   ├── TransitCandidate.cs        # Candidate with parameters + score
│   ├── TransitFeatureSet.cs       # Named feature dictionary
│   └── StellarProperties.cs       # Teff, radius, mass, logg
├── Enums/
│   ├── CadenceType.cs             # Short / Long / Fast
│   ├── DetrendingMethod.cs        # MedianFilter / Polynomial / SavitzkyGolay
│   ├── PeriodSearchMethod.cs      # BLS / LombScargle / Both
│   ├── SpectralType.cs            # O B A F G K M
│   └── TransitDisposition.cs      # Unknown / Candidate / Confirmed / FalsePositive
├── Features/
│   ├── TransitFeatureExtractor.cs # 12 transit-specific features
│   └── WindowedFeatureExtractor.cs# Phase-folded windows + feature columns
├── Pipeline/
│   ├── TransitDetectionConfig.cs  # Detection thresholds and methods
│   ├── TransitDetectionPipeline.cs# Classical detection orchestrator
│   ├── LightCurvePreprocessor.cs  # Detrend → outlier removal → normalize
│   ├── PeriodSearcher.cs          # BLS / Lomb-Scargle wrapper
│   ├── TransitFitter.cs           # Non-linear transit model fitting
│   ├── TransitValidator.cs        # SNR, odd/even, secondary eclipse checks
│   ├── TrainerConfig.cs           # ML hyperparameters and model selection
│   ├── TransitClassifierTrainer.cs# Grid-search training façade
│   ├── TrainedTransitModel.cs     # Trained model + metrics wrapper
│   ├── TransitPrediction.cs       # Prediction result (candidate + probability)
│   ├── TransitInferencePipeline.cs# Stateless thread-safe inference
│   └── ModelSerializer.cs         # Binary serialization for deployment
├── ExoplanetEngine.cs             # ISimulationEngine implementation
├── ExoplanetEngineConfig.cs       # Engine configuration
└── TransitDetectedEvent.cs        # EventBus event type
```

### Cross-Section Dependencies

```
ExoplanetEngine
├── Engines.Common          # ISimulationEngine, EventBus
├── Physics.Astro           # TransitModel, LimbDarkening, TransitGeometry, KeplerOrbit
├── Statistics
│   ├── TimeSeriesAnalysis  # BLS, Lomb-Scargle, PhaseFolding, Detrending
│   ├── Fitting             # NonlinearLeastSquaresFitter, GoodnessOfFit
│   └── Robust              # SigmaClipping
├── ML
│   ├── Sequence            # CNN1D, LSTM, BiLSTM classifiers
│   └── Experiment          # SupervisedExperiment, PipelineGrid
└── Numerics.Objects        # Matrix, VectorN
```
