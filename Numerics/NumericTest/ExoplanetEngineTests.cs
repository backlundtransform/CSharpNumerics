using System;
using System.Collections.Generic;
using CSharpNumerics.Engines.Common;
using CSharpNumerics.Engines.Exoplanet;
using CSharpNumerics.Engines.Exoplanet.Data;
using CSharpNumerics.Engines.Exoplanet.Enums;
using CSharpNumerics.Engines.Exoplanet.Features;
using CSharpNumerics.Engines.Exoplanet.Pipeline;
using CSharpNumerics.ML.Experiment;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace NumericTest;

[TestClass]
public class ExoplanetEngineTests
{
    #region Helpers

    private static LightCurve CreateSyntheticTransitCurve(
        double period, double epoch, double depth, double durationDays,
        int numPoints = 2000, double cadenceDays = 0.02, double noiseLevel = 0.0002,
        int seed = 42)
    {
        var rng = new System.Random(seed);
        var time = new double[numPoints];
        var flux = new double[numPoints];
        var err = new double[numPoints];
        var qual = new int[numPoints];

        double halfDur = durationDays / 2.0;

        for (int i = 0; i < numPoints; i++)
        {
            time[i] = i * cadenceDays;
            flux[i] = 1.0;

            double phase = ((time[i] - epoch) % period + period) % period;
            if (phase > period / 2.0) phase -= period;
            if (Math.Abs(phase) <= halfDur)
                flux[i] = 1.0 - depth;

            flux[i] += (rng.NextDouble() * 2.0 - 1.0) * noiseLevel;
            err[i] = noiseLevel;
        }

        return new LightCurve(time, flux, err, qual, new LightCurveMetadata());
    }

    private static LightCurve CreateFlatCurve(
        int numPoints = 2000, double cadenceDays = 0.02, double noiseLevel = 0.0002,
        int seed = 99)
    {
        var rng = new System.Random(seed);
        var time = new double[numPoints];
        var flux = new double[numPoints];
        var err = new double[numPoints];
        var qual = new int[numPoints];

        for (int i = 0; i < numPoints; i++)
        {
            time[i] = i * cadenceDays;
            flux[i] = 1.0 + (rng.NextDouble() * 2.0 - 1.0) * noiseLevel;
            err[i] = noiseLevel;
        }

        return new LightCurve(time, flux, err, qual, new LightCurveMetadata());
    }

    #endregion

    // ── Engine lifecycle ──────────────────────────────────────────────

    [TestMethod]
    public void ExoplanetEngine_InitSetsInitialized()
    {
        var engine = new ExoplanetEngine(new ExoplanetEngineConfig());
        Assert.IsFalse(engine.IsInitialized);
        engine.Init();
        Assert.IsTrue(engine.IsInitialized);
    }

    [TestMethod]
    public void ExoplanetEngine_ResetClearsState()
    {
        var engine = new ExoplanetEngine(new ExoplanetEngineConfig());
        engine.Init();
        engine.Enqueue(CreateFlatCurve());
        engine.Step(1.0);

        engine.Reset();

        Assert.IsFalse(engine.IsInitialized);
        Assert.AreEqual(0.0, engine.Time);
        Assert.AreEqual(0, engine.ProcessedCount);
        Assert.AreEqual(0, engine.Detections.Count);
    }

    [TestMethod]
    [ExpectedException(typeof(InvalidOperationException))]
    public void ExoplanetEngine_StepBeforeInit_Throws()
    {
        var engine = new ExoplanetEngine(new ExoplanetEngineConfig());
        engine.Step(1.0);
    }

    [TestMethod]
    public void ExoplanetEngine_TimeAdvances()
    {
        var engine = new ExoplanetEngine(new ExoplanetEngineConfig());
        engine.Init();
        engine.Step(1.5);
        engine.Step(2.5);
        Assert.AreEqual(4.0, engine.Time, 1e-10);
    }

    // ── Classical detection path ──────────────────────────────────────

    [TestMethod]
    public void ExoplanetEngine_ClassicalDetection_FindsTransit()
    {
        var config = new ExoplanetEngineConfig
        {
            Detection = new TransitDetectionConfig
            {
                MinPeriodDays = 1.0,
                MaxPeriodDays = 10.0,
                SnrThreshold = 3.0,
                MinTransitDepthPpm = 50,
                MaxPlanets = 1
            }
        };

        var engine = new ExoplanetEngine(config);
        engine.Init();

        // Inject a transit with period=3 days, depth=0.01
        var lc = CreateSyntheticTransitCurve(
            period: 3.0, epoch: 0.5, depth: 0.01, durationDays: 0.15,
            numPoints: 3000, cadenceDays: 0.02, noiseLevel: 0.0002);

        engine.Enqueue(lc);
        engine.Step(1.0);

        Assert.AreEqual(1, engine.ProcessedCount);
        Assert.IsTrue(engine.Detections.Count >= 1, "Expected at least 1 transit detection");

        // Verify detected period is close to injected period
        var detected = engine.Detections[0];
        Assert.AreEqual(3.0, detected.Parameters.Period, 0.15,
            $"Detected period {detected.Parameters.Period:F4} should be ~3.0 days");
    }

    [TestMethod]
    public void ExoplanetEngine_FlatCurve_NoDetections()
    {
        var config = new ExoplanetEngineConfig
        {
            Detection = new TransitDetectionConfig
            {
                MinPeriodDays = 1.0,
                MaxPeriodDays = 10.0,
                SnrThreshold = 7.0,
                MinTransitDepthPpm = 100,
                MaxPlanets = 1
            }
        };

        var engine = new ExoplanetEngine(config);
        engine.Init();

        engine.Enqueue(CreateFlatCurve());
        engine.Step(1.0);

        Assert.AreEqual(0, engine.Detections.Count, "Flat curve should produce no detections");
    }

    // ── EventBus integration ──────────────────────────────────────────

    [TestMethod]
    public void ExoplanetEngine_PublishesTransitDetectedEvents()
    {
        var config = new ExoplanetEngineConfig
        {
            Detection = new TransitDetectionConfig
            {
                MinPeriodDays = 1.0,
                MaxPeriodDays = 10.0,
                SnrThreshold = 3.0,
                MinTransitDepthPpm = 50,
                MaxPlanets = 1
            }
        };

        var engine = new ExoplanetEngine(config);
        var events = new List<TransitDetectedEvent>();
        engine.Bus.Subscribe<TransitDetectedEvent>(e => events.Add(e));
        engine.Init();

        var lc = CreateSyntheticTransitCurve(
            period: 3.0, epoch: 0.5, depth: 0.01, durationDays: 0.15,
            numPoints: 3000, cadenceDays: 0.02, noiseLevel: 0.0002);

        engine.Enqueue(lc);
        engine.Step(5.0);

        Assert.IsTrue(events.Count >= 1, "Expected at least 1 TransitDetectedEvent");
        Assert.AreEqual(5.0, events[0].Timestamp, 1e-10);
        Assert.IsNotNull(events[0].Candidate);
    }

    // ── ML-assisted detection path ────────────────────────────────────

    [TestMethod]
    public void ExoplanetEngine_WithTrainedModel_RunsMLPath()
    {
        var detectionConfig = new TransitDetectionConfig
        {
            MinPeriodDays = 1.0,
            MaxPeriodDays = 10.0,
            SnrThreshold = 3.0,
            MinTransitDepthPpm = 50,
            MaxPlanets = 1
        };

        var trainerConfig = new TrainerConfig
        {
            WindowSize = 20,
            CandidateModels = new[] { TransitModelType.CNN1D },
            Epochs = 30,
            LearningRate = 0.02,
            Filters = 4,
            KernelSize = 3,
            HiddenUnits = 4,
            BatchSize = 4,
            Cv = CrossValidatorConfig.KFold(folds: 2)
        };

        // Create training data: 4 transit curves, 4 flat curves
        var curves = new LightCurve[8];
        var labels = new TransitDisposition[8];
        for (int i = 0; i < 4; i++)
        {
            curves[i] = CreateSyntheticTransitCurve(
                period: 3.0, epoch: 0.5, depth: 0.01, durationDays: 0.15,
                numPoints: 2000, cadenceDays: 0.02, noiseLevel: 0.0002, seed: 100 + i);
            labels[i] = TransitDisposition.Confirmed;
        }
        for (int i = 4; i < 8; i++)
        {
            curves[i] = CreateFlatCurve(numPoints: 2000, seed: 200 + i);
            labels[i] = TransitDisposition.FalsePositive;
        }

        // Train model
        var model = TransitClassifierTrainer.Train(curves, labels, trainerConfig, detectionConfig);
        Assert.IsNotNull(model);

        // Init engine with trained model
        var engineConfig = new ExoplanetEngineConfig
        {
            Detection = detectionConfig,
            Training = trainerConfig
        };

        var engine = new ExoplanetEngine(engineConfig);
        engine.Init(model);
        Assert.IsTrue(engine.IsInitialized);

        // Process a new transit curve
        var testCurve = CreateSyntheticTransitCurve(
            period: 3.0, epoch: 0.5, depth: 0.01, durationDays: 0.15,
            numPoints: 2000, cadenceDays: 0.02, noiseLevel: 0.0002, seed: 999);

        engine.Enqueue(testCurve);
        engine.Step(1.0);

        Assert.AreEqual(1, engine.ProcessedCount);
        // ML path should produce predictions (may or may not detect depending on model quality)
    }

    // ── Multiple light curves ────────────────────────────────────────

    [TestMethod]
    public void ExoplanetEngine_ProcessesMultipleCurvesInOneStep()
    {
        var config = new ExoplanetEngineConfig
        {
            Detection = new TransitDetectionConfig
            {
                MinPeriodDays = 1.0,
                MaxPeriodDays = 10.0,
                SnrThreshold = 3.0,
                MinTransitDepthPpm = 50,
                MaxPlanets = 1
            }
        };

        var engine = new ExoplanetEngine(config);
        engine.Init();

        engine.Enqueue(CreateSyntheticTransitCurve(
            period: 3.0, epoch: 0.5, depth: 0.01, durationDays: 0.15, seed: 1));
        engine.Enqueue(CreateSyntheticTransitCurve(
            period: 5.0, epoch: 0.3, depth: 0.008, durationDays: 0.20, seed: 2));

        engine.Step(1.0);

        Assert.AreEqual(2, engine.ProcessedCount);
    }

    // ── Config defaults ──────────────────────────────────────────────

    [TestMethod]
    public void ExoplanetEngineConfig_HasDefaults()
    {
        var config = new ExoplanetEngineConfig();
        Assert.IsNotNull(config.Detection);
        Assert.IsNotNull(config.Training);
        Assert.IsNull(config.ModelPath);
    }

    [TestMethod]
    public void TransitDetectedEvent_StoresValues()
    {
        var candidate = new TransitCandidate(
            new TransitParameters { Period = 3.0, Depth = 0.01 },
            0.95, TransitDisposition.Candidate,
            CreateFlatCurve(numPoints: 50),
            new TransitFeatureSet());

        var evt = new TransitDetectedEvent(candidate, 42.0);
        Assert.AreSame(candidate, evt.Candidate);
        Assert.AreEqual(42.0, evt.Timestamp, 1e-10);
    }

    // ── End-to-end demo test ─────────────────────────────────────────

    [TestMethod]
    public void Demo_EndToEnd_TrainAndInfer()
    {
        // ── 1. Configure ──────────────────────────────────────────
        var detectionConfig = new TransitDetectionConfig
        {
            MinPeriodDays = 1.0,
            MaxPeriodDays = 10.0,
            SnrThreshold = 3.0,
            MinTransitDepthPpm = 50,
            MaxPlanets = 1
        };

        var trainerConfig = new TrainerConfig
        {
            WindowSize = 20,
            CandidateModels = new[] { TransitModelType.CNN1D },
            Epochs = 30,
            LearningRate = 0.02,
            Filters = 4,
            KernelSize = 3,
            HiddenUnits = 4,
            BatchSize = 4,
            Cv = CrossValidatorConfig.KFold(folds: 2)
        };

        // ── 2. Synthesize training data ───────────────────────────
        // Simulate TESS-like cadence: 2-min exposures (~0.00139 days) over 27 days
        // For test speed we use coarser cadence
        var trainCurves = new LightCurve[6];
        var trainLabels = new TransitDisposition[6];

        // 3 transit curves with different parameters
        trainCurves[0] = CreateSyntheticTransitCurve(3.0, 0.5, 0.01, 0.15, 2000, 0.02, 0.0002, 10);
        trainLabels[0] = TransitDisposition.Confirmed;
        trainCurves[1] = CreateSyntheticTransitCurve(3.0, 0.5, 0.01, 0.15, 2000, 0.02, 0.0002, 11);
        trainLabels[1] = TransitDisposition.Confirmed;
        trainCurves[2] = CreateSyntheticTransitCurve(3.0, 0.5, 0.01, 0.15, 2000, 0.02, 0.0002, 12);
        trainLabels[2] = TransitDisposition.Confirmed;

        // 3 flat (no transit) curves
        trainCurves[3] = CreateFlatCurve(2000, 0.02, 0.0002, 20);
        trainLabels[3] = TransitDisposition.FalsePositive;
        trainCurves[4] = CreateFlatCurve(2000, 0.02, 0.0002, 21);
        trainLabels[4] = TransitDisposition.FalsePositive;
        trainCurves[5] = CreateFlatCurve(2000, 0.02, 0.0002, 22);
        trainLabels[5] = TransitDisposition.FalsePositive;

        // ── 3. Train model ────────────────────────────────────────
        var model = TransitClassifierTrainer.Train(
            trainCurves, trainLabels, trainerConfig, detectionConfig);

        Assert.IsNotNull(model, "Training should produce a model");
        Assert.IsNotNull(model.ModelName, "Model name should be set");

        // ── 4. Serialize → Deserialize (roundtrip) ────────────────
        byte[] serialized = ModelSerializer.Serialize(model);
        Assert.IsTrue(serialized.Length > 0, "Serialized model should not be empty");

        var restored = ModelSerializer.Deserialize(serialized, model);
        Assert.IsNotNull(restored);

        // ── 5. Run inference via engine ────────────────────────────
        var engineConfig = new ExoplanetEngineConfig
        {
            Detection = detectionConfig,
            Training = trainerConfig
        };

        var engine = new ExoplanetEngine(engineConfig);
        var detectedEvents = new List<TransitDetectedEvent>();
        engine.Bus.Subscribe<TransitDetectedEvent>(e => detectedEvents.Add(e));
        engine.Init(restored);

        // Inference on a new transit curve
        var newTransitCurve = CreateSyntheticTransitCurve(
            3.0, 0.5, 0.01, 0.15, 2000, 0.02, 0.0002, 555);
        engine.Enqueue(newTransitCurve);
        engine.Step(1.0);

        // Inference on a flat curve
        var newFlatCurve = CreateFlatCurve(2000, 0.02, 0.0002, 666);
        engine.Enqueue(newFlatCurve);
        engine.Step(1.0);

        Assert.AreEqual(2, engine.ProcessedCount, "Both curves should be processed");

        // ── 6. Results format (what a web service would return) ───
        // Each detection has: candidate parameters + timestamp
        foreach (var evt in detectedEvents)
        {
            Assert.IsNotNull(evt.Candidate.Parameters);
            Assert.IsTrue(evt.Timestamp > 0);
        }
    }
}
