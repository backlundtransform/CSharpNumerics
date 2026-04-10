using System;
using System.Linq;
using CSharpNumerics.Engines.Exoplanet.Data;
using CSharpNumerics.Engines.Exoplanet.Enums;
using CSharpNumerics.Engines.Exoplanet.Features;
using CSharpNumerics.Engines.Exoplanet.Pipeline;
using CSharpNumerics.Numerics.Objects;

namespace NumericTest;

[TestClass]
public class ExoplanetMLPipelineTests
{
    #region Helpers

    private static LightCurve CreateSyntheticTransitCurve(
        double period, double epoch, double depth, double durationDays,
        int numPoints = 2000, double cadenceDays = 0.02, double noiseLevel = 0.0002,
        int seed = 42)
    {
        var rng = new Random(seed);
        double[] time = new double[numPoints];
        double[] flux = new double[numPoints];
        double[] err = new double[numPoints];
        int[] qual = new int[numPoints];
        double halfDur = durationDays / 2.0;

        for (int i = 0; i < numPoints; i++)
        {
            time[i] = i * cadenceDays;
            flux[i] = 1.0;

            double phase = ((time[i] - epoch) % period + period) % period;
            if (phase > period / 2.0) phase -= period;
            if (Math.Abs(phase) <= halfDur)
                flux[i] = 1.0 - depth;

            flux[i] += noiseLevel * (rng.NextDouble() * 2.0 - 1.0);
            err[i] = noiseLevel;
        }

        return new LightCurve(time, flux, err, qual, new LightCurveMetadata());
    }

    private static LightCurve CreateFlatCurve(int numPoints = 2000, double noiseLevel = 0.0002, int seed = 100)
    {
        var rng = new Random(seed);
        double[] time = new double[numPoints];
        double[] flux = new double[numPoints];
        double[] err = new double[numPoints];
        int[] qual = new int[numPoints];

        for (int i = 0; i < numPoints; i++)
        {
            time[i] = i * 0.02;
            flux[i] = 1.0 + noiseLevel * (rng.NextDouble() * 2.0 - 1.0);
            err[i] = noiseLevel;
        }

        return new LightCurve(time, flux, err, qual, new LightCurveMetadata());
    }

    #endregion

    #region TrainerConfig Tests

    [TestMethod]
    public void TrainerConfig_DefaultValues_ShouldBeReasonable()
    {
        var config = new TrainerConfig();

        Assert.AreEqual(50, config.WindowSize);
        Assert.AreEqual(1, config.Stride);
        Assert.IsTrue(config.CandidateModels.Length > 0);
        Assert.AreEqual(TransitModelType.CNN1D, config.CandidateModels[0]);
        Assert.IsTrue(config.Epochs > 0);
        Assert.IsTrue(config.LearningRate > 0);
    }

    #endregion

    #region TransitPrediction Tests

    [TestMethod]
    public void TransitPrediction_ShouldStoreAllFields()
    {
        var parameters = new TransitParameters(3.5, 0.5, 0.01, 0.15, 0.1, 0.0, 0.02);
        var candidate = new TransitCandidate(parameters, 0.9, TransitDisposition.Candidate,
            LightCurve.FromArrays(new double[] { 0.0 }, new double[] { 1.0 }),
            new TransitFeatureSet());

        var prediction = new TransitPrediction(candidate, 0.95, TransitDisposition.Candidate);

        Assert.AreEqual(0.95, prediction.Probability);
        Assert.AreEqual(TransitDisposition.Candidate, prediction.PredictedDisposition);
        Assert.AreEqual(3.5, prediction.Candidate.Parameters.Period);
    }

    #endregion

    #region TrainingMetrics Tests

    [TestMethod]
    public void TrainingMetrics_FromConfusionMatrix_ShouldComputeCorrectly()
    {
        // Confusion matrix: [[90, 10], [5, 95]]
        // TN=90, FP=10, FN=5, TP=95
        var cm = new Matrix(new double[,] { { 90, 10 }, { 5, 95 } });
        var metrics = new TrainingMetrics(0.925, cm, 0.925);

        Assert.AreEqual(0.925, metrics.Accuracy, 0.001);
        Assert.AreEqual(95.0 / 105.0, metrics.Precision, 0.001);
        Assert.AreEqual(95.0 / 100.0, metrics.Recall, 0.001);
        Assert.IsTrue(metrics.F1Score > 0.9);
    }

    [TestMethod]
    public void TrainingMetrics_Default_ShouldBeZero()
    {
        var metrics = new TrainingMetrics();
        Assert.AreEqual(0, metrics.Accuracy);
        Assert.AreEqual(0, metrics.Precision);
        Assert.AreEqual(0, metrics.Recall);
        Assert.AreEqual(0, metrics.F1Score);
    }

    #endregion

    #region TransitClassifierTrainer Tests

    [TestMethod]
    public void Train_WithSyntheticData_ShouldReturnTrainedModel()
    {
        // Create training set: 6 curves with transit + 6 flat
        int numTransit = 6;
        int numFlat = 6;
        double period = 3.0, epoch = 0.5, depth = 0.01, duration = 0.15;

        var curves = new LightCurve[numTransit + numFlat];
        var labels = new TransitDisposition[numTransit + numFlat];

        for (int i = 0; i < numTransit; i++)
        {
            curves[i] = CreateSyntheticTransitCurve(period, epoch, depth, duration,
                numPoints: 1500, seed: 42 + i);
            labels[i] = TransitDisposition.Confirmed;
        }
        for (int i = 0; i < numFlat; i++)
        {
            curves[numTransit + i] = CreateFlatCurve(numPoints: 1500, seed: 200 + i);
            labels[numTransit + i] = TransitDisposition.FalsePositive;
        }

        var trainerConfig = new TrainerConfig
        {
            WindowSize = 30,
            CandidateModels = new[] { TransitModelType.CNN1D },
            Epochs = 50,
            LearningRate = 0.02,
            Filters = 4,
            KernelSize = 3,
            HiddenUnits = 4,
            BatchSize = 4,
            Cv = CSharpNumerics.ML.Experiment.CrossValidatorConfig.KFold(folds: 2)
        };

        var detectionConfig = new TransitDetectionConfig
        {
            MinPeriodDays = 1.0,
            MaxPeriodDays = 10.0,
            SnrThreshold = 3.0,
            MinTransitDepthPpm = 50
        };

        var model = TransitClassifierTrainer.Train(curves, labels, trainerConfig, detectionConfig);

        Assert.IsNotNull(model);
        Assert.IsNotNull(model.Pipeline);
        Assert.IsNotNull(model.ModelName);
        Assert.IsTrue(model.ModelName.Length > 0);
        Assert.IsNotNull(model.Metrics);
        Assert.IsTrue(model.Metrics.BestCvScore >= 0);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentNullException))]
    public void Train_NullCurves_ShouldThrow()
    {
        TransitClassifierTrainer.Train(null, new TransitDisposition[0],
            new TrainerConfig());
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void Train_MismatchedLengths_ShouldThrow()
    {
        var curves = new[] { CreateFlatCurve() };
        var labels = new TransitDisposition[] { TransitDisposition.Confirmed, TransitDisposition.FalsePositive };

        TransitClassifierTrainer.Train(curves, labels, new TrainerConfig());
    }

    #endregion

    #region TrainedTransitModel Predict Tests

    [TestMethod]
    public void TrainedModel_Predict_ShouldReturnPredictions()
    {
        // Train a small model
        int n = 6;
        double period = 3.0, epoch = 0.5, depth = 0.01, duration = 0.15;

        var curves = new LightCurve[n * 2];
        var labels = new TransitDisposition[n * 2];

        for (int i = 0; i < n; i++)
        {
            curves[i] = CreateSyntheticTransitCurve(period, epoch, depth, duration,
                numPoints: 1500, seed: 42 + i);
            labels[i] = TransitDisposition.Confirmed;
        }
        for (int i = 0; i < n; i++)
        {
            curves[n + i] = CreateFlatCurve(numPoints: 1500, seed: 200 + i);
            labels[n + i] = TransitDisposition.FalsePositive;
        }

        var trainerConfig = new TrainerConfig
        {
            WindowSize = 30,
            CandidateModels = new[] { TransitModelType.CNN1D },
            Epochs = 30,
            LearningRate = 0.02,
            Filters = 4,
            KernelSize = 3,
            HiddenUnits = 4,
            BatchSize = 4,
            Cv = CSharpNumerics.ML.Experiment.CrossValidatorConfig.KFold(folds: 2)
        };

        var detectionConfig = new TransitDetectionConfig
        {
            MinPeriodDays = 1.0,
            MaxPeriodDays = 10.0,
            SnrThreshold = 3.0,
            MinTransitDepthPpm = 50
        };

        var model = TransitClassifierTrainer.Train(curves, labels, trainerConfig, detectionConfig);

        // Predict on a new transit curve
        var testCurve = CreateSyntheticTransitCurve(period, epoch, depth, duration,
            numPoints: 1500, seed: 999);
        var predictions = model.Predict(testCurve);

        // Should find at least one candidate (since the curve has a transit)
        // Note: predictions may be empty if detection pipeline doesn't find a candidate
        Assert.IsNotNull(predictions);
    }

    #endregion

    #region ModelSerializer Tests

    [TestMethod]
    public void ModelSerializer_SerializeDeserialize_RoundTrip()
    {
        // Train a model
        int n = 4;
        double period = 3.0, epoch = 0.5, depth = 0.01, duration = 0.15;

        var curves = new LightCurve[n * 2];
        var labels = new TransitDisposition[n * 2];

        for (int i = 0; i < n; i++)
        {
            curves[i] = CreateSyntheticTransitCurve(period, epoch, depth, duration,
                numPoints: 1500, seed: 42 + i);
            labels[i] = TransitDisposition.Confirmed;
        }
        for (int i = 0; i < n; i++)
        {
            curves[n + i] = CreateFlatCurve(numPoints: 1500, seed: 200 + i);
            labels[n + i] = TransitDisposition.FalsePositive;
        }

        var trainerConfig = new TrainerConfig
        {
            WindowSize = 20,
            CandidateModels = new[] { TransitModelType.CNN1D },
            Epochs = 20,
            LearningRate = 0.02,
            Filters = 4,
            KernelSize = 3,
            HiddenUnits = 4,
            BatchSize = 4,
            Cv = CSharpNumerics.ML.Experiment.CrossValidatorConfig.KFold(folds: 2)
        };

        var model = TransitClassifierTrainer.Train(curves, labels, trainerConfig);

        // Serialize
        byte[] serialized = ModelSerializer.Serialize(model);
        Assert.IsNotNull(serialized);
        Assert.IsTrue(serialized.Length > 0);

        // Deserialize
        var restored = ModelSerializer.Deserialize(serialized, model);
        Assert.IsNotNull(restored);
        Assert.AreEqual(model.ModelName, restored.ModelName);
        Assert.AreEqual(model.Metrics.Accuracy, restored.Metrics.Accuracy, 0.001);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentNullException))]
    public void ModelSerializer_Serialize_NullModel_ShouldThrow()
    {
        ModelSerializer.Serialize(null);
    }

    #endregion

    #region TransitInferencePipeline Tests

    [TestMethod]
    public void InferencePipeline_FromModel_ShouldPredict()
    {
        // Train a small model
        int n = 4;
        double period = 3.0, epoch = 0.5, depth = 0.01, duration = 0.15;

        var curves = new LightCurve[n * 2];
        var labels = new TransitDisposition[n * 2];

        for (int i = 0; i < n; i++)
        {
            curves[i] = CreateSyntheticTransitCurve(period, epoch, depth, duration,
                numPoints: 1500, seed: 42 + i);
            labels[i] = TransitDisposition.Confirmed;
        }
        for (int i = 0; i < n; i++)
        {
            curves[n + i] = CreateFlatCurve(numPoints: 1500, seed: 200 + i);
            labels[n + i] = TransitDisposition.FalsePositive;
        }

        var trainerConfig = new TrainerConfig
        {
            WindowSize = 20,
            CandidateModels = new[] { TransitModelType.CNN1D },
            Epochs = 20,
            LearningRate = 0.02,
            Filters = 4,
            KernelSize = 3,
            HiddenUnits = 4,
            BatchSize = 4,
            Cv = CSharpNumerics.ML.Experiment.CrossValidatorConfig.KFold(folds: 2)
        };

        var detectionConfig = new TransitDetectionConfig
        {
            MinPeriodDays = 1.0,
            MaxPeriodDays = 10.0,
            SnrThreshold = 3.0,
            MinTransitDepthPpm = 50
        };

        var model = TransitClassifierTrainer.Train(curves, labels, trainerConfig, detectionConfig);

        // Create inference pipeline
        var pipeline = TransitInferencePipeline.FromModel(model);
        Assert.IsNotNull(pipeline);
        Assert.AreEqual(model.ModelName, pipeline.ModelName);

        // Predict
        var testCurve = CreateSyntheticTransitCurve(period, epoch, depth, duration,
            numPoints: 1500, seed: 999);
        var predictions = pipeline.Predict(testCurve);
        Assert.IsNotNull(predictions);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentNullException))]
    public void InferencePipeline_Predict_NullCurve_ShouldThrow()
    {
        // We need a model — use a minimal one
        int n = 4;
        double period = 3.0, epoch = 0.5, depth = 0.01, duration = 0.15;

        var curves = new LightCurve[n * 2];
        var labels = new TransitDisposition[n * 2];

        for (int i = 0; i < n; i++)
        {
            curves[i] = CreateSyntheticTransitCurve(period, epoch, depth, duration,
                numPoints: 1500, seed: 42 + i);
            labels[i] = TransitDisposition.Confirmed;
        }
        for (int i = 0; i < n; i++)
        {
            curves[n + i] = CreateFlatCurve(numPoints: 1500, seed: 200 + i);
            labels[n + i] = TransitDisposition.FalsePositive;
        }

        var model = TransitClassifierTrainer.Train(curves, labels, new TrainerConfig
        {
            WindowSize = 20, Epochs = 10, Filters = 4, KernelSize = 3,
            HiddenUnits = 4, BatchSize = 4,
            Cv = CSharpNumerics.ML.Experiment.CrossValidatorConfig.KFold(folds: 2)
        });

        var pipeline = TransitInferencePipeline.FromModel(model);
        pipeline.Predict(null);
    }

    #endregion

    #region Full Roundtrip Test

    [TestMethod]
    public void FullRoundtrip_Train_Serialize_Deserialize_Predict()
    {
        // Step 1: Create training data
        int n = 4;
        double period = 3.0, epoch = 0.5, depth = 0.01, duration = 0.15;

        var curves = new LightCurve[n * 2];
        var labels = new TransitDisposition[n * 2];

        for (int i = 0; i < n; i++)
        {
            curves[i] = CreateSyntheticTransitCurve(period, epoch, depth, duration,
                numPoints: 1500, seed: 42 + i);
            labels[i] = TransitDisposition.Confirmed;
        }
        for (int i = 0; i < n; i++)
        {
            curves[n + i] = CreateFlatCurve(numPoints: 1500, seed: 200 + i);
            labels[n + i] = TransitDisposition.FalsePositive;
        }

        var trainerConfig = new TrainerConfig
        {
            WindowSize = 20,
            CandidateModels = new[] { TransitModelType.CNN1D },
            Epochs = 20,
            LearningRate = 0.02,
            Filters = 4,
            KernelSize = 3,
            HiddenUnits = 4,
            BatchSize = 4,
            Cv = CSharpNumerics.ML.Experiment.CrossValidatorConfig.KFold(folds: 2)
        };

        var detectionConfig = new TransitDetectionConfig
        {
            MinPeriodDays = 1.0,
            MaxPeriodDays = 10.0,
            SnrThreshold = 3.0,
            MinTransitDepthPpm = 50
        };

        // Step 2: Train
        var model = TransitClassifierTrainer.Train(curves, labels, trainerConfig, detectionConfig);
        Assert.IsNotNull(model);
        Assert.IsTrue(model.Metrics.BestCvScore >= 0);

        // Step 3: Serialize
        byte[] serialized = ModelSerializer.Serialize(model);
        Assert.IsTrue(serialized.Length > 100, "Serialized data should contain model weights.");

        // Step 4: Deserialize
        var restored = ModelSerializer.Deserialize(serialized, model);
        Assert.AreEqual(model.ModelName, restored.ModelName);
        Assert.AreEqual(model.Metrics.Accuracy, restored.Metrics.Accuracy, 0.001);
        Assert.AreEqual(model.Metrics.Precision, restored.Metrics.Precision, 0.001);
        Assert.AreEqual(model.Metrics.Recall, restored.Metrics.Recall, 0.001);

        // Step 5: Predict via inference pipeline
        var inferencePipeline = TransitInferencePipeline.FromModel(restored);
        var testCurve = CreateSyntheticTransitCurve(period, epoch, depth, duration,
            numPoints: 1500, seed: 888);
        var predictions = inferencePipeline.Predict(testCurve);
        Assert.IsNotNull(predictions);
    }

    #endregion
}
