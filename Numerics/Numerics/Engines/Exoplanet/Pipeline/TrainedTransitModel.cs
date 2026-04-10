using System;
using System.Collections.Generic;
using CSharpNumerics.Engines.Exoplanet.Data;
using CSharpNumerics.Engines.Exoplanet.Enums;
using CSharpNumerics.Engines.Exoplanet.Features;
using CSharpNumerics.Numerics.Objects;
using MLPipeline = CSharpNumerics.ML.Pipeline;

namespace CSharpNumerics.Engines.Exoplanet.Pipeline;

/// <summary>
/// Wraps a trained ML pipeline with the transit detection configuration and training metrics.
/// Provides a high-level <see cref="Predict"/> method for classifying new light curves.
/// </summary>
public class TrainedTransitModel
{
    /// <summary>The trained ML pipeline (model + optional scaler/selector).</summary>
    public MLPipeline Pipeline { get; }

    /// <summary>Detection config used during training (period range, SNR thresholds, etc.).</summary>
    public TransitDetectionConfig DetectionConfig { get; }

    /// <summary>Trainer config used during training (window size, hyperparameters).</summary>
    public TrainerConfig TrainerCfg { get; }

    /// <summary>Name of the best model architecture (e.g. "CNN1DClassifier").</summary>
    public string ModelName { get; }

    /// <summary>Training metrics from cross-validation.</summary>
    public TrainingMetrics Metrics { get; }

    public TrainedTransitModel(MLPipeline pipeline, TransitDetectionConfig detectionConfig,
        TrainerConfig trainerCfg, string modelName, TrainingMetrics metrics)
    {
        Pipeline = pipeline ?? throw new ArgumentNullException(nameof(pipeline));
        DetectionConfig = detectionConfig ?? throw new ArgumentNullException(nameof(detectionConfig));
        TrainerCfg = trainerCfg ?? throw new ArgumentNullException(nameof(trainerCfg));
        ModelName = modelName ?? string.Empty;
        Metrics = metrics ?? new TrainingMetrics();
    }

    /// <summary>
    /// Runs the full prediction pipeline on a light curve:
    /// detect candidates → extract features → classify with ML model.
    /// </summary>
    public TransitPrediction[] Predict(LightCurve lc)
    {
        if (lc == null) throw new ArgumentNullException(nameof(lc));

        // Step 1: Detect candidates using the classical pipeline
        var candidates = TransitDetectionPipeline.Detect(lc, DetectionConfig);
        if (candidates.Length == 0)
            return Array.Empty<TransitPrediction>();

        var predictions = new List<TransitPrediction>();

        foreach (var candidate in candidates)
        {
            // Step 2: Extract features
            candidate.Features = TransitFeatureExtractor.Extract(candidate, lc);

            // Step 3: Create windowed input for ML model
            var X = WindowedFeatureExtractor.CreateInferenceData(
                lc, candidate, TrainerCfg.WindowSize, TrainerCfg.PhaseBins);

            // Step 4: Classify
            var prediction = Pipeline.Predict(X);
            double predictedClass = prediction[0];

            // Binary classification: 0 = no transit, 1 = transit
            bool isTransit = predictedClass >= 0.5;
            double probability = isTransit ? 1.0 : 0.0;

            var disposition = isTransit
                ? TransitDisposition.Candidate
                : TransitDisposition.FalsePositive;

            predictions.Add(new TransitPrediction(candidate, probability, disposition));
        }

        return predictions.ToArray();
    }
}

/// <summary>
/// Training metrics from the transit classifier training process.
/// </summary>
public class TrainingMetrics
{
    public double Accuracy { get; set; }
    public double Precision { get; set; }
    public double Recall { get; set; }
    public double F1Score { get; set; }
    public Matrix ConfusionMatrix { get; set; }
    public double BestCvScore { get; set; }

    public TrainingMetrics() { }

    public TrainingMetrics(double accuracy, Matrix confusionMatrix, double bestCvScore)
    {
        Accuracy = accuracy;
        ConfusionMatrix = confusionMatrix;
        BestCvScore = bestCvScore;
        ComputeFromConfusionMatrix(confusionMatrix);
    }

    private void ComputeFromConfusionMatrix(Matrix cm)
    {
        if (cm.rowLength < 2 || cm.columnLength < 2) return;

        double tp = cm.values[1, 1];
        double fp = cm.values[0, 1];
        double fn = cm.values[1, 0];
        double tn = cm.values[0, 0];

        double total = tp + fp + fn + tn;
        Accuracy = total > 0 ? (tp + tn) / total : 0;
        Precision = (tp + fp) > 0 ? tp / (tp + fp) : 0;
        Recall = (tp + fn) > 0 ? tp / (tp + fn) : 0;
        F1Score = (Precision + Recall) > 0
            ? 2.0 * Precision * Recall / (Precision + Recall)
            : 0;
    }
}
