using System;
using System.Collections.Generic;
using CSharpNumerics.Engines.Exoplanet.Data;
using CSharpNumerics.Engines.Exoplanet.Enums;
using CSharpNumerics.Engines.Exoplanet.Features;
using CSharpNumerics.ML;
using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.Experiment;
using CSharpNumerics.ML.Sequence.Models.Classification;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Engines.Exoplanet.Pipeline;

/// <summary>
/// Façade that combines feature extraction with <see cref="SupervisedExperiment"/>
/// to train a transit classifier via grid search over CNN1D, LSTM, and BiLSTM architectures.
/// </summary>
public static class TransitClassifierTrainer
{
    /// <summary>
    /// Trains a transit classification model from labeled light curves.
    /// </summary>
    /// <param name="curves">Input light curves.</param>
    /// <param name="labels">Disposition labels (one per light curve). Mapped to binary: Confirmed/Candidate → 1, others → 0.</param>
    /// <param name="trainerConfig">Training hyperparameters and model selection config.</param>
    /// <param name="detectionConfig">Transit detection config (used during inference).</param>
    /// <returns>A trained model ready for prediction.</returns>
    public static TrainedTransitModel Train(
        LightCurve[] curves,
        TransitDisposition[] labels,
        TrainerConfig trainerConfig,
        TransitDetectionConfig detectionConfig = null)
    {
        if (curves == null) throw new ArgumentNullException(nameof(curves));
        if (labels == null) throw new ArgumentNullException(nameof(labels));
        if (curves.Length != labels.Length) throw new ArgumentException("Curves and labels must have the same length.");
        if (curves.Length == 0) throw new ArgumentException("At least one training sample is required.");
        if (trainerConfig == null) throw new ArgumentNullException(nameof(trainerConfig));

        detectionConfig ??= new TransitDetectionConfig();

        // Step 1: Run detection pipeline on each curve to get candidates + features
        var candidateList = new List<TransitCandidate>();
        var curveList = new List<LightCurve>();
        var binaryLabels = new List<double>();

        for (int i = 0; i < curves.Length; i++)
        {
            TransitCandidate candidate;
            try
            {
                var candidates = TransitDetectionPipeline.Detect(curves[i], detectionConfig);
                if (candidates.Length > 0)
                {
                    candidate = candidates[0];
                }
                else
                {
                    // Create a dummy candidate for curves with no detected transit
                    candidate = CreateDummyCandidate(curves[i], detectionConfig);
                }
            }
            catch
            {
                candidate = CreateDummyCandidate(curves[i], detectionConfig);
            }

            // Extract features
            candidate.Features = TransitFeatureExtractor.Extract(candidate, curves[i]);

            candidateList.Add(candidate);
            curveList.Add(curves[i]);
            binaryLabels.Add(IsPositiveLabel(labels[i]) ? 1.0 : 0.0);
        }

        // Step 2: Create windowed training data
        var (X, y) = WindowedFeatureExtractor.CreateTrainingData(
            curveList.ToArray(),
            candidateList.ToArray(),
            binaryLabels.ToArray(),
            trainerConfig.WindowSize,
            trainerConfig.PhaseBins);

        // Step 3: Build pipeline grid
        int numFeatures = TransitFeatureSet.FeatureNames.All.Length;
        int totalFeatures = trainerConfig.WindowSize + numFeatures;

        var grid = new PipelineGrid();

        foreach (var modelType in trainerConfig.CandidateModels)
        {
            switch (modelType)
            {
                case TransitModelType.CNN1D:
                    grid.AddModel<CNN1DClassifier>(g => g
                        .Add("TimeSteps", totalFeatures)
                        .Add("Features", 1)
                        .Add("Filters", trainerConfig.Filters)
                        .Add("KernelSize", trainerConfig.KernelSize)
                        .Add("HiddenUnits", trainerConfig.HiddenUnits)
                        .Add("LearningRate", trainerConfig.LearningRate)
                        .Add("Epochs", trainerConfig.Epochs)
                        .Add("BatchSize", trainerConfig.BatchSize)
                        .Add("ValidationSplit", trainerConfig.ValidationSplit)
                        .Add("UseGlobalAveragePooling", true)
                        .Add("Activation", ActivationType.ReLU));
                    break;

                case TransitModelType.LSTM:
                    grid.AddModel<LSTMClassifier>(g => g
                        .Add("TimeSteps", totalFeatures)
                        .Add("Features", 1)
                        .Add("HiddenSize", trainerConfig.HiddenSize)
                        .Add("HiddenUnits", trainerConfig.HiddenUnits)
                        .Add("LearningRate", trainerConfig.LearningRate)
                        .Add("Epochs", trainerConfig.Epochs)
                        .Add("BatchSize", trainerConfig.BatchSize)
                        .Add("ValidationSplit", trainerConfig.ValidationSplit));
                    break;

                case TransitModelType.BiLSTM:
                    grid.AddModel<BiLSTMClassifier>(g => g
                        .Add("TimeSteps", totalFeatures)
                        .Add("Features", 1)
                        .Add("HiddenSize", trainerConfig.HiddenSize)
                        .Add("HiddenUnits", trainerConfig.HiddenUnits)
                        .Add("LearningRate", trainerConfig.LearningRate)
                        .Add("Epochs", trainerConfig.Epochs)
                        .Add("BatchSize", trainerConfig.BatchSize)
                        .Add("ValidationSplit", trainerConfig.ValidationSplit));
                    break;
            }
        }

        // Step 4: Run experiment
        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(grid)
            .WithCrossValidator(trainerConfig.Cv)
            .Run();

        // Step 5: Extract metrics
        Matrix confusionMatrix = result.BestConfusionMatrix;
        var metrics = new TrainingMetrics(result.BestScore, confusionMatrix, result.BestScore);

        // Step 6: Retrain best pipeline on full dataset
        var bestPipeline = result.BestPipeline.Clone();
        bestPipeline.Fit(X, y);

        return new TrainedTransitModel(bestPipeline, detectionConfig,
            trainerConfig, result.BestModelName, metrics);
    }

    private static bool IsPositiveLabel(TransitDisposition disposition)
    {
        return disposition == TransitDisposition.Confirmed ||
               disposition == TransitDisposition.Candidate;
    }

    private static TransitCandidate CreateDummyCandidate(LightCurve lc, TransitDetectionConfig config)
    {
        double midPeriod = (config.MinPeriodDays + config.MaxPeriodDays) / 2.0;
        double epoch = lc.Length > 0 ? lc.Time[0] : 0.0;
        var parameters = new TransitParameters(midPeriod, epoch, 0.0, 0.0, 0.0, 0.0, 0.0);
        var phaseFolded = LightCurve.FromArrays(
            new double[] { 0.0 }, new double[] { 1.0 });

        return new TransitCandidate(parameters, 0.0, TransitDisposition.Unknown,
            phaseFolded, new TransitFeatureSet());
    }
}
