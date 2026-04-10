using System;
using CSharpNumerics.Engines.Exoplanet.Data;

namespace CSharpNumerics.Engines.Exoplanet.Pipeline;

/// <summary>
/// Stateless inference pipeline for transit classification.
/// Designed for use in web services: load a serialized model once, then call
/// <see cref="Predict"/> for each incoming light curve. Thread-safe (no mutable state).
/// </summary>
public class TransitInferencePipeline
{
    private readonly TrainedTransitModel _model;

    private TransitInferencePipeline(TrainedTransitModel model)
    {
        _model = model ?? throw new ArgumentNullException(nameof(model));
    }

    /// <summary>
    /// Loads a trained model from a serialized byte array.
    /// The template provides the pipeline structure (model architecture) for deserialization.
    /// </summary>
    public static TransitInferencePipeline LoadModel(byte[] serialized, TrainedTransitModel template)
    {
        if (serialized == null) throw new ArgumentNullException(nameof(serialized));
        if (template == null) throw new ArgumentNullException(nameof(template));

        var model = ModelSerializer.Deserialize(serialized, template);
        return new TransitInferencePipeline(model);
    }

    /// <summary>
    /// Loads a trained model directly (already deserialized or just trained).
    /// </summary>
    public static TransitInferencePipeline FromModel(TrainedTransitModel model)
    {
        return new TransitInferencePipeline(model);
    }

    /// <summary>
    /// Classifies transit candidates in a light curve.
    /// Thread-safe: this method has no side effects on the pipeline state.
    /// </summary>
    public TransitPrediction[] Predict(LightCurve lc)
    {
        if (lc == null) throw new ArgumentNullException(nameof(lc));
        return _model.Predict(lc);
    }

    /// <summary>
    /// Returns the model name (architecture) used by this pipeline.
    /// </summary>
    public string ModelName => _model.ModelName;

    /// <summary>
    /// Returns the training metrics for the loaded model.
    /// </summary>
    public TrainingMetrics Metrics => _model.Metrics;
}
