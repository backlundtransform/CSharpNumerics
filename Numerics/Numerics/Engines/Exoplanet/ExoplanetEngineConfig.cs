using CSharpNumerics.Engines.Exoplanet.Pipeline;

namespace CSharpNumerics.Engines.Exoplanet;

/// <summary>
/// Configuration for <see cref="ExoplanetEngine"/>.
/// Bundles detection and training settings together with an optional pre-trained model path.
/// </summary>
public class ExoplanetEngineConfig
{
    /// <summary>Transit detection configuration (period range, SNR thresholds, detrending).</summary>
    public TransitDetectionConfig Detection { get; set; } = new TransitDetectionConfig();

    /// <summary>ML training configuration (model type, hyperparameters, cross-validation).</summary>
    public TrainerConfig Training { get; set; } = new TrainerConfig();

    /// <summary>
    /// Optional path to a serialized <see cref="TrainedTransitModel"/> file.
    /// When set, the engine loads the model at <see cref="ExoplanetEngine.Init"/> and uses
    /// ML-assisted classification in addition to classical detection.
    /// </summary>
    public string ModelPath { get; set; }
}
