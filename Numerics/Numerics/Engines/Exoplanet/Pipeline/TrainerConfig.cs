using CSharpNumerics.ML.Experiment;

namespace CSharpNumerics.Engines.Exoplanet.Pipeline;

/// <summary>
/// Specifies which sequence model architecture to try during transit classifier training.
/// </summary>
public enum TransitModelType
{
    CNN1D,
    LSTM,
    BiLSTM
}

/// <summary>
/// Configuration for <see cref="TransitClassifierTrainer"/>.
/// Controls windowing, model selection, cross-validation, and training hyperparameters.
/// </summary>
public class TrainerConfig
{
    /// <summary>Number of phase-folded bins per ML input window.</summary>
    public int WindowSize { get; set; } = 50;

    /// <summary>Stride between consecutive windows (used when generating training data from raw arrays).</summary>
    public int Stride { get; set; } = 1;

    /// <summary>Which model architectures to evaluate during grid search.</summary>
    public TransitModelType[] CandidateModels { get; set; } = new[]
    {
        TransitModelType.CNN1D
    };

    /// <summary>Cross-validation configuration for model selection.</summary>
    public CrossValidatorConfig Cv { get; set; } = CrossValidatorConfig.KFold(folds: 3);

    /// <summary>Maximum training epochs per model.</summary>
    public int Epochs { get; set; } = 150;

    /// <summary>Fraction of training data held out for early stopping validation.</summary>
    public double ValidationSplit { get; set; } = 0.0;

    /// <summary>Learning rate for gradient descent.</summary>
    public double LearningRate { get; set; } = 0.02;

    /// <summary>Mini-batch size for training.</summary>
    public int BatchSize { get; set; } = 16;

    /// <summary>Number of convolutional filters (CNN1D only).</summary>
    public int Filters { get; set; } = 8;

    /// <summary>Kernel size for 1D convolutions (CNN1D only).</summary>
    public int KernelSize { get; set; } = 5;

    /// <summary>Hidden units in the dense layer.</summary>
    public int HiddenUnits { get; set; } = 8;

    /// <summary>LSTM/BiLSTM hidden state size.</summary>
    public int HiddenSize { get; set; } = 16;

    /// <summary>Number of phase bins for phase folding.</summary>
    public int PhaseBins { get; set; } = 200;
}
