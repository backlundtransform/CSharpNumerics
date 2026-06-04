namespace CSharpNumerics.ML.Sequence.Enums;

public enum ConvolutionPaddingMode
{
    Valid,
    Same,

    /// <summary>
    /// Causal padding: the input is left-padded by (kernelSize - 1) × dilation so that
    /// output[t] depends only on inputs at or before t. Required for temporal models (TCN)
    /// where leakage from future timesteps must be avoided.
    /// </summary>
    Causal
}