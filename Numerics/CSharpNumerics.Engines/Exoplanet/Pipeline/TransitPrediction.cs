using CSharpNumerics.Engines.Exoplanet.Data;
using CSharpNumerics.Engines.Exoplanet.Enums;

namespace CSharpNumerics.Engines.Exoplanet.Pipeline;

/// <summary>
/// Prediction result from a trained transit classification model.
/// </summary>
public class TransitPrediction
{
    /// <summary>The transit candidate with parameters and features.</summary>
    public TransitCandidate Candidate { get; }

    /// <summary>Model-estimated probability that this is a genuine transit (0–1).</summary>
    public double Probability { get; }

    /// <summary>Predicted disposition based on the model output.</summary>
    public TransitDisposition PredictedDisposition { get; }

    public TransitPrediction(TransitCandidate candidate, double probability, TransitDisposition predictedDisposition)
    {
        Candidate = candidate;
        Probability = probability;
        PredictedDisposition = predictedDisposition;
    }
}
