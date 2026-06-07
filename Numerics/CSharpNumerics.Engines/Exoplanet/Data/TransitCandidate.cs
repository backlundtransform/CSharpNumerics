using CSharpNumerics.Engines.Exoplanet.Enums;

namespace CSharpNumerics.Engines.Exoplanet.Data;

public class TransitCandidate
{
    public TransitParameters Parameters { get; set; }
    public double Score { get; set; }
    public TransitDisposition Disposition { get; set; }
    public LightCurve PhaseFoldedCurve { get; set; }
    public TransitFeatureSet Features { get; set; }

    public TransitCandidate() { }

    public TransitCandidate(TransitParameters parameters, double score,
        TransitDisposition disposition, LightCurve phaseFoldedCurve, TransitFeatureSet features)
    {
        Parameters = parameters;
        Score = score;
        Disposition = disposition;
        PhaseFoldedCurve = phaseFoldedCurve;
        Features = features;
    }
}
