using CSharpNumerics.Engines.Exoplanet.Data;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Fitting;

namespace CSharpNumerics.Engines.Exoplanet.Pipeline;

public class TransitFitResult
{
    public TransitParameters Parameters { get; }
    public VectorN Uncertainties { get; }
    public double ChiSquared { get; }
    public double BIC { get; }
    public double ReducedChiSquared { get; }
    public VectorN FittedFlux { get; }
    public VectorN Residuals { get; }

    public TransitFitResult(TransitParameters parameters, VectorN uncertainties,
        double chiSquared, double bic, double reducedChiSquared,
        VectorN fittedFlux, VectorN residuals)
    {
        Parameters = parameters;
        Uncertainties = uncertainties;
        ChiSquared = chiSquared;
        BIC = bic;
        ReducedChiSquared = reducedChiSquared;
        FittedFlux = fittedFlux;
        Residuals = residuals;
    }
}
