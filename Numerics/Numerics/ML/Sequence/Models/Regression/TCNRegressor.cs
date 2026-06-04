using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.ML.Sequence.Models.Internal;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.ML.Sequence.Models.Regression;

/// <summary>
/// Temporal Convolutional Network regressor: <c>TCNBlock → GlobalAvgPool → Dense → output</c>.
/// Suitable for long-range sequence regression (e.g. residual flow separation).
/// </summary>
public class TCNRegressor : TCNModelBase, IRegressionModel
{
    protected override int ResolveOutputSize(VectorN y)
    {
        return 1;
    }

    protected override VectorN ComputeLossGradient(VectorN output, double target)
    {
        return new VectorN(new[] { output[0] - target });
    }

    protected override double ComputeValidationLoss(VectorN output, double target)
    {
        double error = output[0] - target;
        return error * error;
    }

    protected override double MapPrediction(VectorN output)
    {
        return output[0];
    }

    public override IModel Clone()
    {
        var clone = new TCNRegressor();
        CopySharedParametersTo(clone);
        return clone;
    }
}
