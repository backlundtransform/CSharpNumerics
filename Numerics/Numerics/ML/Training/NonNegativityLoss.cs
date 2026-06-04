using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.ML.Training;

/// <summary>
/// Penalises negative predictions: <c>λ Σ max(-yᵢ, 0)²</c>. Zero when every component is
/// non-negative, growing quadratically with the magnitude of any negative entry. Used to push
/// flow-decomposition components towards physically valid (non-negative) values. The target is
/// ignored.
/// </summary>
public class NonNegativityLoss : ILoss
{
    private readonly double _weight;

    public NonNegativityLoss(double weight = 1.0)
    {
        if (weight < 0.0) throw new ArgumentOutOfRangeException(nameof(weight), "Weight must be non-negative.");
        _weight = weight;
    }

    public double Compute(VectorN prediction, VectorN target)
    {
        double sum = 0.0;
        for (int i = 0; i < prediction.Length; i++)
        {
            if (prediction[i] < 0.0)
                sum += prediction[i] * prediction[i];
        }
        return _weight * sum;
    }

    public VectorN Gradient(VectorN prediction, VectorN target)
    {
        var grad = new double[prediction.Length];
        for (int i = 0; i < prediction.Length; i++)
        {
            // d/dyᵢ [max(-yᵢ,0)²] = 2yᵢ when yᵢ < 0, else 0.
            grad[i] = prediction[i] < 0.0 ? 2.0 * _weight * prediction[i] : 0.0;
        }
        return new VectorN(grad);
    }
}
