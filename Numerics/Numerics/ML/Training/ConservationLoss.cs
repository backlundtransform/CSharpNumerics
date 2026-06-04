using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.ML.Training;

/// <summary>
/// Enforces a conservation law on the decomposition: the predicted components must sum to the
/// same total as the target. Penalty <c>λ (Σtargetᵢ − Σpredᵢ)²</c>, i.e. the predicted total
/// <c>Σpred</c> is pulled towards the reference total <c>Σtarget</c> (the measured <c>Q_total</c>).
/// Zero when the sums match.
/// </summary>
public class ConservationLoss : ILoss
{
    private readonly double _weight;

    public ConservationLoss(double weight = 1.0)
    {
        if (weight < 0.0) throw new ArgumentOutOfRangeException(nameof(weight), "Weight must be non-negative.");
        _weight = weight;
    }

    public double Compute(VectorN prediction, VectorN target)
    {
        double diff = Sum(prediction) - Sum(target);
        return _weight * diff * diff;
    }

    public VectorN Gradient(VectorN prediction, VectorN target)
    {
        double diff = Sum(prediction) - Sum(target);
        double g = 2.0 * _weight * diff;     // ∂/∂yᵢ (Σy − Σt)² = 2(Σy − Σt), identical for every i

        var grad = new double[prediction.Length];
        for (int i = 0; i < prediction.Length; i++)
            grad[i] = g;
        return new VectorN(grad);
    }

    private static double Sum(VectorN v)
    {
        double s = 0.0;
        for (int i = 0; i < v.Length; i++) s += v[i];
        return s;
    }
}
