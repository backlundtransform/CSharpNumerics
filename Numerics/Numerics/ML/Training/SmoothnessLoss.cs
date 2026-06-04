using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.ML.Training;

/// <summary>
/// Penalises curvature via the discrete second difference: <c>λ Σ (y_{t-1} − 2y_t + y_{t+1})²</c>.
/// Encourages a smooth signal (bounded second derivative), e.g. a slowly varying SRC component.
/// Operates on the prediction vector interpreted as an ordered series; the target is ignored.
/// </summary>
public class SmoothnessLoss : ILoss
{
    private readonly double _weight;

    public SmoothnessLoss(double weight = 1.0)
    {
        if (weight < 0.0) throw new ArgumentOutOfRangeException(nameof(weight), "Weight must be non-negative.");
        _weight = weight;
    }

    public double Compute(VectorN prediction, VectorN target)
    {
        int n = prediction.Length;
        if (n < 3) return 0.0;

        double sum = 0.0;
        for (int t = 1; t < n - 1; t++)
        {
            double d = prediction[t - 1] - (2.0 * prediction[t]) + prediction[t + 1];
            sum += d * d;
        }
        return _weight * sum;
    }

    public VectorN Gradient(VectorN prediction, VectorN target)
    {
        int n = prediction.Length;
        var grad = new double[n];
        if (n < 3) return new VectorN(grad);

        for (int t = 1; t < n - 1; t++)
        {
            double d = prediction[t - 1] - (2.0 * prediction[t]) + prediction[t + 1];
            double scaled = 2.0 * _weight * d;
            // Scatter through ∂d_t/∂y: +1 at t-1, -2 at t, +1 at t+1.
            grad[t - 1] += scaled;
            grad[t] += -2.0 * scaled;
            grad[t + 1] += scaled;
        }
        return new VectorN(grad);
    }
}
