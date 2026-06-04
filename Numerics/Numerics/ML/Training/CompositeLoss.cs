using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.Training;

/// <summary>
/// A weighted sum of loss terms: <c>Σ wₖ Lₖ</c>. Both the value and the gradient are the
/// weighted combination of the constituents, so physics constraints (non-negativity,
/// conservation, smoothness) can be layered on top of a data-fidelity term. Build it fluently
/// with <see cref="Add"/>.
/// </summary>
public class CompositeLoss : ILoss
{
    private readonly List<(ILoss Loss, double Weight)> _terms = new();

    /// <summary>
    /// Adds a weighted loss term. Returns <c>this</c> for fluent chaining.
    /// </summary>
    public CompositeLoss Add(ILoss loss, double weight = 1.0)
    {
        if (loss == null) throw new ArgumentNullException(nameof(loss));
        _terms.Add((loss, weight));
        return this;
    }

    /// <summary>Number of loss terms.</summary>
    public int TermCount => _terms.Count;

    public double Compute(VectorN prediction, VectorN target)
    {
        double total = 0.0;
        foreach (var (loss, weight) in _terms)
            total += weight * loss.Compute(prediction, target);
        return total;
    }

    public VectorN Gradient(VectorN prediction, VectorN target)
    {
        var sum = new double[prediction.Length];
        foreach (var (loss, weight) in _terms)
        {
            var g = loss.Gradient(prediction, target);
            for (int i = 0; i < sum.Length; i++)
                sum[i] += weight * g[i];
        }
        return new VectorN(sum);
    }
}
