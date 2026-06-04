using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.ML.Training;

/// <summary>
/// A differentiable loss term over a prediction vector. Implementations return both the
/// scalar penalty and its gradient with respect to the prediction so they can be combined
/// (see <see cref="CompositeLoss"/>) and fed into a training loop.
/// <para>
/// The meaning of <c>target</c> is loss-specific: data-fidelity and conservation losses use
/// it as the reference component vector, while shape penalties such as
/// <see cref="NonNegativityLoss"/> and <see cref="SmoothnessLoss"/> ignore it.
/// </para>
/// </summary>
public interface ILoss
{
    /// <summary>Scalar loss value for the prediction.</summary>
    double Compute(VectorN prediction, VectorN target);

    /// <summary>Gradient of the loss with respect to the prediction.</summary>
    VectorN Gradient(VectorN prediction, VectorN target);
}
