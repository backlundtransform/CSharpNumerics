using System;

namespace CSharpNumerics.Statistics.Robust;

/// <summary>
/// Huber loss function: a robust loss that is quadratic for small residuals
/// and linear for large residuals.
///
///   L(a) = ½ a²            if |a| ≤ δ
///   L(a) = δ(|a| − ½δ)    if |a| &gt; δ
///
/// where δ is the transition threshold. This provides a smooth blend between
/// MSE (sensitive to outliers) and MAE (robust to outliers).
/// </summary>
public static class HuberLoss
{
    /// <summary>Default transition threshold.</summary>
    public const double DefaultDelta = 1.35;

    /// <summary>
    /// Computes the Huber loss for a single residual.
    /// </summary>
    /// <param name="residual">The residual value a = y − ŷ.</param>
    /// <param name="delta">Transition threshold between quadratic and linear regions.</param>
    /// <returns>The Huber loss value.</returns>
    public static double Loss(double residual, double delta = DefaultDelta)
    {
        if (delta <= 0) throw new ArgumentOutOfRangeException(nameof(delta), "Must be positive.");

        double abs = Math.Abs(residual);
        if (abs <= delta)
            return 0.5 * residual * residual;
        return delta * (abs - 0.5 * delta);
    }

    /// <summary>
    /// Computes the derivative (gradient) of the Huber loss for a single residual.
    /// </summary>
    /// <param name="residual">The residual value a = y − ŷ.</param>
    /// <param name="delta">Transition threshold.</param>
    /// <returns>The gradient of the Huber loss.</returns>
    public static double Gradient(double residual, double delta = DefaultDelta)
    {
        if (delta <= 0) throw new ArgumentOutOfRangeException(nameof(delta), "Must be positive.");

        if (Math.Abs(residual) <= delta)
            return residual;
        return delta * Math.Sign(residual);
    }

    /// <summary>
    /// Computes the mean Huber loss over an array of residuals.
    /// </summary>
    /// <param name="residuals">Array of residual values.</param>
    /// <param name="delta">Transition threshold.</param>
    /// <returns>The mean Huber loss.</returns>
    public static double MeanLoss(double[] residuals, double delta = DefaultDelta)
    {
        if (residuals == null) throw new ArgumentNullException(nameof(residuals));
        if (residuals.Length == 0) throw new ArgumentException("Residuals must not be empty.", nameof(residuals));
        if (delta <= 0) throw new ArgumentOutOfRangeException(nameof(delta), "Must be positive.");

        double sum = 0;
        for (int i = 0; i < residuals.Length; i++)
            sum += Loss(residuals[i], delta);

        return sum / residuals.Length;
    }

    /// <summary>
    /// Computes the Huber weight for a residual, suitable for IRLS:
    /// w = 1 if |a| ≤ δ, w = δ / |a| otherwise.
    /// </summary>
    /// <param name="residual">The residual value.</param>
    /// <param name="delta">Transition threshold.</param>
    /// <returns>The Huber weight in (0, 1].</returns>
    public static double Weight(double residual, double delta = DefaultDelta)
    {
        if (delta <= 0) throw new ArgumentOutOfRangeException(nameof(delta), "Must be positive.");

        double abs = Math.Abs(residual);
        if (abs <= delta)
            return 1.0;
        return delta / abs;
    }
}
