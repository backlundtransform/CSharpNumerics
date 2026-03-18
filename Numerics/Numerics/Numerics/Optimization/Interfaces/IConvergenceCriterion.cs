namespace CSharpNumerics.Numerics.Optimization.Interfaces;

/// <summary>
/// Evaluates whether an iterative optimisation has converged.
/// </summary>
public interface IConvergenceCriterion
{
    /// <summary>
    /// Returns true when the optimisation should stop.
    /// </summary>
    /// <param name="iteration">Current iteration (0-based).</param>
    /// <param name="currentLoss">Current objective value.</param>
    /// <param name="gradientNorm">Norm of the current gradient (optional context).</param>
    bool HasConverged(int iteration, double currentLoss, double gradientNorm);

    /// <summary>
    /// Reset internal state for a new optimisation run.
    /// </summary>
    void Reset();
}
