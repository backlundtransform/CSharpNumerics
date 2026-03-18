using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Numerics.Optimization.Interfaces;

/// <summary>
/// General-purpose single-objective optimizer.
/// Minimises a scalar function f(x) given its gradient ∇f(x).
/// </summary>
public interface IOptimizer
{
    /// <summary>
    /// Perform a single parameter update step.
    /// </summary>
    /// <param name="parameters">Current parameter vector (modified in-place conceptually; returns updated).</param>
    /// <param name="gradient">Gradient of the objective w.r.t. parameters.</param>
    /// <returns>Updated parameter vector.</returns>
    VectorN Step(VectorN parameters, VectorN gradient);

    /// <summary>
    /// Reset any internal state (e.g. momentum buffers) for a fresh optimisation run.
    /// </summary>
    void Reset();
}
