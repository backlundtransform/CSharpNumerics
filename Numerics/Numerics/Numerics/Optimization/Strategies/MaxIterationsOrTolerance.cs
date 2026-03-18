using CSharpNumerics.Numerics.Optimization.Interfaces;
using System;

namespace CSharpNumerics.Numerics.Optimization.Strategies;

/// <summary>
/// Stops optimisation when a maximum number of iterations is reached
/// or the gradient norm drops below a tolerance.
/// </summary>
public class MaxIterationsOrTolerance : IConvergenceCriterion
{
    public int MaxIterations { get; set; }
    public double Tolerance { get; set; }

    public MaxIterationsOrTolerance(int maxIterations = 1000, double tolerance = 1e-7)
    {
        MaxIterations = maxIterations;
        Tolerance = tolerance;
    }

    public bool HasConverged(int iteration, double currentLoss, double gradientNorm)
    {
        if (iteration >= MaxIterations) return true;
        if (gradientNorm < Tolerance) return true;
        return false;
    }

    public void Reset() { }
}
