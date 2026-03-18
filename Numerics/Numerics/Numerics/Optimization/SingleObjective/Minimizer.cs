using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;
using System;

namespace CSharpNumerics.Numerics.Optimization.SingleObjective;

/// <summary>
/// Convenience runner that minimises an <see cref="IObjectiveFunction"/>
/// using an <see cref="IOptimizer"/> and an <see cref="IConvergenceCriterion"/>.
/// </summary>
public class Minimizer
{
    public IOptimizer Optimizer { get; }
    public IConvergenceCriterion Convergence { get; }

    public int IterationsUsed { get; private set; }
    public double FinalLoss { get; private set; }

    public Minimizer(IOptimizer optimizer, IConvergenceCriterion convergence)
    {
        Optimizer = optimizer;
        Convergence = convergence;
    }

    /// <summary>
    /// Minimise f starting from x0.
    /// </summary>
    public VectorN Minimize(IObjectiveFunction f, VectorN x0)
    {
        Optimizer.Reset();
        Convergence.Reset();

        VectorN x = new VectorN(x0.Values);
        double loss = 0;

        for (int i = 0; ; i++)
        {
            var grad = f.Gradient(x.Values);
            loss = f.Evaluate(x.Values);
            var gradVec = new VectorN(grad);
            double gradNorm = gradVec.Norm();

            if (Convergence.HasConverged(i, loss, gradNorm))
            {
                IterationsUsed = i;
                FinalLoss = loss;
                return x;
            }

            x = Optimizer.Step(x, gradVec);
        }
    }
}
