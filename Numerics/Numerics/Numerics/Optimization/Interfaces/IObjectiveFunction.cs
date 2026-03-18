namespace CSharpNumerics.Numerics.Optimization.Interfaces;

/// <summary>
/// A scalar objective function and its gradient, for use with <see cref="IOptimizer"/>.
/// </summary>
public interface IObjectiveFunction
{
    /// <summary>Number of parameters.</summary>
    int Dimension { get; }

    /// <summary>Evaluate the objective value at x.</summary>
    double Evaluate(double[] x);

    /// <summary>Evaluate the gradient at x.</summary>
    double[] Gradient(double[] x);
}
