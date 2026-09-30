using System;

namespace CSharpNumerics.Numerics.RootFinding;

/// <summary>
/// Extension-method facade over <see cref="RootFinder"/>, in the same style as the other
/// numerics extensions.
/// </summary>
public static class RootFindingExtensions
{
    /// <summary>
    /// Finds a root of the function in the bracketing interval [<paramref name="lower"/>,
    /// <paramref name="upper"/>] using Brent's method.
    /// </summary>
    /// <param name="function">The function whose root is sought.</param>
    /// <param name="lower">Lower end of a bracketing interval.</param>
    /// <param name="upper">Upper end of a bracketing interval.</param>
    /// <param name="tolerance">Convergence tolerance on the root position.</param>
    /// <param name="maxIterations">Maximum number of iterations.</param>
    /// <returns>The outcome, including whether the run converged.</returns>
    /// <example>
    /// <code>
    /// Func&lt;double, double&gt; f = x =&gt; x * x - 4;
    /// var root = f.FindRoot(0, 5);   // root.Value ≈ 2
    /// </code>
    /// </example>
    public static RootResult FindRoot(
        this Func<double, double> function,
        double lower,
        double upper,
        double tolerance = RootFinder.DefaultTolerance,
        int maxIterations = RootFinder.DefaultMaxIterations) =>
        RootFinder.Brent(function, lower, upper, tolerance, maxIterations);
}
