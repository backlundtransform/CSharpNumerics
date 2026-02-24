using System;

namespace CSharpNumerics.Numerics;

/// <summary>
/// Provides extension methods for computing one-sided and two-sided limits of functions.
/// </summary>
public static class LimitExtensions
{
    /// <summary>
    /// Approximates the left-hand limit of a function at x by evaluating f(x - ε).
    /// </summary>
    /// <param name="function">The function to evaluate.</param>
    /// <param name="x">Point at which to evaluate the left-hand limit.</param>
    /// <returns>An approximation of the left-hand limit at <paramref name="x"/>.</returns>
    public static double LeftLimit(this Func<double, double> function, double x)
    {
        return function(x - double.Epsilon);
    }

    /// <summary>
    /// Approximates the right-hand limit of a function at x by evaluating f(x + ε).
    /// </summary>
    /// <param name="function">The function to evaluate.</param>
    /// <param name="x">Point at which to evaluate the right-hand limit.</param>
    /// <returns>An approximation of the right-hand limit at <paramref name="x"/>.</returns>
    public static double RightLimit(this Func<double, double> function, double x)
    {
        return function(x + double.Epsilon);
    }

    /// <summary>
    /// Approximates the (two-sided) limit of a function at x by comparing left and right evaluations.
    /// If the left and right evaluations are equal, returns that value; otherwise returns NaN.
    /// </summary>
    /// <param name="function">The function to evaluate.</param>
    /// <param name="x">Point at which to evaluate the limit.</param>
    /// <returns>The limit if the left and right evaluations match; otherwise <see cref="double.NaN"/>.</returns>
    public static double Limit(this Func<double, double> function, double x)
    {
        var right = function.LeftLimit(x);

        var left = function.RightLimit(x);

        return (right == left) ? right : double.NaN;
    }
}
