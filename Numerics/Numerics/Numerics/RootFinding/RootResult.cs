using System;

namespace CSharpNumerics.Numerics.RootFinding;

/// <summary>
/// Outcome of a root-finding run: the approximate root together with the information
/// needed to judge whether it can be trusted.
/// </summary>
/// <remarks>
/// A root finder that stops without converging still returns its best estimate, so the
/// value alone says nothing about success — always check <see cref="Converged"/> before
/// using it. Use <see cref="Value"/> for the raw estimate or <see cref="EnsureConverged"/>
/// to turn a failure into an exception.
/// </remarks>
public sealed class RootResult
{
    public RootResult(double value, bool converged, int iterations, double residual, string method)
    {
        Value = value;
        Converged = converged;
        Iterations = iterations;
        Residual = residual;
        Method = method;
    }

    /// <summary>The approximate root. Only meaningful when <see cref="Converged"/> is true.</summary>
    public double Value { get; }

    /// <summary>True if the method reached its tolerance before exhausting its iteration budget.</summary>
    public bool Converged { get; }

    /// <summary>Number of iterations actually performed.</summary>
    public int Iterations { get; }

    /// <summary>The residual |f(<see cref="Value"/>)| at the returned point.</summary>
    public double Residual { get; }

    /// <summary>Name of the method that produced the result, e.g. "Brent".</summary>
    public string Method { get; }

    /// <summary>
    /// Returns <see cref="Value"/> if the run converged, otherwise throws.
    /// Use when a non-converged result is a programming error rather than something to handle.
    /// </summary>
    /// <exception cref="InvalidOperationException">The run did not converge.</exception>
    public double EnsureConverged()
    {
        if (!Converged)
        {
            throw new InvalidOperationException(
                $"{Method} did not converge after {Iterations} iterations (residual {Residual:G6}).");
        }

        return Value;
    }

    public override string ToString() =>
        Converged
            ? $"{Method}: x = {Value:G10} after {Iterations} iterations (residual {Residual:G6})"
            : $"{Method}: did not converge after {Iterations} iterations (best x = {Value:G10}, residual {Residual:G6})";
}
