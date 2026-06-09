using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Extension methods for quantum-mechanical observables of a 1-D wavefunction sampled on a
/// uniform grid: normalisation, probability density, and expectation values / uncertainties of
/// position and momentum. The wavefunction is represented as a <see cref="ComplexNumber"/> array
/// with a known grid spacing Δx.
/// </summary>
public static class WaveFunctionExtensions
{
    /// <summary>
    /// Squared norm ∫|Ψ|² dx (the total probability). Equals 1 for a normalised state.
    /// </summary>
    public static double NormSquared(this ComplexNumber[] psi, double dx)
    {
        if (psi == null) throw new ArgumentNullException(nameof(psi));
        double sum = 0.0;
        foreach (var c in psi)
        {
            double mag = c.GetMagnitude();
            sum += mag * mag;
        }
        return sum * dx;
    }

    /// <summary>
    /// Returns a normalised copy of the wavefunction so that ∫|Ψ|² dx = 1.
    /// </summary>
    public static ComplexNumber[] Normalize(this ComplexNumber[] psi, double dx)
    {
        double norm = Math.Sqrt(psi.NormSquared(dx));
        if (norm == 0.0) throw new InvalidOperationException("Cannot normalise a zero wavefunction.");

        var result = new ComplexNumber[psi.Length];
        for (int i = 0; i < psi.Length; i++)
            result[i] = psi[i] / norm;
        return result;
    }

    /// <summary>
    /// Probability density |Ψ(x)|² at every grid point.
    /// </summary>
    public static double[] ProbabilityDensity(this ComplexNumber[] psi)
    {
        var density = new double[psi.Length];
        for (int i = 0; i < psi.Length; i++)
        {
            double mag = psi[i].GetMagnitude();
            density[i] = mag * mag;
        }
        return density;
    }

    /// <summary>
    /// Position expectation value ⟨x⟩ = ∫ x|Ψ|² dx / ∫|Ψ|² dx.
    /// </summary>
    public static double ExpectationPosition(this ComplexNumber[] psi, double[] grid, double dx)
    {
        if (grid.Length != psi.Length) throw new ArgumentException("Grid and wavefunction lengths must match.");
        double weighted = 0.0, norm = 0.0;
        for (int i = 0; i < psi.Length; i++)
        {
            double mag = psi[i].GetMagnitude();
            double p = mag * mag;
            weighted += grid[i] * p;
            norm += p;
        }
        return weighted / norm;
    }

    /// <summary>
    /// Position spread (standard deviation) Δx = √(⟨x²⟩ − ⟨x⟩²).
    /// </summary>
    public static double PositionUncertainty(this ComplexNumber[] psi, double[] grid, double dx)
    {
        double mean = psi.ExpectationPosition(grid, dx);
        double weighted = 0.0, norm = 0.0;
        for (int i = 0; i < psi.Length; i++)
        {
            double mag = psi[i].GetMagnitude();
            double p = mag * mag;
            weighted += grid[i] * grid[i] * p;
            norm += p;
        }
        double meanSquare = weighted / norm;
        return Math.Sqrt(Math.Max(meanSquare - (mean * mean), 0.0));
    }

    /// <summary>
    /// Momentum expectation value ⟨p⟩ = ∫ Ψ* (−iħ ∂/∂x) Ψ dx / ∫|Ψ|² dx, using a central
    /// difference for the derivative. Real for any valid state (zero for a real wavefunction).
    /// </summary>
    public static double ExpectationMomentum(this ComplexNumber[] psi, double dx, double hbar = 1.0)
    {
        int n = psi.Length;
        var minusIHbar = new ComplexNumber(0.0, -hbar);

        double accumulated = 0.0;
        for (int i = 0; i < n; i++)
        {
            // Central difference with periodic-free (clamped) boundaries.
            ComplexNumber next = i < n - 1 ? psi[i + 1] : psi[i];
            ComplexNumber prev = i > 0 ? psi[i - 1] : psi[i];
            int span = (i < n - 1 ? 1 : 0) + (i > 0 ? 1 : 0);
            if (span == 0) continue;

            ComplexNumber derivative = (next - prev) / (span * dx);
            ComplexNumber term = psi[i].GetConjugate() * (minusIHbar * derivative);
            accumulated += term.realPart;   // imaginary part integrates to ~0 for valid ψ
        }

        double normSq = psi.NormSquared(dx) / dx;   // Σ|ψ|² (dx cancels with the integral below)
        return (accumulated) / normSq;
    }

    /// <summary>
    /// Lifts a real-valued wavefunction (e.g. a stationary state) to a complex array.
    /// </summary>
    public static ComplexNumber[] ToComplexWaveFunction(this double[] realWaveFunction)
    {
        var result = new ComplexNumber[realWaveFunction.Length];
        for (int i = 0; i < realWaveFunction.Length; i++)
            result[i] = new ComplexNumber(realWaveFunction[i], 0.0);
        return result;
    }
}
