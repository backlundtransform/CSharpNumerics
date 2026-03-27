using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Quantum.NoiseModels;

/// <summary>
/// Single-qubit depolarizing channel. With probability p the qubit is replaced
/// by the maximally mixed state (equivalent to applying I, X, Y, Z each with
/// probability p/4), and with probability 1−p it is left unchanged.
///
/// Kraus operators:
///   E₀ = √(1 − 3p/4) · I
///   E₁ = √(p/4) · X
///   E₂ = √(p/4) · Y
///   E₃ = √(p/4) · Z
/// </summary>
public class DepolarizingNoise : INoiseChannel
{
    /// <summary>Depolarizing probability ∈ [0, 1].</summary>
    public double Probability { get; }

    public int QubitCount => 1;

    public DepolarizingNoise(double probability)
    {
        if (probability < 0 || probability > 1)
            throw new ArgumentOutOfRangeException(nameof(probability), "Probability must be in [0, 1].");
        Probability = probability;
    }

    public ComplexMatrix[] GetKrausOperators()
    {
        double p = Probability;
        double s0 = Math.Sqrt(1.0 - 3.0 * p / 4.0);
        double sp = Math.Sqrt(p / 4.0);

        var zero = new ComplexNumber(0, 0);

        // E₀ = √(1 - 3p/4) · I
        var e0 = new ComplexMatrix(new ComplexNumber[,]
        {
            { new ComplexNumber(s0, 0), zero },
            { zero, new ComplexNumber(s0, 0) }
        });

        // E₁ = √(p/4) · X
        var e1 = new ComplexMatrix(new ComplexNumber[,]
        {
            { zero, new ComplexNumber(sp, 0) },
            { new ComplexNumber(sp, 0), zero }
        });

        // E₂ = √(p/4) · Y
        var e2 = new ComplexMatrix(new ComplexNumber[,]
        {
            { zero, new ComplexNumber(0, -sp) },
            { new ComplexNumber(0, sp), zero }
        });

        // E₃ = √(p/4) · Z
        var e3 = new ComplexMatrix(new ComplexNumber[,]
        {
            { new ComplexNumber(sp, 0), zero },
            { zero, new ComplexNumber(-sp, 0) }
        });

        return new[] { e0, e1, e2, e3 };
    }
}
