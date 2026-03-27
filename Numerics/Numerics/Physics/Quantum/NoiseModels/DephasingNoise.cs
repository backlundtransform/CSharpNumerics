using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Quantum.NoiseModels;

/// <summary>
/// Single-qubit dephasing (phase-flip) channel. With probability p the relative
/// phase of the qubit is randomised, modelling T₂ decoherence.
///
/// Kraus operators:
///   E₀ = √(1 − p) · I
///   E₁ = √p · Z
/// </summary>
public class DephasingNoise : INoiseChannel
{
    /// <summary>Dephasing probability ∈ [0, 1].</summary>
    public double Probability { get; }

    public int QubitCount => 1;

    public DephasingNoise(double probability)
    {
        if (probability < 0 || probability > 1)
            throw new ArgumentOutOfRangeException(nameof(probability), "Probability must be in [0, 1].");
        Probability = probability;
    }

    public ComplexMatrix[] GetKrausOperators()
    {
        double p = Probability;
        double s0 = Math.Sqrt(1.0 - p);
        double s1 = Math.Sqrt(p);
        var zero = new ComplexNumber(0, 0);

        // E₀ = √(1-p) · I
        var e0 = new ComplexMatrix(new ComplexNumber[,]
        {
            { new ComplexNumber(s0, 0), zero },
            { zero, new ComplexNumber(s0, 0) }
        });

        // E₁ = √p · Z
        var e1 = new ComplexMatrix(new ComplexNumber[,]
        {
            { new ComplexNumber(s1, 0), zero },
            { zero, new ComplexNumber(-s1, 0) }
        });

        return new[] { e0, e1 };
    }
}
