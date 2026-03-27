using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Quantum.NoiseModels;

/// <summary>
/// Single-qubit amplitude damping channel, modelling energy dissipation (T₁ decay).
/// The excited state |1⟩ decays to the ground state |0⟩ with probability γ.
///
/// Kraus operators:
///   E₀ = [[1, 0], [0, √(1−γ)]]
///   E₁ = [[0, √γ], [0, 0]]
/// </summary>
public class AmplitudeDampingNoise : INoiseChannel
{
    /// <summary>Damping parameter γ ∈ [0, 1]. γ = 1 means complete decay.</summary>
    public double Gamma { get; }

    public int QubitCount => 1;

    public AmplitudeDampingNoise(double gamma)
    {
        if (gamma < 0 || gamma > 1)
            throw new ArgumentOutOfRangeException(nameof(gamma), "Gamma must be in [0, 1].");
        Gamma = gamma;
    }

    public ComplexMatrix[] GetKrausOperators()
    {
        double g = Gamma;
        var zero = new ComplexNumber(0, 0);
        var one = new ComplexNumber(1, 0);

        // E₀ = [[1, 0], [0, √(1−γ)]]
        var e0 = new ComplexMatrix(new ComplexNumber[,]
        {
            { one, zero },
            { zero, new ComplexNumber(Math.Sqrt(1.0 - g), 0) }
        });

        // E₁ = [[0, √γ], [0, 0]]
        var e1 = new ComplexMatrix(new ComplexNumber[,]
        {
            { zero, new ComplexNumber(Math.Sqrt(g), 0) },
            { zero, zero }
        });

        return new[] { e0, e1 };
    }
}
