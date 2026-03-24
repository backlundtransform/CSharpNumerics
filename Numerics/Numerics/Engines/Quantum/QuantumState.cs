using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Quantum;
using System;

namespace CSharpNumerics.Engines.Quantum;

/// <summary>
/// Represents the quantum state of a multi-qubit system as a vector of complex amplitudes.
/// For n qubits the state has 2^n amplitudes in the computational basis.
/// </summary>
public class QuantumState
{
    /// <summary>Complex amplitudes in the computational basis.</summary>
    public ComplexVectorN Amplitudes { get; }

    /// <summary>Number of qubits in the system.</summary>
    public int QubitCount { get; }

    public QuantumState(ComplexVectorN amplitudes)
    {
        int n = amplitudes.Length;
        if ((n & (n - 1)) != 0)
            throw new ArgumentException("Amplitude count must be a power of 2.", nameof(amplitudes));

        Amplitudes = amplitudes;
        QubitCount = BitLength(n) - 1;
    }

    /// <summary>
    /// Returns the measurement probability of a specific computational basis state.
    /// </summary>
    public double GetProbability(int basisStateIndex)
    {
        if (basisStateIndex < 0 || basisStateIndex >= Amplitudes.Length)
            throw new ArgumentOutOfRangeException(nameof(basisStateIndex));

        var a = Amplitudes[basisStateIndex];
        return a.realPart * a.realPart + a.imaginaryPart * a.imaginaryPart;
    }

    /// <summary>
    /// Returns measurement probabilities for all computational basis states as a VectorN.
    /// </summary>
    public VectorN GetProbabilities()
    {
        var probs = new double[Amplitudes.Length];
        for (int i = 0; i < Amplitudes.Length; i++)
        {
            var a = Amplitudes[i];
            probs[i] = a.realPart * a.realPart + a.imaginaryPart * a.imaginaryPart;
        }
        return new VectorN(probs);
    }

    /// <summary>
    /// Returns the Bloch vector for a single-qubit state.
    /// Throws if the state has more than one qubit.
    /// </summary>
    public BlochVector GetBlochVector()
    {
        if (QubitCount != 1)
            throw new InvalidOperationException(
                $"Bloch vector is only defined for single-qubit states (this state has {QubitCount} qubits).");

        return BlochVector.FromAmplitudes(Amplitudes[0], Amplitudes[1]);
    }

    private static int BitLength(int value)
    {
        int bits = 0;
        while (value > 0) { value >>= 1; bits++; }
        return bits;
    }
}
