using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// The quantum Fourier transform on n qubits:
/// QFT|x⟩ = 2^(−n/2) Σ_y e^(2πi·xy/2ⁿ) |y⟩, and its inverse with e^(−2πi·xy/2ⁿ).
///
/// The register is read as an integer with qubitIndices[0] as the least significant bit,
/// the same convention as <see cref="ModularMultiplyGate"/>.
///
/// <see cref="QuantumGate.GetMatrix"/> builds the dense 2ⁿ×2ⁿ matrix, which is fine for a few
/// qubits. <see cref="ApplyAsCircuit"/> applies the same transform as n Hadamard gates,
/// n(n − 1)/2 controlled-phase gates and ⌊n/2⌋ swaps. That is the textbook circuit, and it is
/// much faster on large registers.
/// </summary>
public class QuantumFourierTransform : QuantumGate
{
    private readonly int _qubitCount;

    /// <summary>True for the inverse transform QFT†.</summary>
    public bool Inverse { get; }

    /// <summary>Creates the QFT (or its inverse) on <paramref name="qubitCount"/> qubits.</summary>
    /// <param name="qubitCount">Number of qubits n ≥ 1.</param>
    /// <param name="inverse">True for the inverse transform.</param>
    public QuantumFourierTransform(int qubitCount, bool inverse = false)
    {
        if (qubitCount < 1 || qubitCount > 30)
            throw new ArgumentOutOfRangeException(nameof(qubitCount), qubitCount, "qubitCount must be between 1 and 30.");
        _qubitCount = qubitCount;
        Inverse = inverse;
    }

    public override int QubitCount => _qubitCount;

    public override ComplexMatrix GetMatrix()
    {
        int size = 1 << _qubitCount;
        double scale = 1.0 / Math.Sqrt(size);
        double sign = Inverse ? -1.0 : 1.0;

        var m = new ComplexNumber[size, size];
        for (int y = 0; y < size; y++)
        {
            for (int x = 0; x < size; x++)
            {
                // Reduce x·y mod 2ⁿ before converting to an angle to keep the phase accurate.
                long k = (long)x * y % size;
                m[y, x] = ComplexNumber.FromPolarCoordinates(scale, sign * 2.0 * Math.PI * k / size);
            }
        }
        return new ComplexMatrix(m);
    }

    /// <summary>
    /// Applies the transform with its gate decomposition instead of the dense matrix.
    /// The result equals <see cref="QuantumGate.Apply"/> up to rounding.
    /// </summary>
    /// <param name="amplitudes">Full state vector (length 2^totalQubits).</param>
    /// <param name="qubitIndices">The register, least significant qubit first. Length must equal <see cref="QubitCount"/>.</param>
    /// <param name="totalQubits">Total number of qubits in the system.</param>
    /// <returns>New state vector after the transform.</returns>
    public ComplexVectorN ApplyAsCircuit(ComplexVectorN amplitudes, int[] qubitIndices, int totalQubits)
    {
        if (qubitIndices == null) throw new ArgumentNullException(nameof(qubitIndices));
        if (qubitIndices.Length != _qubitCount)
            throw new ArgumentException($"The transform acts on {_qubitCount} qubit(s) but received {qubitIndices.Length} index/indices.");

        int n = _qubitCount;
        var hadamard = new HadamardGate();
        var swap = new SWAPGate();
        var state = amplitudes;

        if (!Inverse)
        {
            for (int j = n - 1; j >= 0; j--)
            {
                state = hadamard.Apply(state, new[] { qubitIndices[j] }, totalQubits);
                for (int k = j - 1; k >= 0; k--)
                    state = new CPhaseGate(Math.PI / (1L << (j - k)))
                        .Apply(state, new[] { qubitIndices[k], qubitIndices[j] }, totalQubits);
            }
            for (int i = 0; i < n / 2; i++)
                state = swap.Apply(state, new[] { qubitIndices[i], qubitIndices[n - 1 - i] }, totalQubits);
        }
        else
        {
            // The forward circuit in reverse order, with every phase negated.
            for (int i = 0; i < n / 2; i++)
                state = swap.Apply(state, new[] { qubitIndices[i], qubitIndices[n - 1 - i] }, totalQubits);
            for (int j = 0; j < n; j++)
            {
                for (int k = 0; k < j; k++)
                    state = new CPhaseGate(-Math.PI / (1L << (j - k)))
                        .Apply(state, new[] { qubitIndices[k], qubitIndices[j] }, totalQubits);
                state = hadamard.Apply(state, new[] { qubitIndices[j] }, totalQubits);
            }
        }

        return state;
    }
}
