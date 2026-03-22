using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Abstract base class for quantum gates. A gate is a unitary operator
/// represented by a complex matrix that acts on one or more qubits.
/// </summary>
public abstract class QuantumGate
{
    /// <summary>Number of qubits this gate acts on.</summary>
    public abstract int QubitCount { get; }

    /// <summary>Returns the unitary matrix representation of this gate.</summary>
    public abstract ComplexMatrix GetMatrix();

    /// <summary>
    /// Applies this gate to a quantum state vector.
    /// </summary>
    /// <param name="amplitudes">Full state vector (length 2^totalQubits).</param>
    /// <param name="qubitIndices">Physical qubit indices this gate targets. Length must equal <see cref="QubitCount"/>.</param>
    /// <param name="totalQubits">Total number of qubits in the system.</param>
    /// <returns>New state vector after gate application.</returns>
    public ComplexVectorN Apply(ComplexVectorN amplitudes, int[] qubitIndices, int totalQubits)
    {
        if (qubitIndices.Length != QubitCount)
            throw new ArgumentException($"Gate requires {QubitCount} qubit(s) but received {qubitIndices.Length} index/indices.");

        int stateSize = amplitudes.Length;
        int gateSize = 1 << QubitCount;
        var matrix = GetMatrix();

        var result = new ComplexVectorN(stateSize);
        for (int i = 0; i < stateSize; i++)
            result[i] = new ComplexNumber(amplitudes[i].realPart, amplitudes[i].imaginaryPart);

        int targetMask = 0;
        for (int q = 0; q < QubitCount; q++)
            targetMask |= 1 << qubitIndices[q];

        for (int baseIdx = 0; baseIdx < stateSize; baseIdx++)
        {
            if ((baseIdx & targetMask) != 0) continue;

            var groupIndices = new int[gateSize];
            for (int g = 0; g < gateSize; g++)
            {
                int idx = baseIdx;
                for (int q = 0; q < QubitCount; q++)
                {
                    if ((g & (1 << q)) != 0)
                        idx |= 1 << qubitIndices[q];
                }
                groupIndices[g] = idx;
            }

            var groupAmps = new ComplexVectorN(gateSize);
            for (int g = 0; g < gateSize; g++)
                groupAmps[g] = amplitudes[groupIndices[g]];

            for (int row = 0; row < gateSize; row++)
            {
                var sum = new ComplexNumber(0, 0);
                for (int col = 0; col < gateSize; col++)
                    sum = sum + matrix.values[row, col] * groupAmps[col];
                result[groupIndices[row]] = sum;
            }
        }

        return result;
    }
}
