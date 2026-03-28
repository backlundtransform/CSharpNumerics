using CSharpNumerics.Physics.Quantum;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Quantum.Algorithms;

/// <summary>
/// Inverse Quantum Fourier Transform (QFT†).
/// Generates a <see cref="QuantumCircuit"/> that reverses the QFT:
/// SWAP cascade first, then reversed controlled-phase rotations with
/// negated angles, and Hadamard gates in reverse order.
/// </summary>
public static class InverseQFT
{
    /// <summary>
    /// Creates an inverse QFT circuit operating on the given <paramref name="qubits"/>.
    /// </summary>
    /// <param name="totalQubits">Total number of qubits in the circuit.</param>
    /// <param name="qubits">Qubit indices on which to apply QFT† (same order as QFT).</param>
    public static QuantumCircuit CreateCircuit(int totalQubits, params int[] qubits)
    {
        if (qubits == null || qubits.Length == 0)
            throw new ArgumentException("At least one qubit must be specified.", nameof(qubits));

        var circuit = new QuantumCircuit(totalQubits);
        int n = qubits.Length;

        // SWAP cascade to reverse qubit order (done first in inverse)
        for (int i = 0; i < n / 2; i++)
        {
            circuit.AddInstruction(new QuantumInstruction(
                new SWAPGate(), new List<int> { qubits[i], qubits[n - 1 - i] }));
        }

        // Reverse order: j from n-1 down to 0
        for (int j = n - 1; j >= 0; j--)
        {
            // Controlled-phase rotations with negated angles, k from n-1 down to j+1
            for (int k = n - 1; k > j; k--)
            {
                double angle = -Math.PI / (1 << (k - j)); // -π / 2^(k-j)
                circuit.AddInstruction(new QuantumInstruction(
                    new CPhaseGate(angle), new List<int> { qubits[k], qubits[j] }));
            }

            // Hadamard on qubit j
            circuit.AddInstruction(new QuantumInstruction(
                new HadamardGate(), new List<int> { qubits[j] }));
        }

        return circuit;
    }
}
