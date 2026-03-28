using CSharpNumerics.Physics.Quantum;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Quantum.Algorithms;

/// <summary>
/// Quantum Fourier Transform (QFT).
/// Generates a <see cref="QuantumCircuit"/> implementing the standard QFT
/// on the specified qubits: H gates, controlled-phase rotations, and a
/// final SWAP cascade to reverse qubit order.
/// </summary>
public static class QFT
{
    /// <summary>
    /// Creates a QFT circuit operating on the given <paramref name="qubits"/>.
    /// </summary>
    /// <param name="totalQubits">Total number of qubits in the circuit.</param>
    /// <param name="qubits">Qubit indices on which to apply QFT (order matters).</param>
    public static QuantumCircuit CreateCircuit(int totalQubits, params int[] qubits)
    {
        if (qubits == null || qubits.Length == 0)
            throw new ArgumentException("At least one qubit must be specified.", nameof(qubits));

        var circuit = new QuantumCircuit(totalQubits);
        int n = qubits.Length;

        for (int j = 0; j < n; j++)
        {
            // Hadamard on qubit j
            circuit.AddInstruction(new QuantumInstruction(
                new HadamardGate(), new List<int> { qubits[j] }));

            // Controlled-phase rotations
            for (int k = j + 1; k < n; k++)
            {
                double angle = Math.PI / (1 << (k - j)); // π / 2^(k-j)
                circuit.AddInstruction(new QuantumInstruction(
                    new CPhaseGate(angle), new List<int> { qubits[k], qubits[j] }));
            }
        }

        // SWAP cascade to reverse qubit order
        for (int i = 0; i < n / 2; i++)
        {
            circuit.AddInstruction(new QuantumInstruction(
                new SWAPGate(), new List<int> { qubits[i], qubits[n - 1 - i] }));
        }

        return circuit;
    }
}
