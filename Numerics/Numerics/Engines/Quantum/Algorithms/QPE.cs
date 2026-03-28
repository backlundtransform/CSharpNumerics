using CSharpNumerics.Physics.Quantum;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Quantum.Algorithms;

/// <summary>
/// Quantum Phase Estimation (QPE).
/// Given a unitary gate U and an eigenstate |ψ⟩ such that U|ψ⟩ = e^(2πiφ)|ψ⟩,
/// QPE estimates the phase φ to t bits of precision using t counting qubits.
///
/// Circuit structure:
/// 1. H on all counting qubits
/// 2. Controlled-U^(2^k) from counting qubit k to target qubits, for k = 0..t-1
/// 3. Inverse QFT on counting qubits
///
/// After measurement of the counting register, the result ≈ φ · 2^t.
/// </summary>
public static class QPE
{
    /// <summary>
    /// Creates a QPE circuit.
    /// </summary>
    /// <param name="totalQubits">Total number of qubits in the circuit.</param>
    /// <param name="countingQubits">Qubit indices for the counting register (t qubits).</param>
    /// <param name="targetQubits">Qubit indices for the target register (must hold an eigenstate of U).</param>
    /// <param name="unitaryGate">The unitary gate U whose eigenvalue phase is to be estimated.</param>
    public static QuantumCircuit CreateCircuit(int totalQubits, int[] countingQubits,
        int[] targetQubits, QuantumGate unitaryGate)
    {
        if (countingQubits == null || countingQubits.Length == 0)
            throw new ArgumentException("At least one counting qubit is required.", nameof(countingQubits));
        if (targetQubits == null || targetQubits.Length == 0)
            throw new ArgumentException("At least one target qubit is required.", nameof(targetQubits));
        if (unitaryGate == null)
            throw new ArgumentNullException(nameof(unitaryGate));
        if (unitaryGate.QubitCount != targetQubits.Length)
            throw new ArgumentException(
                $"Unitary gate acts on {unitaryGate.QubitCount} qubit(s) but {targetQubits.Length} target qubit(s) provided.",
                nameof(targetQubits));

        var circuit = new QuantumCircuit(totalQubits);
        int t = countingQubits.Length;

        // 1. Hadamard on all counting qubits
        foreach (int q in countingQubits)
            circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { q }));

        // 2. Controlled-U^(2^(t-1-k)): countingQubits[0] = MSB controls most repetitions,
        //    matching the QFT convention where qubits[0] is the most significant bit.
        var controlledU = new ControlledGate(unitaryGate);
        for (int k = 0; k < t; k++)
        {
            int repetitions = 1 << (t - 1 - k);

            // Build qubit index list: [countingQubits[k], targetQubits[0], targetQubits[1], ...]
            var qubits = new List<int> { countingQubits[k] };
            qubits.AddRange(targetQubits);

            for (int rep = 0; rep < repetitions; rep++)
                circuit.AddInstruction(new QuantumInstruction(controlledU, qubits));
        }

        // 3. Inverse QFT on counting qubits
        var iqft = InverseQFT.CreateCircuit(totalQubits, countingQubits);
        foreach (var instruction in iqft.Instructions)
            circuit.AddInstruction(instruction);

        return circuit;
    }
}
