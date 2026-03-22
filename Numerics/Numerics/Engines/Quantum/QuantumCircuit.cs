using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Quantum;

/// <summary>
/// A quantum circuit consisting of a fixed number of qubits and a sequence of
/// <see cref="QuantumInstruction"/>s to be executed in order.
/// </summary>
public class QuantumCircuit
{
    /// <summary>Number of qubits in the circuit.</summary>
    public int QubitCount { get; }

    /// <summary>Ordered list of instructions (gates + target qubits).</summary>
    public List<QuantumInstruction> Instructions { get; }

    public QuantumCircuit(int qubitCount)
    {
        if (qubitCount <= 0)
            throw new ArgumentException("Qubit count must be positive.", nameof(qubitCount));

        QubitCount = qubitCount;
        Instructions = new List<QuantumInstruction>();
    }

    /// <summary>Appends an instruction to the circuit.</summary>
    public void AddInstruction(QuantumInstruction instruction)
    {
        if (instruction == null) throw new ArgumentNullException(nameof(instruction));

        foreach (var idx in instruction.QubitIndices)
        {
            if (idx < 0 || idx >= QubitCount)
                throw new ArgumentOutOfRangeException(
                    nameof(instruction),
                    $"Qubit index {idx} is out of range for a {QubitCount}-qubit circuit.");
        }

        Instructions.Add(instruction);
    }
}
