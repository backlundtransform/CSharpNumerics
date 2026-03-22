using CSharpNumerics.Physics.Quantum;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Quantum;

/// <summary>
/// Pairs a <see cref="QuantumGate"/> with the qubit indices it targets.
/// </summary>
public class QuantumInstruction
{
    /// <summary>Physical qubit indices this instruction acts on.</summary>
    public List<int> QubitIndices { get; }

    /// <summary>The gate to apply.</summary>
    public QuantumGate Gate { get; }

    public QuantumInstruction(QuantumGate gate, List<int> qubitIndices)
    {
        if (gate == null) throw new ArgumentNullException(nameof(gate));
        if (qubitIndices == null) throw new ArgumentNullException(nameof(qubitIndices));
        if (qubitIndices.Count != gate.QubitCount)
            throw new ArgumentException($"Gate requires {gate.QubitCount} qubit(s) but {qubitIndices.Count} index/indices were provided.");

        Gate = gate;
        QubitIndices = qubitIndices;
    }
}
