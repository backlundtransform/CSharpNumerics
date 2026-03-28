using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Phase oracle gate — an n-qubit diagonal unitary that flips the phase of
/// specified basis states: |w⟩ → −|w⟩ for marked states, |x⟩ → |x⟩ otherwise.
///
/// Used as the oracle in Grover's search algorithm and as the zero-state
/// phase flip inside the Grover diffusion operator.
/// </summary>
public class PhaseOracle : QuantumGate
{
    private readonly int _qubitCount;
    private readonly HashSet<int> _markedStates;

    /// <summary>
    /// Creates a phase oracle for the given qubit count and marked basis states.
    /// </summary>
    /// <param name="qubitCount">Number of qubits (search space = 2^qubitCount).</param>
    /// <param name="markedStates">Basis-state indices whose phase will be flipped.</param>
    public PhaseOracle(int qubitCount, params int[] markedStates)
    {
        if (qubitCount <= 0)
            throw new ArgumentException("Qubit count must be positive.", nameof(qubitCount));

        int stateCount = 1 << qubitCount;
        _markedStates = new HashSet<int>();
        foreach (var s in markedStates)
        {
            if (s < 0 || s >= stateCount)
                throw new ArgumentOutOfRangeException(nameof(markedStates),
                    $"Marked state {s} is out of range for {qubitCount} qubits (0..{stateCount - 1}).");
            _markedStates.Add(s);
        }

        _qubitCount = qubitCount;
    }

    public override int QubitCount => _qubitCount;

    public override ComplexMatrix GetMatrix()
    {
        int size = 1 << _qubitCount;
        var zero = new ComplexNumber(0, 0);
        var one = new ComplexNumber(1, 0);
        var negOne = new ComplexNumber(-1, 0);

        var m = new ComplexNumber[size, size];
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                m[i, j] = zero;

        for (int i = 0; i < size; i++)
            m[i, i] = _markedStates.Contains(i) ? negOne : one;

        return new ComplexMatrix(m);
    }
}
