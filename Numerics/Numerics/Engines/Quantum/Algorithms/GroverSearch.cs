using CSharpNumerics.Physics.Quantum;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Quantum.Algorithms;

/// <summary>
/// Grover's search algorithm.
/// Generates a <see cref="QuantumCircuit"/> that amplifies the amplitude of marked
/// basis states in an unstructured search space. Given M marked states in a space
/// of N = 2^n, the algorithm finds a marked state with high probability after
/// ≈ (π/4)√(N/M) iterations.
///
/// The circuit consists of:
/// 1. Hadamard on all search qubits (uniform superposition)
/// 2. Repeated Grover iterations:
///    a. Phase oracle — flips phase of marked states
///    b. Diffusion operator — reflects about the mean amplitude
/// </summary>
public static class GroverSearch
{
    /// <summary>
    /// Creates a Grover search circuit.
    /// </summary>
    /// <param name="totalQubits">Total number of qubits in the circuit.</param>
    /// <param name="searchQubits">Qubit indices forming the search space.</param>
    /// <param name="markedStates">
    /// Basis-state indices (relative to the search space) to search for.
    /// For example, if searchQubits = {0,1,2} and you want to find |101⟩ → markedStates = {5}.
    /// </param>
    /// <param name="iterations">
    /// Number of Grover iterations. If null, uses the optimal count ⌊π/4 · √(N/M)⌋.
    /// </param>
    public static QuantumCircuit CreateCircuit(int totalQubits, int[] searchQubits,
        int[] markedStates, int? iterations = null)
    {
        if (searchQubits == null || searchQubits.Length == 0)
            throw new ArgumentException("At least one search qubit must be specified.", nameof(searchQubits));
        if (markedStates == null || markedStates.Length == 0)
            throw new ArgumentException("At least one marked state must be specified.", nameof(markedStates));

        int n = searchQubits.Length;
        int iterCount = iterations ?? OptimalIterations(n, markedStates.Length);

        var circuit = new QuantumCircuit(totalQubits);
        var qubits = new List<int>(searchQubits);

        // 1. Initial superposition: H on all search qubits
        foreach (int q in searchQubits)
            circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { q }));

        // Pre-build oracle and zero-state phase flip gates
        var oracle = new PhaseOracle(n, markedStates);
        var zeroFlip = new PhaseOracle(n, 0);

        // 2. Grover iterations
        for (int t = 0; t < iterCount; t++)
        {
            // Oracle: flip phase of marked states
            circuit.AddInstruction(new QuantumInstruction(oracle, qubits));

            // Diffusion: H → phase-flip |0⟩ → H
            foreach (int q in searchQubits)
                circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { q }));

            circuit.AddInstruction(new QuantumInstruction(zeroFlip, qubits));

            foreach (int q in searchQubits)
                circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { q }));
        }

        return circuit;
    }

    /// <summary>
    /// Computes the optimal number of Grover iterations: ⌊π/4 · √(N/M)⌋,
    /// where N = 2^<paramref name="searchSpaceQubits"/> and M = <paramref name="markedCount"/>.
    /// </summary>
    /// <param name="searchSpaceQubits">Number of qubits in the search space.</param>
    /// <param name="markedCount">Number of marked (target) states.</param>
    public static int OptimalIterations(int searchSpaceQubits, int markedCount = 1)
    {
        if (searchSpaceQubits <= 0)
            throw new ArgumentException("Search space must have at least one qubit.", nameof(searchSpaceQubits));
        if (markedCount <= 0)
            throw new ArgumentException("Must have at least one marked state.", nameof(markedCount));

        int N = 1 << searchSpaceQubits;
        return Math.Max(1, (int)Math.Floor(Math.PI / 4.0 * Math.Sqrt((double)N / markedCount)));
    }
}
