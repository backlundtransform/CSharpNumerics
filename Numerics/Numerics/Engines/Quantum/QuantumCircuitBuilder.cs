using CSharpNumerics.Engines.Quantum.Algorithms;
using CSharpNumerics.Physics.Quantum;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Quantum;

/// <summary>
/// Fluent builder for constructing quantum circuits.
/// Usage: <c>QuantumCircuitBuilder.New(2).H(0).CNOT(0, 1).Build()</c>
/// </summary>
public class QuantumCircuitBuilder
{
    private readonly QuantumCircuit _circuit;

    private QuantumCircuitBuilder(int qubitCount)
    {
        _circuit = new QuantumCircuit(qubitCount);
    }

    /// <summary>Creates a new builder for a circuit with the given number of qubits.</summary>
    public static QuantumCircuitBuilder New(int qubitCount) => new QuantumCircuitBuilder(qubitCount);

    /// <summary>Builds and returns the constructed circuit.</summary>
    public QuantumCircuit Build() => _circuit;

    // ── Single-qubit gates ─────────────────────────────────────

    /// <summary>Hadamard gate on qubit <paramref name="qubit"/>.</summary>
    public QuantumCircuitBuilder H(int qubit)
    {
        _circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { qubit }));
        return this;
    }

    /// <summary>Pauli-X (NOT) gate on qubit <paramref name="qubit"/>.</summary>
    public QuantumCircuitBuilder X(int qubit)
    {
        _circuit.AddInstruction(new QuantumInstruction(new PauliXGate(), new List<int> { qubit }));
        return this;
    }

    /// <summary>Pauli-Z gate on qubit <paramref name="qubit"/>.</summary>
    public QuantumCircuitBuilder Z(int qubit)
    {
        _circuit.AddInstruction(new QuantumInstruction(new PauliZGate(), new List<int> { qubit }));
        return this;
    }

    /// <summary>S (√Z) gate on qubit <paramref name="qubit"/>.</summary>
    public QuantumCircuitBuilder S(int qubit)
    {
        _circuit.AddInstruction(new QuantumInstruction(new SGate(), new List<int> { qubit }));
        return this;
    }

    /// <summary>T (π/8) gate on qubit <paramref name="qubit"/>.</summary>
    public QuantumCircuitBuilder T(int qubit)
    {
        _circuit.AddInstruction(new QuantumInstruction(new TGate(), new List<int> { qubit }));
        return this;
    }

    /// <summary>Rotation about X-axis by <paramref name="theta"/> radians.</summary>
    public QuantumCircuitBuilder Rx(int qubit, double theta)
    {
        _circuit.AddInstruction(new QuantumInstruction(new RxGate(theta), new List<int> { qubit }));
        return this;
    }

    /// <summary>Rotation about Y-axis by <paramref name="theta"/> radians.</summary>
    public QuantumCircuitBuilder Ry(int qubit, double theta)
    {
        _circuit.AddInstruction(new QuantumInstruction(new RyGate(theta), new List<int> { qubit }));
        return this;
    }

    /// <summary>Rotation about Z-axis by <paramref name="theta"/> radians.</summary>
    public QuantumCircuitBuilder Rz(int qubit, double theta)
    {
        _circuit.AddInstruction(new QuantumInstruction(new RzGate(theta), new List<int> { qubit }));
        return this;
    }

    // ── Two-qubit gates ────────────────────────────────────────

    /// <summary>CNOT (controlled-X) with <paramref name="control"/> and <paramref name="target"/>.</summary>
    public QuantumCircuitBuilder CNOT(int control, int target)
    {
        _circuit.AddInstruction(new QuantumInstruction(new CNOTGate(), new List<int> { control, target }));
        return this;
    }

    /// <summary>CZ (controlled-Z) with <paramref name="control"/> and <paramref name="target"/>.</summary>
    public QuantumCircuitBuilder CZ(int control, int target)
    {
        _circuit.AddInstruction(new QuantumInstruction(new CZGate(), new List<int> { control, target }));
        return this;
    }

    /// <summary>Controlled-Phase gate CP(θ) with <paramref name="control"/> and <paramref name="target"/>.</summary>
    public QuantumCircuitBuilder CPhase(int control, int target, double theta)
    {
        _circuit.AddInstruction(new QuantumInstruction(new CPhaseGate(theta), new List<int> { control, target }));
        return this;
    }

    /// <summary>SWAP gate between qubit <paramref name="qubit1"/> and <paramref name="qubit2"/>.</summary>
    public QuantumCircuitBuilder SWAP(int qubit1, int qubit2)
    {
        _circuit.AddInstruction(new QuantumInstruction(new SWAPGate(), new List<int> { qubit1, qubit2 }));
        return this;
    }

    // ── Three-qubit gates ──────────────────────────────────────

    /// <summary>Toffoli (CCNOT) gate with two controls and one target.</summary>
    public QuantumCircuitBuilder Toffoli(int control1, int control2, int target)
    {
        _circuit.AddInstruction(new QuantumInstruction(new ToffoliGate(), new List<int> { control1, control2, target }));
        return this;
    }

    /// <summary>Fredkin (CSWAP) gate: swaps target1 and target2 when control is |1⟩.</summary>
    public QuantumCircuitBuilder Fredkin(int control, int target1, int target2)
    {
        _circuit.AddInstruction(new QuantumInstruction(new FredkinGate(), new List<int> { control, target1, target2 }));
        return this;
    }

    // ── Additional single-qubit gates ──────────────────────────

    /// <summary>Pauli-Y gate on qubit <paramref name="qubit"/>.</summary>
    public QuantumCircuitBuilder Y(int qubit)
    {
        _circuit.AddInstruction(new QuantumInstruction(new PauliYGate(), new List<int> { qubit }));
        return this;
    }

    /// <summary>General phase gate P(θ) on qubit <paramref name="qubit"/>.</summary>
    public QuantumCircuitBuilder Phase(int qubit, double theta)
    {
        _circuit.AddInstruction(new QuantumInstruction(new PhaseGate(theta), new List<int> { qubit }));
        return this;
    }

    // ── Composite algorithms ───────────────────────────────────

    /// <summary>Applies the Quantum Fourier Transform to the specified qubits.</summary>
    public QuantumCircuitBuilder ApplyQFT(params int[] qubits)
    {
        var qftCircuit = QFT.CreateCircuit(_circuit.QubitCount, qubits);
        foreach (var instruction in qftCircuit.Instructions)
            _circuit.AddInstruction(instruction);
        return this;
    }

    /// <summary>Applies the inverse Quantum Fourier Transform to the specified qubits.</summary>
    public QuantumCircuitBuilder ApplyInverseQFT(params int[] qubits)
    {
        var iqftCircuit = InverseQFT.CreateCircuit(_circuit.QubitCount, qubits);
        foreach (var instruction in iqftCircuit.Instructions)
            _circuit.AddInstruction(instruction);
        return this;
    }

    /// <summary>
    /// Applies Grover's search algorithm to the specified qubits, searching for the given marked states.
    /// </summary>
    /// <param name="markedStates">Basis-state indices to search for (relative to the search qubit space).</param>
    /// <param name="searchQubits">Qubit indices forming the search space.</param>
    public QuantumCircuitBuilder ApplyGrover(int[] markedStates, params int[] searchQubits)
    {
        var groverCircuit = GroverSearch.CreateCircuit(_circuit.QubitCount, searchQubits, markedStates);
        foreach (var instruction in groverCircuit.Instructions)
            _circuit.AddInstruction(instruction);
        return this;
    }

    /// <summary>
    /// Applies Grover's search algorithm with a specified number of iterations.
    /// </summary>
    public QuantumCircuitBuilder ApplyGrover(int[] markedStates, int iterations, params int[] searchQubits)
    {
        var groverCircuit = GroverSearch.CreateCircuit(_circuit.QubitCount, searchQubits, markedStates, iterations);
        foreach (var instruction in groverCircuit.Instructions)
            _circuit.AddInstruction(instruction);
        return this;
    }

    /// <summary>
    /// Applies Quantum Phase Estimation: H on counting qubits, controlled-U^(2^k), inverse QFT.
    /// </summary>
    /// <param name="countingQubits">Qubit indices for the counting register.</param>
    /// <param name="targetQubits">Qubit indices for the target register (must hold eigenstate).</param>
    /// <param name="unitaryGate">The unitary gate whose eigenvalue phase is estimated.</param>
    public QuantumCircuitBuilder ApplyQPE(int[] countingQubits, int[] targetQubits, QuantumGate unitaryGate)
    {
        var qpeCircuit = QPE.CreateCircuit(_circuit.QubitCount, countingQubits, targetQubits, unitaryGate);
        foreach (var instruction in qpeCircuit.Instructions)
            _circuit.AddInstruction(instruction);
        return this;
    }

    /// <summary>
    /// Applies a controlled version of <paramref name="innerGate"/>.
    /// First qubit index is the control, remaining are the inner gate's targets.
    /// </summary>
    public QuantumCircuitBuilder Controlled(QuantumGate innerGate, int control, params int[] targets)
    {
        var controlled = new ControlledGate(innerGate);
        var qubits = new List<int> { control };
        qubits.AddRange(targets);
        _circuit.AddInstruction(new QuantumInstruction(controlled, qubits));
        return this;
    }

    // ── Custom gate ────────────────────────────────────────────

    /// <summary>Applies an arbitrary gate to the specified qubits.</summary>
    public QuantumCircuitBuilder Gate(QuantumGate gate, params int[] qubits)
    {
        _circuit.AddInstruction(new QuantumInstruction(gate, new List<int>(qubits)));
        return this;
    }
}
