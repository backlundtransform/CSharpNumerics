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

    /// <summary>SWAP gate between qubit <paramref name="qubit1"/> and <paramref name="qubit2"/>.</summary>
    public QuantumCircuitBuilder SWAP(int qubit1, int qubit2)
    {
        _circuit.AddInstruction(new QuantumInstruction(new SWAPGate(), new List<int> { qubit1, qubit2 }));
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
