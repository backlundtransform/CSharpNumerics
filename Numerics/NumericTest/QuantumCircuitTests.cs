using Microsoft.VisualStudio.TestTools.UnitTesting;
using CSharpNumerics.Physics.Quantum;
using CSharpNumerics.Engines.Quantum;
using CSharpNumerics.Numerics.Objects;
using System.Collections.Generic;

namespace NumericTest;

[TestClass]
public class QuantumCircuitTests
{
    private const double Tolerance = 1e-10;

    // ── Gate matrices ──────────────────────────────────────────

    [TestMethod]
    public void HadamardGate_MatrixIs2x2()
    {
        var h = new HadamardGate();
        Assert.AreEqual(1, h.QubitCount);
        var m = h.GetMatrix();
        Assert.AreEqual(2, m.rowLength);
        Assert.AreEqual(2, m.columnLength);
    }

    [TestMethod]
    public void PauliXGate_MatrixIs2x2()
    {
        var x = new PauliXGate();
        Assert.AreEqual(1, x.QubitCount);
        var m = x.GetMatrix();
        Assert.AreEqual(2, m.rowLength);
        Assert.AreEqual(2, m.columnLength);
    }

    [TestMethod]
    public void CNOTGate_MatrixIs4x4()
    {
        var cnot = new CNOTGate();
        Assert.AreEqual(2, cnot.QubitCount);
        var m = cnot.GetMatrix();
        Assert.AreEqual(4, m.rowLength);
        Assert.AreEqual(4, m.columnLength);
    }

    // ── Single-qubit circuits ──────────────────────────────────

    [TestMethod]
    public void Hadamard_CreatesEqualSuperposition()
    {
        var circuit = new QuantumCircuit(1);
        circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));

        var state = new QuantumSimulator().Run(circuit);

        Assert.AreEqual(0.5, state.GetProbability(0), Tolerance);
        Assert.AreEqual(0.5, state.GetProbability(1), Tolerance);
    }

    [TestMethod]
    public void PauliX_FlipsBit()
    {
        var circuit = new QuantumCircuit(1);
        circuit.AddInstruction(new QuantumInstruction(new PauliXGate(), new List<int> { 0 }));

        var state = new QuantumSimulator().Run(circuit);

        Assert.AreEqual(0.0, state.GetProbability(0), Tolerance);
        Assert.AreEqual(1.0, state.GetProbability(1), Tolerance);
    }

    [TestMethod]
    public void DoubleHadamard_ReturnsToOriginal()
    {
        var circuit = new QuantumCircuit(1);
        circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));
        circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));

        var state = new QuantumSimulator().Run(circuit);

        Assert.AreEqual(1.0, state.GetProbability(0), Tolerance);
        Assert.AreEqual(0.0, state.GetProbability(1), Tolerance);
    }

    [TestMethod]
    public void DoublePauliX_ReturnsToOriginal()
    {
        var circuit = new QuantumCircuit(1);
        circuit.AddInstruction(new QuantumInstruction(new PauliXGate(), new List<int> { 0 }));
        circuit.AddInstruction(new QuantumInstruction(new PauliXGate(), new List<int> { 0 }));

        var state = new QuantumSimulator().Run(circuit);

        Assert.AreEqual(1.0, state.GetProbability(0), Tolerance);
        Assert.AreEqual(0.0, state.GetProbability(1), Tolerance);
    }

    // ── Multi-qubit / entanglement ─────────────────────────────

    [TestMethod]
    public void BellState_HadamardThenCNOT()
    {
        // H|0⟩ ⊗ |0⟩  →  (|00⟩ + |10⟩)/√2
        // CNOT(control=0, target=1)  →  (|00⟩ + |11⟩)/√2
        var circuit = new QuantumCircuit(2);
        circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));
        circuit.AddInstruction(new QuantumInstruction(new CNOTGate(), new List<int> { 0, 1 }));

        var state = new QuantumSimulator().Run(circuit);

        Assert.AreEqual(0.5, state.GetProbability(0), Tolerance);  // |00⟩
        Assert.AreEqual(0.0, state.GetProbability(1), Tolerance);  // |01⟩
        Assert.AreEqual(0.0, state.GetProbability(2), Tolerance);  // |10⟩
        Assert.AreEqual(0.5, state.GetProbability(3), Tolerance);  // |11⟩
    }

    [TestMethod]
    public void ThreeQubit_ProbabilitiesSumToOne()
    {
        var circuit = new QuantumCircuit(3);
        circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));
        circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 1 }));
        circuit.AddInstruction(new QuantumInstruction(new PauliXGate(), new List<int> { 2 }));

        var state = new QuantumSimulator().Run(circuit);

        var probs = state.GetProbabilities();
        double sum = 0;
        for (int i = 0; i < probs.Length; i++)
            sum += probs[i];

        Assert.AreEqual(1.0, sum, Tolerance);
    }

    [TestMethod]
    public void ThreeQubit_HadamardOnAll_UniformDistribution()
    {
        var circuit = new QuantumCircuit(3);
        circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));
        circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 1 }));
        circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 2 }));

        var state = new QuantumSimulator().Run(circuit);

        // Each of the 8 basis states should have probability 1/8
        for (int i = 0; i < 8; i++)
            Assert.AreEqual(0.125, state.GetProbability(i), Tolerance);
    }

    // ── QuantumState ───────────────────────────────────────────

    [TestMethod]
    public void QuantumState_GetProbabilities_ReturnsVectorN()
    {
        var circuit = new QuantumCircuit(2);
        circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));

        var state = new QuantumSimulator().Run(circuit);
        var probs = state.GetProbabilities();

        Assert.AreEqual(4, probs.Length);
    }

    [TestMethod]
    public void QuantumState_QubitCount()
    {
        var amps = new ComplexVectorN(new ComplexNumber[]
        {
            new ComplexNumber(1, 0),
            new ComplexNumber(0, 0),
            new ComplexNumber(0, 0),
            new ComplexNumber(0, 0)
        });
        var state = new QuantumState(amps);
        Assert.AreEqual(2, state.QubitCount);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void QuantumState_NonPowerOfTwo_Throws()
    {
        new QuantumState(new ComplexVectorN(new ComplexNumber[]
        {
            new ComplexNumber(1, 0),
            new ComplexNumber(0, 0),
            new ComplexNumber(0, 0)
        }));
    }

    // ── QuantumCircuit ─────────────────────────────────────────

    [TestMethod]
    public void QuantumCircuit_Properties()
    {
        var circuit = new QuantumCircuit(5);
        Assert.AreEqual(5, circuit.QubitCount);
        Assert.AreEqual(0, circuit.Instructions.Count);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentOutOfRangeException))]
    public void QuantumCircuit_OutOfRangeQubit_Throws()
    {
        var circuit = new QuantumCircuit(2);
        circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 3 }));
    }

    // ── QuantumInstruction validation ──────────────────────────

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void QuantumInstruction_WrongQubitCount_Throws()
    {
        // CNOT needs 2 qubits, passing 1
        new QuantumInstruction(new CNOTGate(), new List<int> { 0 });
    }
}
