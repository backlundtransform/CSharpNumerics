using Microsoft.VisualStudio.TestTools.UnitTesting;
using CSharpNumerics.Physics.Quantum;
using CSharpNumerics.Engines.Quantum;
using CSharpNumerics.Numerics.Objects;
using System;
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

    // ── PauliZ gate ────────────────────────────────────────────

    [TestMethod]
    public void PauliZ_LeavesZeroUnchanged()
    {
        var circuit = new QuantumCircuit(1);
        circuit.AddInstruction(new QuantumInstruction(new PauliZGate(), new List<int> { 0 }));

        var state = new QuantumSimulator().Run(circuit);

        Assert.AreEqual(1.0, state.GetProbability(0), Tolerance);
        Assert.AreEqual(0.0, state.GetProbability(1), Tolerance);
    }

    [TestMethod]
    public void PauliZ_FlipsPhaseOfOne()
    {
        // X|0⟩ = |1⟩, then Z|1⟩ = −|1⟩ (same probability, different phase)
        var circuit = new QuantumCircuit(1);
        circuit.AddInstruction(new QuantumInstruction(new PauliXGate(), new List<int> { 0 }));
        circuit.AddInstruction(new QuantumInstruction(new PauliZGate(), new List<int> { 0 }));

        var state = new QuantumSimulator().Run(circuit);

        Assert.AreEqual(0.0, state.GetProbability(0), Tolerance);
        Assert.AreEqual(1.0, state.GetProbability(1), Tolerance);
        // Amplitude should be -1
        Assert.AreEqual(-1.0, state.Amplitudes[1].realPart, Tolerance);
    }

    [TestMethod]
    public void DoublePauliZ_ReturnsToOriginal()
    {
        // Z² = I
        var circuit = new QuantumCircuit(1);
        circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));
        circuit.AddInstruction(new QuantumInstruction(new PauliZGate(), new List<int> { 0 }));
        circuit.AddInstruction(new QuantumInstruction(new PauliZGate(), new List<int> { 0 }));

        var state = new QuantumSimulator().Run(circuit);

        // Should be back to H|0⟩ = (|0⟩+|1⟩)/√2
        Assert.AreEqual(0.5, state.GetProbability(0), Tolerance);
        Assert.AreEqual(0.5, state.GetProbability(1), Tolerance);
    }

    // ── S gate ─────────────────────────────────────────────────

    [TestMethod]
    public void SGate_SquaredEqualsZ()
    {
        // S²|ψ⟩ should equal Z|ψ⟩ for any state
        // Test on H|0⟩ = (|0⟩+|1⟩)/√2
        var circuitS2 = new QuantumCircuit(1);
        circuitS2.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));
        circuitS2.AddInstruction(new QuantumInstruction(new SGate(), new List<int> { 0 }));
        circuitS2.AddInstruction(new QuantumInstruction(new SGate(), new List<int> { 0 }));

        var circuitZ = new QuantumCircuit(1);
        circuitZ.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));
        circuitZ.AddInstruction(new QuantumInstruction(new PauliZGate(), new List<int> { 0 }));

        var sim = new QuantumSimulator();
        var stateS2 = sim.Run(circuitS2);
        var stateZ = sim.Run(circuitZ);

        for (int i = 0; i < 2; i++)
        {
            Assert.AreEqual(stateZ.Amplitudes[i].realPart, stateS2.Amplitudes[i].realPart, Tolerance);
            Assert.AreEqual(stateZ.Amplitudes[i].imaginaryPart, stateS2.Amplitudes[i].imaginaryPart, Tolerance);
        }
    }

    // ── T gate ─────────────────────────────────────────────────

    [TestMethod]
    public void TGate_SquaredEqualsS()
    {
        // T²|ψ⟩ = S|ψ⟩
        var circuitT2 = new QuantumCircuit(1);
        circuitT2.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));
        circuitT2.AddInstruction(new QuantumInstruction(new TGate(), new List<int> { 0 }));
        circuitT2.AddInstruction(new QuantumInstruction(new TGate(), new List<int> { 0 }));

        var circuitS = new QuantumCircuit(1);
        circuitS.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));
        circuitS.AddInstruction(new QuantumInstruction(new SGate(), new List<int> { 0 }));

        var sim = new QuantumSimulator();
        var stateT2 = sim.Run(circuitT2);
        var stateS = sim.Run(circuitS);

        for (int i = 0; i < 2; i++)
        {
            Assert.AreEqual(stateS.Amplitudes[i].realPart, stateT2.Amplitudes[i].realPart, Tolerance);
            Assert.AreEqual(stateS.Amplitudes[i].imaginaryPart, stateT2.Amplitudes[i].imaginaryPart, Tolerance);
        }
    }

    // ── Rotation gates ─────────────────────────────────────────

    [TestMethod]
    public void RxPi_EquivalentToPauliX_UpToGlobalPhase()
    {
        // Rx(π) = -iX → same probabilities as X
        var circuit = new QuantumCircuit(1);
        circuit.AddInstruction(new QuantumInstruction(new RxGate(Math.PI), new List<int> { 0 }));

        var state = new QuantumSimulator().Run(circuit);

        Assert.AreEqual(0.0, state.GetProbability(0), Tolerance);
        Assert.AreEqual(1.0, state.GetProbability(1), Tolerance);
    }

    [TestMethod]
    public void Rx2Pi_ReturnsToOriginal_UpToGlobalPhase()
    {
        // Rx(2π) = -I → probabilities unchanged
        var circuit = new QuantumCircuit(1);
        circuit.AddInstruction(new QuantumInstruction(new RxGate(2 * Math.PI), new List<int> { 0 }));

        var state = new QuantumSimulator().Run(circuit);

        Assert.AreEqual(1.0, state.GetProbability(0), Tolerance);
        Assert.AreEqual(0.0, state.GetProbability(1), Tolerance);
    }

    [TestMethod]
    public void RyPi_EquivalentToPauliY_UpToGlobalPhase()
    {
        // Ry(π)|0⟩ = |1⟩
        var circuit = new QuantumCircuit(1);
        circuit.AddInstruction(new QuantumInstruction(new RyGate(Math.PI), new List<int> { 0 }));

        var state = new QuantumSimulator().Run(circuit);

        Assert.AreEqual(0.0, state.GetProbability(0), Tolerance);
        Assert.AreEqual(1.0, state.GetProbability(1), Tolerance);
    }

    [TestMethod]
    public void RyHalfPi_CreatesEqualSuperposition()
    {
        // Ry(π/2)|0⟩ = (|0⟩+|1⟩)/√2
        var circuit = new QuantumCircuit(1);
        circuit.AddInstruction(new QuantumInstruction(new RyGate(Math.PI / 2), new List<int> { 0 }));

        var state = new QuantumSimulator().Run(circuit);

        Assert.AreEqual(0.5, state.GetProbability(0), Tolerance);
        Assert.AreEqual(0.5, state.GetProbability(1), Tolerance);
    }

    [TestMethod]
    public void RzPi_EquivalentToPauliZ_UpToGlobalPhase()
    {
        // Rz(π) = -iZ → same probabilities as Z on H|0⟩
        var circuit = new QuantumCircuit(1);
        circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));
        circuit.AddInstruction(new QuantumInstruction(new RzGate(Math.PI), new List<int> { 0 }));

        var state = new QuantumSimulator().Run(circuit);

        Assert.AreEqual(0.5, state.GetProbability(0), Tolerance);
        Assert.AreEqual(0.5, state.GetProbability(1), Tolerance);
    }

    [TestMethod]
    public void Rz_PreservesProbabilities_OnBasisState()
    {
        // Rz on |0⟩ only changes global phase
        var circuit = new QuantumCircuit(1);
        circuit.AddInstruction(new QuantumInstruction(new RzGate(1.23), new List<int> { 0 }));

        var state = new QuantumSimulator().Run(circuit);

        Assert.AreEqual(1.0, state.GetProbability(0), Tolerance);
        Assert.AreEqual(0.0, state.GetProbability(1), Tolerance);
    }

    // ── CZ gate ────────────────────────────────────────────────

    [TestMethod]
    public void CZGate_FlipsPhaseOf11()
    {
        // Prepare |11⟩ via X⊗X, then apply CZ
        var circuit = new QuantumCircuit(2);
        circuit.AddInstruction(new QuantumInstruction(new PauliXGate(), new List<int> { 0 }));
        circuit.AddInstruction(new QuantumInstruction(new PauliXGate(), new List<int> { 1 }));
        circuit.AddInstruction(new QuantumInstruction(new CZGate(), new List<int> { 0, 1 }));

        var state = new QuantumSimulator().Run(circuit);

        // Still |11⟩ with probability 1, but amplitude = -1
        Assert.AreEqual(0.0, state.GetProbability(0), Tolerance);
        Assert.AreEqual(1.0, state.GetProbability(3), Tolerance);
        Assert.AreEqual(-1.0, state.Amplitudes[3].realPart, Tolerance);
    }

    [TestMethod]
    public void CZGate_LeavesOtherStatesUnchanged()
    {
        // |01⟩ → |01⟩ (no phase flip)
        var circuit = new QuantumCircuit(2);
        circuit.AddInstruction(new QuantumInstruction(new PauliXGate(), new List<int> { 1 }));
        circuit.AddInstruction(new QuantumInstruction(new CZGate(), new List<int> { 0, 1 }));

        var state = new QuantumSimulator().Run(circuit);

        Assert.AreEqual(1.0, state.GetProbability(2), Tolerance); // |01⟩
        Assert.AreEqual(1.0, state.Amplitudes[2].realPart, Tolerance);
    }

    // ── SWAP gate ──────────────────────────────────────────────

    [TestMethod]
    public void SWAPGate_Swaps01To10()
    {
        // Prepare |10⟩ (qubit 0 = |1⟩, qubit 1 = |0⟩), SWAP → |01⟩
        var circuit = new QuantumCircuit(2);
        circuit.AddInstruction(new QuantumInstruction(new PauliXGate(), new List<int> { 0 }));
        circuit.AddInstruction(new QuantumInstruction(new SWAPGate(), new List<int> { 0, 1 }));

        var state = new QuantumSimulator().Run(circuit);

        // |10⟩ (idx 1) → |01⟩ (idx 2)
        Assert.AreEqual(0.0, state.GetProbability(0), Tolerance);
        Assert.AreEqual(0.0, state.GetProbability(1), Tolerance);
        Assert.AreEqual(1.0, state.GetProbability(2), Tolerance);
        Assert.AreEqual(0.0, state.GetProbability(3), Tolerance);
    }

    [TestMethod]
    public void SWAPGate_Leaves00Unchanged()
    {
        var circuit = new QuantumCircuit(2);
        circuit.AddInstruction(new QuantumInstruction(new SWAPGate(), new List<int> { 0, 1 }));

        var state = new QuantumSimulator().Run(circuit);

        Assert.AreEqual(1.0, state.GetProbability(0), Tolerance);
    }

    [TestMethod]
    public void SWAPGate_Leaves11Unchanged()
    {
        var circuit = new QuantumCircuit(2);
        circuit.AddInstruction(new QuantumInstruction(new PauliXGate(), new List<int> { 0 }));
        circuit.AddInstruction(new QuantumInstruction(new PauliXGate(), new List<int> { 1 }));
        circuit.AddInstruction(new QuantumInstruction(new SWAPGate(), new List<int> { 0, 1 }));

        var state = new QuantumSimulator().Run(circuit);

        Assert.AreEqual(1.0, state.GetProbability(3), Tolerance);
    }

    [TestMethod]
    public void DoubleSWAP_ReturnsToOriginal()
    {
        // SWAP² = I
        var circuit = new QuantumCircuit(2);
        circuit.AddInstruction(new QuantumInstruction(new PauliXGate(), new List<int> { 0 }));
        circuit.AddInstruction(new QuantumInstruction(new SWAPGate(), new List<int> { 0, 1 }));
        circuit.AddInstruction(new QuantumInstruction(new SWAPGate(), new List<int> { 0, 1 }));

        var state = new QuantumSimulator().Run(circuit);

        // Back to |10⟩
        Assert.AreEqual(0.0, state.GetProbability(0), Tolerance);
        Assert.AreEqual(1.0, state.GetProbability(1), Tolerance);
        Assert.AreEqual(0.0, state.GetProbability(2), Tolerance);
        Assert.AreEqual(0.0, state.GetProbability(3), Tolerance);
    }

    // ── Bloch vector ───────────────────────────────────────────

    [TestMethod]
    public void BlochVector_ZeroState_NorthPole()
    {
        // |0⟩ → Bloch (0, 0, 1)
        var circuit = new QuantumCircuit(1);
        var state = new QuantumSimulator().Run(circuit);
        var bloch = state.GetBlochVector();

        Assert.AreEqual(0.0, bloch.X, Tolerance);
        Assert.AreEqual(0.0, bloch.Y, Tolerance);
        Assert.AreEqual(1.0, bloch.Z, Tolerance);
        Assert.AreEqual(1.0, bloch.Radius, Tolerance);
    }

    [TestMethod]
    public void BlochVector_OneState_SouthPole()
    {
        // X|0⟩ = |1⟩ → Bloch (0, 0, -1)
        var circuit = new QuantumCircuit(1);
        circuit.AddInstruction(new QuantumInstruction(new PauliXGate(), new List<int> { 0 }));
        var state = new QuantumSimulator().Run(circuit);
        var bloch = state.GetBlochVector();

        Assert.AreEqual(0.0, bloch.X, Tolerance);
        Assert.AreEqual(0.0, bloch.Y, Tolerance);
        Assert.AreEqual(-1.0, bloch.Z, Tolerance);
    }

    [TestMethod]
    public void BlochVector_HadamardState_PositiveX()
    {
        // H|0⟩ = (|0⟩+|1⟩)/√2 → Bloch (1, 0, 0)
        var circuit = new QuantumCircuit(1);
        circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));
        var state = new QuantumSimulator().Run(circuit);
        var bloch = state.GetBlochVector();

        Assert.AreEqual(1.0, bloch.X, Tolerance);
        Assert.AreEqual(0.0, bloch.Y, Tolerance);
        Assert.AreEqual(0.0, bloch.Z, Tolerance);
    }

    [TestMethod]
    public void BlochVector_MinusState_NegativeX()
    {
        // H|0⟩ then Z → (|0⟩−|1⟩)/√2 → Bloch (−1, 0, 0)
        var circuit = new QuantumCircuit(1);
        circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));
        circuit.AddInstruction(new QuantumInstruction(new PauliZGate(), new List<int> { 0 }));
        var state = new QuantumSimulator().Run(circuit);
        var bloch = state.GetBlochVector();

        Assert.AreEqual(-1.0, bloch.X, Tolerance);
        Assert.AreEqual(0.0, bloch.Y, Tolerance);
        Assert.AreEqual(0.0, bloch.Z, Tolerance);
    }

    [TestMethod]
    public void BlochVector_PlusI_PositiveY()
    {
        // H|0⟩ then S → (|0⟩+i|1⟩)/√2 → Bloch (0, 1, 0)
        var circuit = new QuantumCircuit(1);
        circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));
        circuit.AddInstruction(new QuantumInstruction(new SGate(), new List<int> { 0 }));
        var state = new QuantumSimulator().Run(circuit);
        var bloch = state.GetBlochVector();

        Assert.AreEqual(0.0, bloch.X, Tolerance);
        Assert.AreEqual(1.0, bloch.Y, Tolerance);
        Assert.AreEqual(0.0, bloch.Z, Tolerance);
    }

    [TestMethod]
    public void BlochVector_RadiusIsOne_ForPureState()
    {
        // Any single-qubit pure state should have radius = 1
        var circuit = new QuantumCircuit(1);
        circuit.AddInstruction(new QuantumInstruction(new RyGate(1.23), new List<int> { 0 }));
        circuit.AddInstruction(new QuantumInstruction(new RzGate(0.77), new List<int> { 0 }));
        var state = new QuantumSimulator().Run(circuit);
        var bloch = state.GetBlochVector();

        Assert.AreEqual(1.0, bloch.Radius, Tolerance);
    }

    [TestMethod]
    public void BlochVector_ThetaPhi_Consistency()
    {
        // Verify spherical ↔ Cartesian roundtrip
        var circuit = new QuantumCircuit(1);
        circuit.AddInstruction(new QuantumInstruction(new RyGate(Math.PI / 3), new List<int> { 0 }));
        circuit.AddInstruction(new QuantumInstruction(new RzGate(Math.PI / 4), new List<int> { 0 }));
        var state = new QuantumSimulator().Run(circuit);
        var bloch = state.GetBlochVector();

        double theta = bloch.Theta;
        double phi = bloch.Phi;
        double r = bloch.Radius;

        Assert.AreEqual(bloch.X, r * Math.Sin(theta) * Math.Cos(phi), Tolerance);
        Assert.AreEqual(bloch.Y, r * Math.Sin(theta) * Math.Sin(phi), Tolerance);
        Assert.AreEqual(bloch.Z, r * Math.Cos(theta), Tolerance);
    }

    [TestMethod]
    public void BlochVector_ToVector_ReturnsCorrectVector()
    {
        var bloch = BlochVector.FromAmplitudes(
            new ComplexNumber(1, 0), new ComplexNumber(0, 0));
        var v = bloch.ToVector();

        Assert.AreEqual(bloch.X, v.x, Tolerance);
        Assert.AreEqual(bloch.Y, v.y, Tolerance);
        Assert.AreEqual(bloch.Z, v.z, Tolerance);
    }

    [TestMethod]
    [ExpectedException(typeof(InvalidOperationException))]
    public void BlochVector_MultiQubit_Throws()
    {
        var circuit = new QuantumCircuit(2);
        var state = new QuantumSimulator().Run(circuit);
        state.GetBlochVector();
    }

    // ── Fluent builder ─────────────────────────────────────────

    [TestMethod]
    public void Builder_BellState()
    {
        var circuit = QuantumCircuitBuilder.New(2)
            .H(0)
            .CNOT(0, 1)
            .Build();

        var state = new QuantumSimulator().Run(circuit);

        Assert.AreEqual(0.5, state.GetProbability(0), Tolerance);  // |00⟩
        Assert.AreEqual(0.0, state.GetProbability(1), Tolerance);  // |01⟩
        Assert.AreEqual(0.0, state.GetProbability(2), Tolerance);  // |10⟩
        Assert.AreEqual(0.5, state.GetProbability(3), Tolerance);  // |11⟩
    }

    [TestMethod]
    public void Builder_AllSingleQubitGates()
    {
        // Ensure each builder method creates a valid circuit
        var circuit = QuantumCircuitBuilder.New(1)
            .H(0).X(0).Z(0).S(0).T(0)
            .Rx(0, 0.5).Ry(0, 0.5).Rz(0, 0.5)
            .Build();

        Assert.AreEqual(8, circuit.Instructions.Count);
        Assert.AreEqual(1, circuit.QubitCount);
    }

    [TestMethod]
    public void Builder_TwoQubitGates()
    {
        var circuit = QuantumCircuitBuilder.New(2)
            .CNOT(0, 1)
            .CZ(0, 1)
            .SWAP(0, 1)
            .Build();

        Assert.AreEqual(3, circuit.Instructions.Count);
    }

    [TestMethod]
    public void Builder_CustomGate()
    {
        var circuit = QuantumCircuitBuilder.New(2)
            .Gate(new HadamardGate(), 1)
            .Gate(new CNOTGate(), 0, 1)
            .Build();

        var state = new QuantumSimulator().Run(circuit);
        var probs = state.GetProbabilities();
        double sum = 0;
        for (int i = 0; i < probs.Length; i++) sum += probs[i];
        Assert.AreEqual(1.0, sum, Tolerance);
    }

    [TestMethod]
    public void Builder_MatchesManualCircuit()
    {
        // Fluent and manual should produce identical states
        var manual = new QuantumCircuit(2);
        manual.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));
        manual.AddInstruction(new QuantumInstruction(new RyGate(1.23), new List<int> { 1 }));
        manual.AddInstruction(new QuantumInstruction(new CNOTGate(), new List<int> { 0, 1 }));

        var fluent = QuantumCircuitBuilder.New(2)
            .H(0).Ry(1, 1.23).CNOT(0, 1)
            .Build();

        var sim = new QuantumSimulator();
        var s1 = sim.Run(manual);
        var s2 = sim.Run(fluent);

        for (int i = 0; i < s1.Amplitudes.Length; i++)
        {
            Assert.AreEqual(s1.Amplitudes[i].realPart, s2.Amplitudes[i].realPart, Tolerance);
            Assert.AreEqual(s1.Amplitudes[i].imaginaryPart, s2.Amplitudes[i].imaginaryPart, Tolerance);
        }
    }

    // ── Measurement ────────────────────────────────────────────

    [TestMethod]
    public void Measure_BasisState_AlwaysReturnsSame()
    {
        // |0⟩ should always measure 0
        var circuit = QuantumCircuitBuilder.New(1).Build();
        var sim = new QuantumSimulator();
        var rng = new Random(42);

        for (int i = 0; i < 20; i++)
        {
            var state = sim.Run(circuit);
            int result = state.Measure(rng);
            Assert.AreEqual(0, result);
        }
    }

    [TestMethod]
    public void Measure_CollapsesToBasisState()
    {
        // After measurement, state should be a basis state (probability = 1)
        var circuit = QuantumCircuitBuilder.New(1).H(0).Build();
        var state = new QuantumSimulator().Run(circuit);
        var rng = new Random(42);

        int result = state.Measure(rng);

        Assert.AreEqual(1.0, state.GetProbability(result), Tolerance);
        Assert.AreEqual(0.0, state.GetProbability(1 - result), Tolerance);
    }

    [TestMethod]
    public void Measure_Superposition_ReturnsValidState()
    {
        // H|0⟩ → should return 0 or 1
        var circuit = QuantumCircuitBuilder.New(1).H(0).Build();
        var state = new QuantumSimulator().Run(circuit);
        var rng = new Random(42);

        int result = state.Measure(rng);
        Assert.IsTrue(result == 0 || result == 1);
    }

    [TestMethod]
    public void MeasureQubit_CollapsesTargetQubit()
    {
        // Bell state: (|00⟩+|11⟩)/√2 — measuring qubit 0 should collapse both
        var circuit = QuantumCircuitBuilder.New(2).H(0).CNOT(0, 1).Build();
        var state = new QuantumSimulator().Run(circuit);
        var rng = new Random(42);

        int outcome = state.MeasureQubit(0, rng);

        if (outcome == 0)
        {
            Assert.AreEqual(1.0, state.GetProbability(0), Tolerance); // |00⟩
            Assert.AreEqual(0.0, state.GetProbability(3), Tolerance); // |11⟩
        }
        else
        {
            Assert.AreEqual(0.0, state.GetProbability(0), Tolerance); // |00⟩
            Assert.AreEqual(1.0, state.GetProbability(3), Tolerance); // |11⟩
        }
    }

    [TestMethod]
    public void MeasureQubit_PreservesOtherQubit()
    {
        // |+⟩⊗|0⟩ — measuring qubit 1 (always |0⟩) should leave qubit 0 in superposition
        var circuit = QuantumCircuitBuilder.New(2).H(0).Build();
        var state = new QuantumSimulator().Run(circuit);
        var rng = new Random(42);

        int outcome = state.MeasureQubit(1, rng);
        Assert.AreEqual(0, outcome); // qubit 1 is always |0⟩

        // Qubit 0 should still be in superposition: |00⟩ and |01⟩ each 0.5
        // Since qubit 1 → |0⟩, only |00⟩ (idx 0) and |10⟩ (idx 1) survive
        Assert.AreEqual(0.5, state.GetProbability(0), Tolerance);
        Assert.AreEqual(0.5, state.GetProbability(1), Tolerance);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentOutOfRangeException))]
    public void MeasureQubit_InvalidIndex_Throws()
    {
        var state = new QuantumSimulator().Run(QuantumCircuitBuilder.New(1).Build());
        state.MeasureQubit(5, new Random());
    }

    [TestMethod]
    public void Sample_ReturnsCorrectCounts()
    {
        // |0⟩ should always return {0: shots}
        var state = new QuantumSimulator().Run(QuantumCircuitBuilder.New(1).Build());
        var rng = new Random(42);

        var counts = state.Sample(1000, rng);

        Assert.AreEqual(1, counts.Count);
        Assert.AreEqual(1000, counts[0]);
    }

    [TestMethod]
    public void Sample_DoesNotCollapseState()
    {
        // After sampling, probabilities should still be 0.5/0.5
        var state = new QuantumSimulator().Run(QuantumCircuitBuilder.New(1).H(0).Build());
        var rng = new Random(42);

        state.Sample(1000, rng);

        Assert.AreEqual(0.5, state.GetProbability(0), Tolerance);
        Assert.AreEqual(0.5, state.GetProbability(1), Tolerance);
    }

    [TestMethod]
    public void Sample_Superposition_ApproximatelyEqualCounts()
    {
        // H|0⟩ → ~50/50 split with enough shots
        var state = new QuantumSimulator().Run(QuantumCircuitBuilder.New(1).H(0).Build());
        var rng = new Random(42);

        var counts = state.Sample(10000, rng);

        Assert.IsTrue(counts.ContainsKey(0));
        Assert.IsTrue(counts.ContainsKey(1));

        // Each should be roughly 5000 ± 300
        Assert.IsTrue(Math.Abs(counts[0] - 5000) < 500, $"|0⟩ count was {counts[0]}");
        Assert.IsTrue(Math.Abs(counts[1] - 5000) < 500, $"|1⟩ count was {counts[1]}");
    }

    [TestMethod]
    public void Sample_TotalCountsEqualShots()
    {
        var state = new QuantumSimulator().Run(
            QuantumCircuitBuilder.New(2).H(0).H(1).Build());
        var rng = new Random(42);

        var counts = state.Sample(5000, rng);

        int total = 0;
        foreach (var kv in counts) total += kv.Value;
        Assert.AreEqual(5000, total);
    }

    // ── Fidelity ───────────────────────────────────────────────

    [TestMethod]
    public void Fidelity_IdenticalStates_ReturnsOne()
    {
        var state = new QuantumSimulator().Run(QuantumCircuitBuilder.New(1).H(0).Build());
        double f = QuantumFidelity.Fidelity(state, state);
        Assert.AreEqual(1.0, f, Tolerance);
    }

    [TestMethod]
    public void Fidelity_OrthogonalStates_ReturnsZero()
    {
        var s0 = new QuantumSimulator().Run(QuantumCircuitBuilder.New(1).Build());       // |0⟩
        var s1 = new QuantumSimulator().Run(QuantumCircuitBuilder.New(1).X(0).Build());  // |1⟩
        double f = QuantumFidelity.Fidelity(s0, s1);
        Assert.AreEqual(0.0, f, Tolerance);
    }

    [TestMethod]
    public void Fidelity_SuperpositionVsBasis_ReturnsHalf()
    {
        // |⟨0|+⟩|² = |1/√2|² = 0.5
        var s0 = new QuantumSimulator().Run(QuantumCircuitBuilder.New(1).Build());
        var sPlus = new QuantumSimulator().Run(QuantumCircuitBuilder.New(1).H(0).Build());
        double f = QuantumFidelity.Fidelity(s0, sPlus);
        Assert.AreEqual(0.5, f, Tolerance);
    }

    [TestMethod]
    public void Fidelity_IsSymmetric()
    {
        var s1 = new QuantumSimulator().Run(QuantumCircuitBuilder.New(1).Ry(0, 1.0).Build());
        var s2 = new QuantumSimulator().Run(QuantumCircuitBuilder.New(1).Ry(0, 2.0).Build());
        double f12 = QuantumFidelity.Fidelity(s1, s2);
        double f21 = QuantumFidelity.Fidelity(s2, s1);
        Assert.AreEqual(f12, f21, Tolerance);
    }

    [TestMethod]
    public void Fidelity_VectorOverload()
    {
        var v1 = new ComplexVectorN(new ComplexNumber[]
        {
            new ComplexNumber(1, 0),
            new ComplexNumber(0, 0)
        });
        var v2 = new ComplexVectorN(new ComplexNumber[]
        {
            new ComplexNumber(1.0 / Math.Sqrt(2), 0),
            new ComplexNumber(1.0 / Math.Sqrt(2), 0)
        });
        double f = QuantumFidelity.Fidelity(v1, v2);
        Assert.AreEqual(0.5, f, Tolerance);
    }

    [TestMethod]
    public void BlochFidelity_IdenticalVectors_ReturnsOne()
    {
        var b = BlochVector.FromAmplitudes(
            new ComplexNumber(1, 0), new ComplexNumber(0, 0));  // |0⟩
        double f = QuantumFidelity.BlochFidelity(b, b);
        Assert.AreEqual(1.0, f, Tolerance);
    }

    [TestMethod]
    public void BlochFidelity_OppositeVectors_ReturnsZero()
    {
        var north = new BlochVector(0, 0, 1);   // |0⟩
        var south = new BlochVector(0, 0, -1);  // |1⟩
        double f = QuantumFidelity.BlochFidelity(north, south);
        Assert.AreEqual(0.0, f, Tolerance);
    }

    [TestMethod]
    public void BlochFidelity_MatchesStateFidelity()
    {
        var sim = new QuantumSimulator();
        var s1 = sim.Run(QuantumCircuitBuilder.New(1).Ry(0, Math.PI / 3).Build());
        var s2 = sim.Run(QuantumCircuitBuilder.New(1).Ry(0, Math.PI / 5).Build());

        double fState = QuantumFidelity.Fidelity(s1, s2);
        double fBloch = QuantumFidelity.BlochFidelity(s1.GetBlochVector(), s2.GetBlochVector());

        Assert.AreEqual(fState, fBloch, 1e-8);
    }

    [TestMethod]
    public void Fidelity_MultiQubit_BellState()
    {
        var sim = new QuantumSimulator();
        var bell = sim.Run(QuantumCircuitBuilder.New(2).H(0).CNOT(0, 1).Build());
        var zero = sim.Run(QuantumCircuitBuilder.New(2).Build());

        double f = QuantumFidelity.Fidelity(bell, zero);
        // |⟨00|Bell⟩|² = |1/√2|² = 0.5
        Assert.AreEqual(0.5, f, Tolerance);
    }
}
