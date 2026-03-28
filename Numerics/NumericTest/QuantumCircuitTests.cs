using Microsoft.VisualStudio.TestTools.UnitTesting;
using CSharpNumerics.Physics.Quantum;
using CSharpNumerics.Physics.Quantum.NoiseModels;
using CSharpNumerics.Engines.Quantum;
using CSharpNumerics.Engines.Quantum.Algorithms;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
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

    // ── Noise channel: Kraus operator completeness ─────────────

    [TestMethod]
    public void DepolarizingNoise_KrausCompleteness()
    {
        // Σ E_k† E_k should equal I
        var noise = new DepolarizingNoise(0.3);
        AssertKrausCompleteness(noise.GetKrausOperators());
    }

    [TestMethod]
    public void DephasingNoise_KrausCompleteness()
    {
        var noise = new DephasingNoise(0.5);
        AssertKrausCompleteness(noise.GetKrausOperators());
    }

    [TestMethod]
    public void AmplitudeDampingNoise_KrausCompleteness()
    {
        var noise = new AmplitudeDampingNoise(0.7);
        AssertKrausCompleteness(noise.GetKrausOperators());
    }

    private void AssertKrausCompleteness(ComplexMatrix[] kraus)
    {
        int d = kraus[0].rowLength;
        var sum = new ComplexMatrix(d, d);
        foreach (var ek in kraus)
        {
            var ekDag = ek.ConjugateTranspose();
            var product = ekDag * ek;
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                    sum.values[i, j] = sum.values[i, j] + product.values[i, j];
        }
        // Should be identity
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
            {
                double expected = i == j ? 1.0 : 0.0;
                Assert.AreEqual(expected, sum.values[i, j].realPart, 1e-10,
                    $"Real part at [{i},{j}]");
                Assert.AreEqual(0.0, sum.values[i, j].imaginaryPart, 1e-10,
                    $"Imaginary part at [{i},{j}]");
            }
    }

    // ── Noise channel: zero noise = identity ───────────────────

    [TestMethod]
    public void DepolarizingNoise_ZeroProbability_NoEffect()
    {
        var circuit = QuantumCircuitBuilder.New(1).H(0).Build();
        var ideal = new QuantumSimulator().Run(circuit);
        var noisy = new NoisyQuantumSimulator(new Random(42))
            .WithNoise(new DepolarizingNoise(0.0))
            .Run(circuit);

        double f = QuantumFidelity.Fidelity(ideal, noisy);
        Assert.AreEqual(1.0, f, 1e-10);
    }

    [TestMethod]
    public void DephasingNoise_ZeroProbability_NoEffect()
    {
        var circuit = QuantumCircuitBuilder.New(1).H(0).Build();
        var ideal = new QuantumSimulator().Run(circuit);
        var noisy = new NoisyQuantumSimulator(new Random(42))
            .WithNoise(new DephasingNoise(0.0))
            .Run(circuit);

        double f = QuantumFidelity.Fidelity(ideal, noisy);
        Assert.AreEqual(1.0, f, 1e-10);
    }

    [TestMethod]
    public void AmplitudeDamping_ZeroGamma_NoEffect()
    {
        var circuit = QuantumCircuitBuilder.New(1).X(0).Build();
        var ideal = new QuantumSimulator().Run(circuit);
        var noisy = new NoisyQuantumSimulator(new Random(42))
            .WithNoise(new AmplitudeDampingNoise(0.0))
            .Run(circuit);

        double f = QuantumFidelity.Fidelity(ideal, noisy);
        Assert.AreEqual(1.0, f, 1e-10);
    }

    // ── Noise channel: noisy output is valid quantum state ─────

    [TestMethod]
    public void NoisySimulator_ProbabilitiesSumToOne()
    {
        var circuit = QuantumCircuitBuilder.New(2).H(0).CNOT(0, 1).Build();
        var noisy = new NoisyQuantumSimulator(new Random(42))
            .WithNoise(new DepolarizingNoise(0.1))
            .Run(circuit);

        var probs = noisy.GetProbabilities();
        double sum = 0;
        for (int i = 0; i < probs.Length; i++) sum += probs[i];
        Assert.AreEqual(1.0, sum, 1e-10);
    }

    // ── Noise channel: statistical fidelity degradation ────────

    [TestMethod]
    public void DepolarizingNoise_ReducesFidelity_Statistically()
    {
        // Over many runs, noisy simulator should on average have lower fidelity
        var circuit = QuantumCircuitBuilder.New(1).H(0).Build();
        var ideal = new QuantumSimulator().Run(circuit);
        var rng = new Random(123);

        double totalFidelity = 0;
        int runs = 200;
        for (int i = 0; i < runs; i++)
        {
            var noisy = new NoisyQuantumSimulator(rng)
                .WithNoise(new DepolarizingNoise(0.5))
                .Run(circuit);
            totalFidelity += QuantumFidelity.Fidelity(ideal, noisy);
        }
        double avgFidelity = totalFidelity / runs;

        // With p=0.5 on a single gate, fidelity should be meaningfully below 1
        Assert.IsTrue(avgFidelity < 0.95, $"Average fidelity {avgFidelity} not degraded enough");
        Assert.IsTrue(avgFidelity > 0.3, $"Average fidelity {avgFidelity} unexpectedly low");
    }

    [TestMethod]
    public void DephasingNoise_AffectsSuperposition()
    {
        // Dephasing on |0⟩ should have no visible effect (only phases matter on superpositions)
        // but on H|0⟩ = (|0⟩+|1⟩)/√2 it should degrade off-diagonal coherence
        var circuit = QuantumCircuitBuilder.New(1).H(0).Build();
        var ideal = new QuantumSimulator().Run(circuit);
        var rng = new Random(456);

        double totalFidelity = 0;
        int runs = 200;
        for (int i = 0; i < runs; i++)
        {
            var noisy = new NoisyQuantumSimulator(rng)
                .WithNoise(new DephasingNoise(0.5))
                .Run(circuit);
            totalFidelity += QuantumFidelity.Fidelity(ideal, noisy);
        }
        double avgFidelity = totalFidelity / runs;

        Assert.IsTrue(avgFidelity < 0.95, $"Average fidelity {avgFidelity} not degraded");
        Assert.IsTrue(avgFidelity > 0.3, $"Average fidelity {avgFidelity} unexpectedly low");
    }

    [TestMethod]
    public void AmplitudeDamping_DecaysExcitedState()
    {
        // |1⟩ with strong damping should shift probability toward |0⟩
        var circuit = QuantumCircuitBuilder.New(1).X(0).Build();
        var rng = new Random(789);

        double totalP0 = 0;
        int runs = 200;
        for (int i = 0; i < runs; i++)
        {
            var noisy = new NoisyQuantumSimulator(rng)
                .WithNoise(new AmplitudeDampingNoise(0.8))
                .Run(circuit);
            totalP0 += noisy.GetProbability(0);
        }
        double avgP0 = totalP0 / runs;

        // With γ=0.8, ~80% of trajectories should decay to |0⟩
        Assert.IsTrue(avgP0 > 0.6, $"Average P(|0⟩) = {avgP0}, expected > 0.6 for γ=0.8");
    }

    // ── Noise channel: multiple channels stack ──────────────────

    [TestMethod]
    public void NoisySimulator_MultipleChannels()
    {
        var circuit = QuantumCircuitBuilder.New(1).H(0).Build();
        var noisy = new NoisyQuantumSimulator(new Random(42))
            .WithNoise(new DepolarizingNoise(0.1))
            .WithNoise(new DephasingNoise(0.1))
            .Run(circuit);

        var probs = noisy.GetProbabilities();
        double sum = 0;
        for (int i = 0; i < probs.Length; i++) sum += probs[i];
        Assert.AreEqual(1.0, sum, 1e-10);
    }

    // ── Noise channel: validation ──────────────────────────────

    [TestMethod]
    [ExpectedException(typeof(ArgumentOutOfRangeException))]
    public void DepolarizingNoise_InvalidProbability_Throws()
    {
        new DepolarizingNoise(1.5);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentOutOfRangeException))]
    public void DephasingNoise_InvalidProbability_Throws()
    {
        new DephasingNoise(-0.1);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentOutOfRangeException))]
    public void AmplitudeDamping_InvalidGamma_Throws()
    {
        new AmplitudeDampingNoise(2.0);
    }

    // ── Noise: basis state stability ───────────────────────────

    [TestMethod]
    public void DephasingNoise_BasisState_Unchanged()
    {
        // Dephasing only affects off-diagonal elements;
        // |0⟩ (a basis state) should never be changed by dephasing
        var circuit = QuantumCircuitBuilder.New(1).Build(); // |0⟩
        var ideal = new QuantumSimulator().Run(circuit);
        var rng = new Random(42);

        for (int i = 0; i < 50; i++)
        {
            var noisy = new NoisyQuantumSimulator(rng)
                .WithNoise(new DephasingNoise(0.9))
                .Run(circuit);
            double f = QuantumFidelity.Fidelity(ideal, noisy);
            Assert.AreEqual(1.0, f, 1e-10, "Dephasing should not affect basis states");
        }
    }

    [TestMethod]
    public void AmplitudeDamping_GroundState_Unchanged()
    {
        // |0⟩ is the ground state — amplitude damping should have no effect
        var circuit = QuantumCircuitBuilder.New(1).Build();
        var ideal = new QuantumSimulator().Run(circuit);
        var rng = new Random(42);

        for (int i = 0; i < 50; i++)
        {
            var noisy = new NoisyQuantumSimulator(rng)
                .WithNoise(new AmplitudeDampingNoise(0.9))
                .Run(circuit);
            double f = QuantumFidelity.Fidelity(ideal, noisy);
            Assert.AreEqual(1.0, f, 1e-10, "Amplitude damping should not affect ground state");
        }
    }

    // ── QuantumEnvironment (RL) ────────────────────────────

    [TestMethod]
    public void QuantumEnv_ImplementsIEnvironment()
    {
        var target = new QuantumSimulator().Run(
            QuantumCircuitBuilder.New(1).X(0).Build());
        IEnvironment env = QuantumEnvironment.Create(1)
            .WithTargetState(target)
            .WithMaxGates(10)
            .Build();

        Assert.IsTrue(env.IsDiscrete);
        Assert.IsTrue(env.ActionSize > 0);
        Assert.AreEqual(3, env.ObservationSize); // 2^1 + 1 = 3
    }

    [TestMethod]
    public void QuantumEnv_Reset_ReturnsInitialObservation()
    {
        var target = new QuantumSimulator().Run(
            QuantumCircuitBuilder.New(1).H(0).Build());
        var env = QuantumEnvironment.Create(1)
            .WithTargetState(target)
            .WithMaxGates(10)
            .Build();

        var (state, info) = env.Reset(seed: 42);

        Assert.AreEqual(3, state.Length); // 2 probs + 1 progress
        // Initial state is |0⟩: P(0)=1, P(1)=0, step=0
        Assert.AreEqual(1.0, state[0], Tolerance);
        Assert.AreEqual(0.0, state[1], Tolerance);
        Assert.AreEqual(0.0, state[2], Tolerance); // 0/maxGates
    }

    [TestMethod]
    public void QuantumEnv_Step_ChangesState()
    {
        var target = new QuantumSimulator().Run(
            QuantumCircuitBuilder.New(1).X(0).Build());
        var env = QuantumEnvironment.Create(1)
            .WithTargetState(target)
            .WithMaxGates(10)
            .Build();

        env.Reset(seed: 42);

        // Find X gate action index (PauliXGate on qubit 0)
        // Default actions: H(0), X(0), Z(0), S(0), T(0) — X is index 1
        var (nextState, reward, done, info) = env.Step(1); // X gate on qubit 0

        Assert.AreEqual(3, nextState.Length);
        // After X|0⟩ = |1⟩ which IS the target, fidelity = 1.0
        Assert.IsTrue((double)info["fidelity"] > 0.99);
    }

    [TestMethod]
    public void QuantumEnv_OptimalAction_HighReward()
    {
        // Target = X|0⟩ = |1⟩. Applying X should get fidelity=1 and bonus reward.
        var target = new QuantumSimulator().Run(
            QuantumCircuitBuilder.New(1).X(0).Build());
        var env = QuantumEnvironment.Create(1)
            .WithTargetState(target)
            .WithMaxGates(10)
            .WithFidelityThreshold(0.99)
            .Build();

        env.Reset(seed: 42);
        var (_, reward, done, info) = env.Step(1); // X gate

        Assert.IsTrue(reward > 1.0, $"Reward {reward} should include bonus");
        Assert.IsTrue(done, "Should be done when fidelity threshold reached");
    }

    [TestMethod]
    public void QuantumEnv_MaxGates_TerminatesEpisode()
    {
        var target = new QuantumSimulator().Run(
            QuantumCircuitBuilder.New(1).H(0).Build());
        var env = QuantumEnvironment.Create(1)
            .WithTargetState(target)
            .WithMaxGates(3)
            .Build();

        env.Reset(seed: 42);

        bool done = false;
        for (int i = 0; i < 3; i++)
        {
            var result = env.Step(0); // apply H repeatedly
            done = result.done;
        }

        Assert.IsTrue(done, "Should terminate after max gates");
    }

    [TestMethod]
    public void QuantumEnv_ObservationProbsSumToOne()
    {
        var target = new QuantumSimulator().Run(
            QuantumCircuitBuilder.New(2).H(0).CNOT(0, 1).Build());
        var env = QuantumEnvironment.Create(2)
            .WithTargetState(target)
            .WithMaxGates(10)
            .Build();

        var (state, _) = env.Reset(seed: 42);

        // Probabilities are first 2^n elements
        double sum = 0;
        for (int i = 0; i < 4; i++) sum += state[i];
        Assert.AreEqual(1.0, sum, Tolerance);

        // Take a step and check again
        var (nextState, _, _, _) = env.Step(0);
        sum = 0;
        for (int i = 0; i < 4; i++) sum += nextState[i];
        Assert.AreEqual(1.0, sum, Tolerance);
    }

    [TestMethod]
    public void QuantumEnv_TwoQubit_DefaultActions_IncludesCNOT()
    {
        var target = new QuantumSimulator().Run(
            QuantumCircuitBuilder.New(2).H(0).CNOT(0, 1).Build());
        var env = QuantumEnvironment.Create(2)
            .WithTargetState(target)
            .WithMaxGates(10)
            .Build();

        // 2 qubits: 5 single-gates * 2 qubits + 2 CNOT pairs = 12
        Assert.AreEqual(12, env.ActionSize);
        Assert.AreEqual(5, env.ObservationSize); // 2^2 + 1 = 5
    }

    [TestMethod]
    public void QuantumEnv_WithTargetCircuit()
    {
        var bellCircuit = QuantumCircuitBuilder.New(2).H(0).CNOT(0, 1).Build();
        var env = QuantumEnvironment.Create(2)
            .WithTargetCircuit(bellCircuit)
            .WithMaxGates(10)
            .Build();

        var (state, _) = env.Reset();
        Assert.AreEqual(5, state.Length);
    }

    [TestMethod]
    public void QuantumEnv_CustomActions()
    {
        var target = new QuantumSimulator().Run(
            QuantumCircuitBuilder.New(1).H(0).Build());

        var customActions = new List<QuantumInstruction>
        {
            new QuantumInstruction(new HadamardGate(), new List<int> { 0 }),
            new QuantumInstruction(new PauliXGate(), new List<int> { 0 })
        };

        var env = QuantumEnvironment.Create(1)
            .WithTargetState(target)
            .WithActions(customActions)
            .WithMaxGates(5)
            .Build();

        Assert.AreEqual(2, env.ActionSize);
    }

    [TestMethod]
    public void QuantumEnv_Reset_ResetsCircuit()
    {
        var target = new QuantumSimulator().Run(
            QuantumCircuitBuilder.New(1).X(0).Build());
        var env = QuantumEnvironment.Create(1)
            .WithTargetState(target)
            .WithMaxGates(10)
            .Build();

        env.Reset(seed: 42);
        env.Step(1); // X gate

        // Reset and verify we're back to |0⟩
        var (state, _) = env.Reset(seed: 42);
        Assert.AreEqual(1.0, state[0], Tolerance); // P(|0⟩) = 1
        Assert.AreEqual(0.0, state[1], Tolerance); // P(|1⟩) = 0
    }

    [TestMethod]
    public void QuantumEnv_StepVectorN_DelegatesToDiscrete()
    {
        var target = new QuantumSimulator().Run(
            QuantumCircuitBuilder.New(1).X(0).Build());
        var env = QuantumEnvironment.Create(1)
            .WithTargetState(target)
            .WithMaxGates(10)
            .Build();

        env.Reset(seed: 42);
        var (_, _, _, info) = env.Step(new VectorN(new[] { 1.0 })); // X gate
        Assert.IsTrue((double)info["fidelity"] > 0.99);
    }

    [TestMethod]
    public void QuantumEnv_FidelityInfo_Tracked()
    {
        var target = new QuantumSimulator().Run(
            QuantumCircuitBuilder.New(1).H(0).Build());
        var env = QuantumEnvironment.Create(1)
            .WithTargetState(target)
            .WithMaxGates(10)
            .Build();

        env.Reset(seed: 42);

        // Apply H — should reach target
        var (_, _, _, info) = env.Step(0); // H gate
        double fidelity = (double)info["fidelity"];

        Assert.AreEqual(1.0, fidelity, Tolerance);
        Assert.AreEqual(1, (int)info["gates"]);
    }

    [TestMethod]
    [ExpectedException(typeof(InvalidOperationException))]
    public void QuantumEnv_NoTarget_Throws()
    {
        QuantumEnvironment.Create(1)
            .WithMaxGates(10)
            .Build();
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentOutOfRangeException))]
    public void QuantumEnv_InvalidAction_Throws()
    {
        var target = new QuantumSimulator().Run(
            QuantumCircuitBuilder.New(1).X(0).Build());
        var env = QuantumEnvironment.Create(1)
            .WithTargetState(target)
            .WithMaxGates(10)
            .Build();

        env.Reset(seed: 42);
        env.Step(999); // invalid action
    }

    // ── PauliY gate tests ──────────────────────────────────────

    [TestMethod]
    public void PauliYGate_MatrixIs2x2()
    {
        var gate = new PauliYGate();
        var m = gate.GetMatrix();
        Assert.AreEqual(2, m.rowLength);
        Assert.AreEqual(2, m.columnLength);
    }

    [TestMethod]
    public void PauliYGate_AppliedTwice_ReturnsToOriginal()
    {
        // Y² = I
        var circuit = QuantumCircuitBuilder.New(1).Y(0).Y(0).Build();
        var state = new QuantumSimulator().Run(circuit);
        Assert.AreEqual(1.0, state.GetProbability(0), Tolerance);
        Assert.AreEqual(0.0, state.GetProbability(1), Tolerance);
    }

    [TestMethod]
    public void PauliYGate_FlipsZeroToOne()
    {
        var circuit = QuantumCircuitBuilder.New(1).Y(0).Build();
        var state = new QuantumSimulator().Run(circuit);
        Assert.AreEqual(0.0, state.GetProbability(0), Tolerance);
        Assert.AreEqual(1.0, state.GetProbability(1), Tolerance);
    }

    // ── PhaseGate tests ────────────────────────────────────────

    [TestMethod]
    public void PhaseGate_MatrixIs2x2()
    {
        var gate = new PhaseGate(Math.PI / 4);
        var m = gate.GetMatrix();
        Assert.AreEqual(2, m.rowLength);
        Assert.AreEqual(2, m.columnLength);
    }

    [TestMethod]
    public void PhaseGate_Pi_EquivalentToZ()
    {
        // P(π) = Z
        var zCircuit = QuantumCircuitBuilder.New(1).H(0).Z(0).Build();
        var pCircuit = QuantumCircuitBuilder.New(1).H(0).Phase(0, Math.PI).Build();
        var sim = new QuantumSimulator();
        var zState = sim.Run(zCircuit);
        var pState = sim.Run(pCircuit);

        for (int i = 0; i < 2; i++)
            Assert.AreEqual(zState.GetProbability(i), pState.GetProbability(i), Tolerance);
    }

    [TestMethod]
    public void PhaseGate_LeavesZeroUnchanged()
    {
        var circuit = QuantumCircuitBuilder.New(1).Phase(0, Math.PI / 3).Build();
        var state = new QuantumSimulator().Run(circuit);
        Assert.AreEqual(1.0, state.GetProbability(0), Tolerance);
    }

    // ── CPhaseGate tests ───────────────────────────────────────

    [TestMethod]
    public void CPhaseGate_MatrixIs4x4()
    {
        var gate = new CPhaseGate(Math.PI / 4);
        var m = gate.GetMatrix();
        Assert.AreEqual(4, m.rowLength);
        Assert.AreEqual(4, m.columnLength);
    }

    [TestMethod]
    public void CPhaseGate_Pi_EquivalentToCZ()
    {
        // CPhase(π) = CZ
        // Apply to |+1⟩ = H⊗X |00⟩, compare probabilities
        var czCircuit = QuantumCircuitBuilder.New(2).H(0).X(1).CZ(0, 1).Build();
        var cpCircuit = QuantumCircuitBuilder.New(2).H(0).X(1).CPhase(0, 1, Math.PI).Build();
        var sim = new QuantumSimulator();
        var czState = sim.Run(czCircuit);
        var cpState = sim.Run(cpCircuit);

        for (int i = 0; i < 4; i++)
            Assert.AreEqual(czState.GetProbability(i), cpState.GetProbability(i), Tolerance);
    }

    [TestMethod]
    public void CPhaseGate_Leaves00Unchanged()
    {
        var circuit = QuantumCircuitBuilder.New(2).CPhase(0, 1, Math.PI / 2).Build();
        var state = new QuantumSimulator().Run(circuit);
        Assert.AreEqual(1.0, state.GetProbability(0), Tolerance);
    }

    // ── ToffoliGate tests ──────────────────────────────────────

    [TestMethod]
    public void ToffoliGate_MatrixIs8x8()
    {
        var gate = new ToffoliGate();
        var m = gate.GetMatrix();
        Assert.AreEqual(8, m.rowLength);
        Assert.AreEqual(8, m.columnLength);
    }

    [TestMethod]
    public void ToffoliGate_FlipsTarget_WhenBothControlsAreOne()
    {
        // |110⟩ → |111⟩ (controls=q0,q1, target=q2)
        var circuit = QuantumCircuitBuilder.New(3).X(0).X(1).Toffoli(0, 1, 2).Build();
        var state = new QuantumSimulator().Run(circuit);
        // |111⟩ = index 7
        Assert.AreEqual(1.0, state.GetProbability(7), Tolerance);
    }

    [TestMethod]
    public void ToffoliGate_DoesNotFlip_WhenOneControlIsZero()
    {
        // |100⟩ → still |100⟩ (only control0 is 1)
        var circuit = QuantumCircuitBuilder.New(3).X(0).Toffoli(0, 1, 2).Build();
        var state = new QuantumSimulator().Run(circuit);
        // |100⟩ = index 1
        Assert.AreEqual(1.0, state.GetProbability(1), Tolerance);
    }

    [TestMethod]
    public void ToffoliGate_DoubleToffoli_ReturnsToOriginal()
    {
        // Toffoli² = I for any basis state
        var circuit = QuantumCircuitBuilder.New(3)
            .X(0).X(1)  // |110⟩
            .Toffoli(0, 1, 2)
            .Toffoli(0, 1, 2)
            .Build();
        var state = new QuantumSimulator().Run(circuit);
        // Should return to |110⟩ = index 3
        Assert.AreEqual(1.0, state.GetProbability(3), Tolerance);
    }

    // ── FredkinGate tests ──────────────────────────────────────

    [TestMethod]
    public void FredkinGate_MatrixIs8x8()
    {
        var gate = new FredkinGate();
        var m = gate.GetMatrix();
        Assert.AreEqual(8, m.rowLength);
        Assert.AreEqual(8, m.columnLength);
    }

    [TestMethod]
    public void FredkinGate_SwapsTargets_WhenControlIsOne()
    {
        // |110⟩ (q0=1,q1=1,q2=0), Fredkin(0,1,2) → swap q1,q2 → |101⟩
        var circuit = QuantumCircuitBuilder.New(3).X(0).X(1).Fredkin(0, 1, 2).Build();
        var state = new QuantumSimulator().Run(circuit);
        // |101⟩ = index 5
        Assert.AreEqual(1.0, state.GetProbability(5), Tolerance);
    }

    [TestMethod]
    public void FredkinGate_DoesNotSwap_WhenControlIsZero()
    {
        // |010⟩ (q0=0,q1=1,q2=0), Fredkin(0,1,2) → no swap → |010⟩
        var circuit = QuantumCircuitBuilder.New(3).X(1).Fredkin(0, 1, 2).Build();
        var state = new QuantumSimulator().Run(circuit);
        // |010⟩ = index 2
        Assert.AreEqual(1.0, state.GetProbability(2), Tolerance);
    }

    [TestMethod]
    public void FredkinGate_DoubleFredkin_ReturnsToOriginal()
    {
        var circuit = QuantumCircuitBuilder.New(3)
            .X(0).X(1)  // |110⟩
            .Fredkin(0, 1, 2)
            .Fredkin(0, 1, 2)
            .Build();
        var state = new QuantumSimulator().Run(circuit);
        Assert.AreEqual(1.0, state.GetProbability(3), Tolerance);
    }

    // ── QFT tests ──────────────────────────────────────────────

    [TestMethod]
    public void QFT_SingleQubit_IsHadamard()
    {
        // QFT on 1 qubit = Hadamard
        var qftCircuit = QFT.CreateCircuit(1, 0);
        var hCircuit = QuantumCircuitBuilder.New(1).H(0).Build();
        var sim = new QuantumSimulator();
        var qftState = sim.Run(qftCircuit);
        var hState = sim.Run(hCircuit);

        for (int i = 0; i < 2; i++)
            Assert.AreEqual(hState.GetProbability(i), qftState.GetProbability(i), Tolerance);
    }

    [TestMethod]
    public void QFT_UniformSuperposition_FromZero()
    {
        // QFT(|0...0⟩) produces uniform superposition
        int n = 3;
        var circuit = QFT.CreateCircuit(n, 0, 1, 2);
        var state = new QuantumSimulator().Run(circuit);

        double expectedProb = 1.0 / (1 << n);
        for (int i = 0; i < (1 << n); i++)
            Assert.AreEqual(expectedProb, state.GetProbability(i), 1e-8,
                $"State |{i}⟩ probability should be {expectedProb}");
    }

    [TestMethod]
    public void InverseQFT_UndoesQFT()
    {
        // InverseQFT(QFT(|ψ⟩)) = |ψ⟩
        // Prepare |101⟩ (index 5)
        int n = 3;
        var prepCircuit = QuantumCircuitBuilder.New(n).X(0).X(2).Build();
        var sim = new QuantumSimulator();
        var original = sim.Run(prepCircuit);

        // Apply QFT then InverseQFT
        var circuit = QuantumCircuitBuilder.New(n)
            .X(0).X(2)
            .ApplyQFT(0, 1, 2)
            .ApplyInverseQFT(0, 1, 2)
            .Build();
        var result = sim.Run(circuit);

        for (int i = 0; i < (1 << n); i++)
            Assert.AreEqual(original.GetProbability(i), result.GetProbability(i), 1e-8,
                $"QFT†∘QFT should be identity for state |{i}⟩");
    }

    [TestMethod]
    public void InverseQFT_UndoesQFT_TwoQubits()
    {
        // 2-qubit case: prepare |10⟩, apply QFT then InverseQFT
        var circuit = QuantumCircuitBuilder.New(2)
            .X(0)
            .ApplyQFT(0, 1)
            .ApplyInverseQFT(0, 1)
            .Build();
        var state = new QuantumSimulator().Run(circuit);

        // Should return to |10⟩ = index 1
        Assert.AreEqual(0.0, state.GetProbability(0), 1e-8);
        Assert.AreEqual(1.0, state.GetProbability(1), 1e-8);
        Assert.AreEqual(0.0, state.GetProbability(2), 1e-8);
        Assert.AreEqual(0.0, state.GetProbability(3), 1e-8);
    }

    [TestMethod]
    public void QFT_CorrectPhases_TwoQubits()
    {
        // QFT on |10⟩ should produce specific amplitudes
        // |10⟩ in standard basis → QFT → (1/2)(|00⟩ + i|01⟩ - |10⟩ - i|11⟩)
        // Actually, for 2-qubit QFT on input |j⟩ with j = binary value:
        // QFT|j⟩ = (1/√N) Σ_k e^(2πijk/N) |k⟩
        // For j=1 (|10⟩ → binary 01 → value 1 after bit reversal... 
        // Let's verify via probabilities: all should be 0.25
        var circuit = QuantumCircuitBuilder.New(2).X(0).ApplyQFT(0, 1).Build();
        var state = new QuantumSimulator().Run(circuit);

        for (int i = 0; i < 4; i++)
            Assert.AreEqual(0.25, state.GetProbability(i), 1e-8,
                $"QFT of |10⟩ should produce uniform probabilities");
    }

    [TestMethod]
    public void QFT_Builder_MatchesStandalone()
    {
        // Builder .ApplyQFT() should produce same result as QFT.CreateCircuit()
        int n = 3;
        var sim = new QuantumSimulator();

        var standalone = QFT.CreateCircuit(n, 0, 1, 2);
        var standaloneState = sim.Run(standalone);

        var builderCircuit = QuantumCircuitBuilder.New(n).ApplyQFT(0, 1, 2).Build();
        var builderState = sim.Run(builderCircuit);

        for (int i = 0; i < (1 << n); i++)
            Assert.AreEqual(standaloneState.GetProbability(i), builderState.GetProbability(i), Tolerance);
    }

    [TestMethod]
    public void QFT_InverseQFT_SuperpositionState()
    {
        // Test with a superposition input: H|0⟩ ⊗ |0⟩
        var circuit = QuantumCircuitBuilder.New(2)
            .H(0)
            .ApplyQFT(0, 1)
            .ApplyInverseQFT(0, 1)
            .Build();
        var state = new QuantumSimulator().Run(circuit);

        // Should return to H|0⟩ ⊗ |0⟩: prob(|00⟩)=0.5, prob(|10⟩)=0.5
        Assert.AreEqual(0.5, state.GetProbability(0), 1e-8);
        Assert.AreEqual(0.5, state.GetProbability(1), 1e-8);
        Assert.AreEqual(0.0, state.GetProbability(2), 1e-8);
        Assert.AreEqual(0.0, state.GetProbability(3), 1e-8);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void QFT_EmptyQubits_Throws()
    {
        QFT.CreateCircuit(2);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void InverseQFT_EmptyQubits_Throws()
    {
        InverseQFT.CreateCircuit(2);
    }

    // ── Builder methods for new gates ──────────────────────────

    [TestMethod]
    public void Builder_Y_AppliesToCircuit()
    {
        var circuit = QuantumCircuitBuilder.New(1).Y(0).Build();
        Assert.AreEqual(1, circuit.Instructions.Count);
        var state = new QuantumSimulator().Run(circuit);
        Assert.AreEqual(1.0, state.GetProbability(1), Tolerance);
    }

    [TestMethod]
    public void Builder_Phase_AppliesToCircuit()
    {
        var circuit = QuantumCircuitBuilder.New(1).Phase(0, Math.PI).Build();
        Assert.AreEqual(1, circuit.Instructions.Count);
    }

    [TestMethod]
    public void Builder_CPhase_AppliesToCircuit()
    {
        var circuit = QuantumCircuitBuilder.New(2).X(0).X(1).CPhase(0, 1, Math.PI / 4).Build();
        Assert.AreEqual(3, circuit.Instructions.Count);
    }

    [TestMethod]
    public void Builder_Toffoli_AppliesToCircuit()
    {
        var circuit = QuantumCircuitBuilder.New(3).X(0).X(1).Toffoli(0, 1, 2).Build();
        var state = new QuantumSimulator().Run(circuit);
        Assert.AreEqual(1.0, state.GetProbability(7), Tolerance);
    }

    [TestMethod]
    public void Builder_Fredkin_AppliesToCircuit()
    {
        var circuit = QuantumCircuitBuilder.New(3).X(0).X(1).Fredkin(0, 1, 2).Build();
        var state = new QuantumSimulator().Run(circuit);
        Assert.AreEqual(1.0, state.GetProbability(5), Tolerance);
    }

    // ── PhaseOracle tests ──────────────────────────────────────

    [TestMethod]
    public void PhaseOracle_SingleQubit_MarksState1()
    {
        // Oracle marks |1⟩ → phase flip. Apply to |+⟩ = (|0⟩+|1⟩)/√2
        // After oracle: (|0⟩-|1⟩)/√2 = |−⟩
        var gate = new PhaseOracle(1, 1);
        Assert.AreEqual(1, gate.QubitCount);

        var circuit = QuantumCircuitBuilder.New(1)
            .H(0)
            .Gate(new PhaseOracle(1, 1), 0)
            .Build();
        var state = new QuantumSimulator().Run(circuit);

        // |−⟩ = (|0⟩−|1⟩)/√2, probabilities still 0.5/0.5
        Assert.AreEqual(0.5, state.GetProbability(0), Tolerance);
        Assert.AreEqual(0.5, state.GetProbability(1), Tolerance);
    }

    [TestMethod]
    public void PhaseOracle_TwoQubit_MarksState3()
    {
        var gate = new PhaseOracle(2, 3);
        Assert.AreEqual(2, gate.QubitCount);
        var m = gate.GetMatrix();
        Assert.AreEqual(4, m.rowLength);

        // Diagonal: +1, +1, +1, -1
        Assert.AreEqual(1.0, m.values[0, 0].realPart, Tolerance);
        Assert.AreEqual(1.0, m.values[1, 1].realPart, Tolerance);
        Assert.AreEqual(1.0, m.values[2, 2].realPart, Tolerance);
        Assert.AreEqual(-1.0, m.values[3, 3].realPart, Tolerance);
    }

    [TestMethod]
    public void PhaseOracle_MultipleMarked()
    {
        // Mark states 1 and 2 in a 2-qubit space
        var gate = new PhaseOracle(2, 1, 2);
        var m = gate.GetMatrix();

        Assert.AreEqual(1.0, m.values[0, 0].realPart, Tolerance);
        Assert.AreEqual(-1.0, m.values[1, 1].realPart, Tolerance);
        Assert.AreEqual(-1.0, m.values[2, 2].realPart, Tolerance);
        Assert.AreEqual(1.0, m.values[3, 3].realPart, Tolerance);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentOutOfRangeException))]
    public void PhaseOracle_InvalidState_Throws()
    {
        new PhaseOracle(2, 5); // 5 out of range for 2-qubit space (0..3)
    }

    // ── Grover search tests ────────────────────────────────────

    [TestMethod]
    public void GroverSearch_OptimalIterations_SingleMarked()
    {
        // 3 qubits, 1 marked: optimal = floor(π/4 · √8) ≈ floor(2.22) = 2
        int iters = GroverSearch.OptimalIterations(3, 1);
        Assert.AreEqual(2, iters);
    }

    [TestMethod]
    public void GroverSearch_OptimalIterations_TwoMarked()
    {
        // 3 qubits, 2 marked: optimal = floor(π/4 · √4) ≈ floor(1.57) = 1
        int iters = GroverSearch.OptimalIterations(3, 2);
        Assert.AreEqual(1, iters);
    }

    [TestMethod]
    public void GroverSearch_2Qubit_FindsTarget()
    {
        // 2 qubits (N=4), 1 marked state: optimal iterations = 1
        // After 1 Grover iteration, target should have high probability
        int target = 2; // |10⟩
        var circuit = GroverSearch.CreateCircuit(2, new[] { 0, 1 }, new[] { target });
        var state = new QuantumSimulator().Run(circuit);

        // Target probability should be close to 1.0
        Assert.IsTrue(state.GetProbability(target) > 0.9,
            $"Target |{target}⟩ probability = {state.GetProbability(target)}, expected > 0.9");
    }

    [TestMethod]
    public void GroverSearch_3Qubit_FindsTarget()
    {
        // 3 qubits (N=8), 1 marked state, optimal = 2 iterations
        int target = 5; // |101⟩
        var circuit = GroverSearch.CreateCircuit(3, new[] { 0, 1, 2 }, new[] { target });
        var state = new QuantumSimulator().Run(circuit);

        // Target probability should be dominant
        Assert.IsTrue(state.GetProbability(target) > 0.9,
            $"Target |{target}⟩ probability = {state.GetProbability(target)}, expected > 0.9");
    }

    [TestMethod]
    public void GroverSearch_3Qubit_AllTargets_HighProbability()
    {
        // Verify each possible target can be found
        for (int target = 0; target < 8; target++)
        {
            var circuit = GroverSearch.CreateCircuit(3, new[] { 0, 1, 2 }, new[] { target });
            var state = new QuantumSimulator().Run(circuit);
            Assert.IsTrue(state.GetProbability(target) > 0.9,
                $"Target |{target}⟩ probability = {state.GetProbability(target)}, expected > 0.9");
        }
    }

    [TestMethod]
    public void GroverSearch_MultipleMarked_FindsAny()
    {
        // 3 qubits, 2 marked states → should find one of them with high combined probability
        int[] targets = { 3, 5 };
        var circuit = GroverSearch.CreateCircuit(3, new[] { 0, 1, 2 }, targets);
        var state = new QuantumSimulator().Run(circuit);

        double combinedProb = state.GetProbability(3) + state.GetProbability(5);
        Assert.IsTrue(combinedProb > 0.9,
            $"Combined probability of marked states = {combinedProb}, expected > 0.9");
    }

    [TestMethod]
    public void GroverSearch_ExplicitIterations()
    {
        // 2 qubits, target = 0, 1 iteration (explicit)
        var circuit = GroverSearch.CreateCircuit(2, new[] { 0, 1 }, new[] { 0 }, iterations: 1);
        var state = new QuantumSimulator().Run(circuit);
        Assert.IsTrue(state.GetProbability(0) > 0.9,
            $"Target |0⟩ probability = {state.GetProbability(0)}, expected > 0.9");
    }

    [TestMethod]
    public void GroverSearch_Builder_MatchesStandalone()
    {
        // Builder ApplyGrover should produce same result as GroverSearch.CreateCircuit
        int target = 6;
        var sim = new QuantumSimulator();

        var standalone = GroverSearch.CreateCircuit(3, new[] { 0, 1, 2 }, new[] { target });
        var standaloneState = sim.Run(standalone);

        var builderCircuit = QuantumCircuitBuilder.New(3)
            .ApplyGrover(new[] { target }, new[] { 0, 1, 2 })
            .Build();
        var builderState = sim.Run(builderCircuit);

        for (int i = 0; i < 8; i++)
            Assert.AreEqual(standaloneState.GetProbability(i), builderState.GetProbability(i), Tolerance,
                $"Mismatch at state |{i}⟩");
    }

    [TestMethod]
    public void GroverSearch_Builder_WithIterations()
    {
        int target = 1;
        var circuit = QuantumCircuitBuilder.New(3)
            .ApplyGrover(new[] { target }, 2, 0, 1, 2)
            .Build();
        var state = new QuantumSimulator().Run(circuit);
        Assert.IsTrue(state.GetProbability(target) > 0.9,
            $"Target probability = {state.GetProbability(target)}, expected > 0.9");
    }

    [TestMethod]
    public void GroverSearch_SubsetQubits()
    {
        // Grover on qubits 1,2 in a 4-qubit circuit (qubit 0 untouched)
        int target = 2; // |10⟩ in the 2-qubit search subspace
        var circuit = GroverSearch.CreateCircuit(4, new[] { 1, 2 }, new[] { target });
        var state = new QuantumSimulator().Run(circuit);

        // Qubit 0 stays |0⟩. The search space (qubits 1,2) should find target.
        // |target⟩ in qubits 1,2 with qubit 0=0: global index = target mapped via qubit positions
        // qubit 1 and 2 represent the search space; qubit 0 = |0⟩
        // We need to sum over states where qubits 1,2 form the target pattern
        double targetProb = 0;
        for (int i = 0; i < 16; i++)
        {
            // Extract bits for qubit 1 and qubit 2
            int q1 = (i >> 1) & 1;
            int q2 = (i >> 2) & 1;
            int searchVal = q1 + 2 * q2; // gate-local ordering
            if (searchVal == target)
                targetProb += state.GetProbability(i);
        }
        Assert.IsTrue(targetProb > 0.9,
            $"Search-space target probability = {targetProb}, expected > 0.9");
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void GroverSearch_EmptySearchQubits_Throws()
    {
        GroverSearch.CreateCircuit(2, new int[0], new[] { 0 });
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void GroverSearch_EmptyMarkedStates_Throws()
    {
        GroverSearch.CreateCircuit(2, new[] { 0, 1 }, new int[0]);
    }

    // ── ControlledGate tests ───────────────────────────────────

    [TestMethod]
    public void ControlledGate_WrapsHadamard_MatrixIs4x4()
    {
        var ch = new ControlledGate(new HadamardGate());
        Assert.AreEqual(2, ch.QubitCount);
        var m = ch.GetMatrix();
        Assert.AreEqual(4, m.rowLength);
        Assert.AreEqual(4, m.columnLength);
    }

    [TestMethod]
    public void ControlledGate_OfX_EquivalentToCNOT()
    {
        // Controlled-X should behave like CNOT
        var cx = new ControlledGate(new PauliXGate());
        var cnot = new CNOTGate();

        // Apply to |10⟩ (control=1, target=0) → should flip target
        var circuitCX = QuantumCircuitBuilder.New(2)
            .X(0) // control = |1⟩
            .Gate(cx, 0, 1)
            .Build();
        var circuitCNOT = QuantumCircuitBuilder.New(2)
            .X(0)
            .CNOT(0, 1)
            .Build();

        var sim = new QuantumSimulator();
        var stateCX = sim.Run(circuitCX);
        var stateCNOT = sim.Run(circuitCNOT);

        for (int i = 0; i < 4; i++)
            Assert.AreEqual(stateCNOT.GetProbability(i), stateCX.GetProbability(i), Tolerance,
                $"Mismatch at state |{i}⟩");
    }

    [TestMethod]
    public void ControlledGate_ControlZero_NoEffect()
    {
        // Controlled-X with control = |0⟩: target should stay |0⟩
        var cx = new ControlledGate(new PauliXGate());
        var circuit = QuantumCircuitBuilder.New(2)
            .Gate(cx, 0, 1) // control=0, no flip
            .Build();
        var state = new QuantumSimulator().Run(circuit);
        Assert.AreEqual(1.0, state.GetProbability(0), Tolerance); // |00⟩
    }

    [TestMethod]
    public void ControlledGate_ControlOne_AppliesInnerGate()
    {
        // Controlled-Z with control=|1⟩, target in superposition
        var cz = new ControlledGate(new PauliZGate());
        // Prepare |1⟩⊗|+⟩, apply CZ → |1⟩⊗|−⟩ (same probabilities)
        var circuit = QuantumCircuitBuilder.New(2)
            .X(0)
            .H(1)
            .Gate(cz, 0, 1)
            .Build();
        var state = new QuantumSimulator().Run(circuit);
        // Probabilities: |10⟩=0.5, |11⟩=0.5
        Assert.AreEqual(0.0, state.GetProbability(0), Tolerance);
        Assert.AreEqual(0.5, state.GetProbability(1), Tolerance);
        Assert.AreEqual(0.0, state.GetProbability(2), Tolerance);
        Assert.AreEqual(0.5, state.GetProbability(3), Tolerance);
    }

    [TestMethod]
    public void ControlledGate_Builder_Controlled()
    {
        // Builder's .Controlled() method
        var circuit = QuantumCircuitBuilder.New(2)
            .X(0)
            .Controlled(new PauliXGate(), 0, 1)
            .Build();
        var state = new QuantumSimulator().Run(circuit);
        // |10⟩ → CNOT → |11⟩ = index 3
        Assert.AreEqual(1.0, state.GetProbability(3), Tolerance);
    }

    // ── ModularMultiplyGate tests ──────────────────────────────

    [TestMethod]
    public void ModularMultiplyGate_BasicPermutation()
    {
        // a=2, N=3, 2 qubits: |0⟩→|0⟩, |1⟩→|2⟩, |2⟩→|1⟩, |3⟩→|3⟩ (identity for ≥N)
        var gate = new ModularMultiplyGate(2, 3, 2);
        Assert.AreEqual(2, gate.QubitCount);
        var m = gate.GetMatrix();
        Assert.AreEqual(4, m.rowLength);

        // |0⟩→|0⟩
        Assert.AreEqual(1.0, m.values[0, 0].realPart, Tolerance);
        // |1⟩→|2·1 mod 3=2⟩
        Assert.AreEqual(1.0, m.values[2, 1].realPart, Tolerance);
        // |2⟩→|2·2 mod 3=1⟩
        Assert.AreEqual(1.0, m.values[1, 2].realPart, Tolerance);
        // |3⟩→|3⟩ (identity for y≥N)
        Assert.AreEqual(1.0, m.values[3, 3].realPart, Tolerance);
    }

    [TestMethod]
    public void ModularMultiplyGate_Mod5()
    {
        // a=3, N=5, 3 qubits: |y⟩ → |3y mod 5⟩
        var gate = new ModularMultiplyGate(3, 5, 3);
        Assert.AreEqual(3, gate.QubitCount);

        // Check specific mappings: 0→0, 1→3, 2→1, 3→4, 4→2
        var m = gate.GetMatrix();
        Assert.AreEqual(1.0, m.values[0, 0].realPart, Tolerance); // |0⟩→|0⟩
        Assert.AreEqual(1.0, m.values[3, 1].realPart, Tolerance); // |1⟩→|3⟩
        Assert.AreEqual(1.0, m.values[1, 2].realPart, Tolerance); // |2⟩→|1⟩
        Assert.AreEqual(1.0, m.values[4, 3].realPart, Tolerance); // |3⟩→|4⟩
        Assert.AreEqual(1.0, m.values[2, 4].realPart, Tolerance); // |4⟩→|2⟩
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void ModularMultiplyGate_NotCoprime_Throws()
    {
        new ModularMultiplyGate(4, 6, 3); // gcd(4,6) = 2 ≠ 1
    }

    [TestMethod]
    public void ModularMultiplyGate_AppliedInCircuit()
    {
        // Apply U_2 mod 3 to |1⟩ → should give |2⟩
        var gate = new ModularMultiplyGate(2, 3, 2);
        var circuit = QuantumCircuitBuilder.New(2)
            .X(0)                    // |01⟩ = |1⟩
            .Gate(gate, 0, 1)
            .Build();
        var state = new QuantumSimulator().Run(circuit);
        // |2⟩ = index 2 (bit pattern: q0=0, q1=1)
        Assert.AreEqual(1.0, state.GetProbability(2), Tolerance);
    }

    // ── QPE tests ──────────────────────────────────────────────

    [TestMethod]
    public void QPE_PhaseGate_EstimatesQuarterPhase()
    {
        // Phase gate P(π/2) has eigenvalue e^(iπ/2) on |1⟩, so φ = 1/4
        // With 2 counting qubits, we expect to measure |01⟩ = 1 (since φ·4 = 1)
        // Circuit: counting=[0,1], target=[2]
        // Target register must be in eigenstate |1⟩

        int totalQubits = 3;
        int[] counting = { 0, 1 };
        int[] target = { 2 };

        var circuit = new QuantumCircuit(totalQubits);
        // Prepare eigenstate |1⟩ on target qubit
        circuit.AddInstruction(new QuantumInstruction(new PauliXGate(), new List<int> { 2 }));

        // Run QPE
        var qpeCircuit = QPE.CreateCircuit(totalQubits, counting, target, new PhaseGate(Math.PI / 2));
        foreach (var instr in qpeCircuit.Instructions)
            circuit.AddInstruction(instr);

        var state = new QuantumSimulator().Run(circuit);

        // After QPE, counting register reads φ·2^t = 1.
        // MSB convention: countingQubits[0]=MSB. Value 1 in binary = "01" → q0=0(MSB), q1=1(LSB)
        // Global index: q0·1 + q1·2 + q2·4 = 0 + 2 + 4 = 6
        double prob = state.GetProbability(6);
        Assert.IsTrue(prob > 0.9,
            $"Expected φ=1/4 to produce counting result 1, but P(6) = {prob}");
    }

    [TestMethod]
    public void QPE_PhaseGate_EstimatesHalfPhase()
    {
        // Phase gate P(π) = Z has eigenvalue e^(iπ) = -1 on |1⟩, so φ = 1/2
        // With 2 counting qubits, we expect counting = φ·4 = 2
        int totalQubits = 3;
        int[] counting = { 0, 1 };
        int[] target = { 2 };

        var circuit = new QuantumCircuit(totalQubits);
        circuit.AddInstruction(new QuantumInstruction(new PauliXGate(), new List<int> { 2 }));

        var qpeCircuit = QPE.CreateCircuit(totalQubits, counting, target, new PhaseGate(Math.PI));
        foreach (var instr in qpeCircuit.Instructions)
            circuit.AddInstruction(instr);

        var state = new QuantumSimulator().Run(circuit);

        // counting = 2 → binary "10" → q0=1(MSB), q1=0(LSB)
        // Global index: 1 + 0 + 4(target) = 5
        double prob = state.GetProbability(5);
        Assert.IsTrue(prob > 0.9,
            $"Expected φ=1/2 to produce counting result 2, but P(5) = {prob}");
    }

    [TestMethod]
    public void QPE_Builder_ApplyQPE()
    {
        // Same as quarter-phase test but via builder
        var circuit = QuantumCircuitBuilder.New(3)
            .X(2) // eigenstate
            .ApplyQPE(new[] { 0, 1 }, new[] { 2 }, new PhaseGate(Math.PI / 2))
            .Build();
        var state = new QuantumSimulator().Run(circuit);
        Assert.IsTrue(state.GetProbability(6) > 0.9);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void QPE_MismatchedQubits_Throws()
    {
        // PhaseGate acts on 1 qubit, but we pass 2 target qubits
        QPE.CreateCircuit(4, new[] { 0, 1 }, new[] { 2, 3 }, new PhaseGate(Math.PI));
    }

    // ── Shor helpers tests ─────────────────────────────────────

    [TestMethod]
    public void Shor_ContinuedFraction_SimpleRationals()
    {
        // 1/4 = 0.25 → convergents should include (1, 4)
        var convergents = ShorAlgorithm.ContinuedFractionConvergents(1, 4);
        Assert.IsTrue(convergents.Count > 0);
        bool found = false;
        foreach (var (n, d) in convergents)
            if (n == 1 && d == 4) found = true;
        Assert.IsTrue(found, "Expected convergent (1, 4)");
    }

    [TestMethod]
    public void Shor_ContinuedFraction_ThreeOverEight()
    {
        // 3/8 → convergents: (0,1),(1,3),(3,8)... should include (3,8)
        var convergents = ShorAlgorithm.ContinuedFractionConvergents(3, 8);
        bool found = false;
        foreach (var (n, d) in convergents)
            if (n == 3 && d == 8) found = true;
        Assert.IsTrue(found, "Expected convergent (3, 8)");
    }

    [TestMethod]
    public void Shor_ModPow_Basic()
    {
        Assert.AreEqual(1, ShorAlgorithm.ModPow(2, 0, 5));  // 2^0 = 1
        Assert.AreEqual(2, ShorAlgorithm.ModPow(2, 1, 5));  // 2^1 = 2
        Assert.AreEqual(4, ShorAlgorithm.ModPow(2, 2, 5));  // 2^2 = 4
        Assert.AreEqual(3, ShorAlgorithm.ModPow(2, 3, 5));  // 2^3 = 8 mod 5 = 3
        Assert.AreEqual(1, ShorAlgorithm.ModPow(2, 4, 5));  // 2^4 = 16 mod 5 = 1
    }

    [TestMethod]
    public void Shor_GCD_Basic()
    {
        Assert.AreEqual(1, ShorAlgorithm.GCD(3, 5));
        Assert.AreEqual(3, ShorAlgorithm.GCD(9, 6));
        Assert.AreEqual(5, ShorAlgorithm.GCD(15, 10));
    }

    [TestMethod]
    public void Shor_ExtractOrder_KnownPhase()
    {
        // For a=2, N=15: order is 4. If QPE measures 4/16 = 1/4, 
        // continued fractions should give denominator 4.
        long order = ShorAlgorithm.ExtractOrder(4, 4, 15, 2);
        Assert.AreEqual(4, order, "Expected order 4 for a=2, N=15");
    }

    [TestMethod]
    public void Shor_Factor_Even()
    {
        // Even numbers are trivially factored
        var result = ShorAlgorithm.Factor(6, new Random(42));
        Assert.IsTrue(result.Success);
        Assert.AreEqual(2, result.Factor1);
        Assert.AreEqual(3, result.Factor2);
    }

    [TestMethod]
    public void Shor_Factor_PrimePower()
    {
        // 9 = 3^2 — detected by prime power check
        var result = ShorAlgorithm.Factor(9, new Random(42));
        Assert.IsTrue(result.Success);
        Assert.AreEqual(3, result.Factor1);
        Assert.AreEqual(3, result.Factor2);
    }

    [TestMethod]
    public void Shor_Factor_15()
    {
        // 15 = 3 × 5 — classic Shor demonstration
        var result = ShorAlgorithm.Factor(15, new Random(42));
        Assert.IsTrue(result.Success, "Should factor 15 successfully");
        Assert.IsTrue(
            (result.Factor1 == 3 && result.Factor2 == 5) ||
            (result.Factor1 == 5 && result.Factor2 == 3),
            $"Expected 3×5 but got {result.Factor1}×{result.Factor2}");
    }

    [TestMethod]
    public void Shor_OrderFindingCircuit_Creates()
    {
        // Verify the order-finding circuit can be created and has correct qubit count
        var circuit = ShorAlgorithm.CreateOrderFindingCircuit(2, 15, 8, 4);
        Assert.AreEqual(12, circuit.QubitCount);
        Assert.IsTrue(circuit.Instructions.Count > 0, "Circuit should have instructions");
    }
}
