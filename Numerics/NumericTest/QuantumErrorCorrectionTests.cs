using Microsoft.VisualStudio.TestTools.UnitTesting;
using CSharpNumerics.Physics.Quantum;
using CSharpNumerics.Physics.Quantum.ErrorCorrection;
using CSharpNumerics.Engines.Quantum;
using CSharpNumerics.Engines.Quantum.ErrorCorrection;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace NumericTest;

[TestClass]
public class QuantumErrorCorrectionTests
{
    private const double Tolerance = 1e-10;
    private const double LooseTolerance = 1e-6;

    // ── BitFlipCode3 — Code properties ──────────────────────────

    [TestMethod]
    public void BitFlipCode3_HasCorrectParameters()
    {
        var code = new BitFlipCode3();
        Assert.AreEqual(3, code.PhysicalQubits);
        Assert.AreEqual(1, code.LogicalQubits);
        Assert.AreEqual(1, code.Distance);
        Assert.AreEqual(2, code.SyndromeQubits);
    }

    [TestMethod]
    public void BitFlipCode3_HasTwoStabilizers()
    {
        var code = new BitFlipCode3();
        var stabs = code.GetStabilizers();
        Assert.AreEqual(2, stabs.Count);
        // Z₀Z₁
        Assert.AreEqual(2, stabs[0].Count);
        Assert.AreEqual((0, 'Z'), stabs[0][0]);
        Assert.AreEqual((1, 'Z'), stabs[0][1]);
        // Z₁Z₂
        Assert.AreEqual(2, stabs[1].Count);
        Assert.AreEqual((1, 'Z'), stabs[1][0]);
        Assert.AreEqual((2, 'Z'), stabs[1][1]);
    }

    [TestMethod]
    public void BitFlipCode3_CorrectionMap_HasFourEntries()
    {
        var code = new BitFlipCode3();
        var map = code.GetCorrectionMap();
        Assert.AreEqual(4, map.Count);
        // No error
        Assert.AreEqual(0, map[0b00].Count);
        // Error on qubit 0
        Assert.AreEqual(1, map[0b01].Count);
        Assert.AreEqual((0, 'X'), map[0b01][0]);
        // Error on qubit 2
        Assert.AreEqual(1, map[0b10].Count);
        Assert.AreEqual((2, 'X'), map[0b10][0]);
        // Error on qubit 1
        Assert.AreEqual(1, map[0b11].Count);
        Assert.AreEqual((1, 'X'), map[0b11][0]);
    }

    // ── PhaseFlipCode3 — Code properties ────────────────────────

    [TestMethod]
    public void PhaseFlipCode3_HasCorrectParameters()
    {
        var code = new PhaseFlipCode3();
        Assert.AreEqual(3, code.PhysicalQubits);
        Assert.AreEqual(1, code.LogicalQubits);
        Assert.AreEqual(1, code.Distance);
        Assert.AreEqual(2, code.SyndromeQubits);
    }

    [TestMethod]
    public void PhaseFlipCode3_HasXStabilizers()
    {
        var code = new PhaseFlipCode3();
        var stabs = code.GetStabilizers();
        Assert.AreEqual(2, stabs.Count);
        // X₀X₁
        Assert.AreEqual((0, 'X'), stabs[0][0]);
        Assert.AreEqual((1, 'X'), stabs[0][1]);
        // X₁X₂
        Assert.AreEqual((1, 'X'), stabs[1][0]);
        Assert.AreEqual((2, 'X'), stabs[1][1]);
    }

    [TestMethod]
    public void PhaseFlipCode3_CorrectionMap_UsesZGates()
    {
        var code = new PhaseFlipCode3();
        var map = code.GetCorrectionMap();
        Assert.AreEqual(4, map.Count);
        Assert.AreEqual(0, map[0b00].Count);
        Assert.AreEqual((0, 'Z'), map[0b01][0]);
        Assert.AreEqual((2, 'Z'), map[0b10][0]);
        Assert.AreEqual((1, 'Z'), map[0b11][0]);
    }

    // ── SyndromeDecoder ─────────────────────────────────────────

    [TestMethod]
    public void SyndromeDecoder_NoError_ReturnsEmptyList()
    {
        var decoder = new SyndromeDecoder(new BitFlipCode3());
        var corrections = decoder.Decode(0);
        Assert.AreEqual(0, corrections.Count);
    }

    [TestMethod]
    public void SyndromeDecoder_Syndrome01_CorrectsBitFlipOnQubit0()
    {
        var decoder = new SyndromeDecoder(new BitFlipCode3());
        var corrections = decoder.Decode(0b01);
        Assert.AreEqual(1, corrections.Count);
        Assert.AreEqual((0, 'X'), corrections[0]);
    }

    [TestMethod]
    public void SyndromeDecoder_Syndrome10_CorrectsBitFlipOnQubit2()
    {
        var decoder = new SyndromeDecoder(new BitFlipCode3());
        var corrections = decoder.Decode(0b10);
        Assert.AreEqual(1, corrections.Count);
        Assert.AreEqual((2, 'X'), corrections[0]);
    }

    [TestMethod]
    public void SyndromeDecoder_Syndrome11_CorrectsBitFlipOnQubit1()
    {
        var decoder = new SyndromeDecoder(new BitFlipCode3());
        var corrections = decoder.Decode(0b11);
        Assert.AreEqual(1, corrections.Count);
        Assert.AreEqual((1, 'X'), corrections[0]);
    }

    [TestMethod]
    public void SyndromeDecoder_FromBitsArray_Works()
    {
        var decoder = new SyndromeDecoder(new BitFlipCode3());
        var corrections = decoder.Decode(new[] { 1, 1 }); // syndrome = 0b11
        Assert.AreEqual(1, corrections.Count);
        Assert.AreEqual((1, 'X'), corrections[0]);
    }

    [TestMethod]
    public void SyndromeDecoder_UnknownSyndrome_ReturnsEmpty()
    {
        // Custom map with only one entry
        var map = new Dictionary<int, List<(int qubit, char pauli)>>
        {
            { 0, new List<(int, char)>() }
        };
        var decoder = new SyndromeDecoder(map);
        var corrections = decoder.Decode(99);
        Assert.AreEqual(0, corrections.Count);
    }

    // ── ErrorCorrectionSimulator — BitFlip ──────────────────────

    [TestMethod]
    public void BitFlip_NoError_RecoversWith100PercentFidelity()
    {
        var code = new BitFlipCode3();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        // |0⟩ state
        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(1, 0);
        initial[1] = new ComplexNumber(0, 0);

        var result = sim.RunBitFlipCorrection(code, initial, null, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Fidelity should be ~1.0 but was {result.Fidelity}");
        Assert.AreEqual(0, result.Syndrome);
    }

    [TestMethod]
    public void BitFlip_NoError_RecoversSuperpositionState()
    {
        var code = new BitFlipCode3();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        // |+⟩ = (|0⟩ + |1⟩)/√2
        double inv = 1.0 / Math.Sqrt(2.0);
        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(inv, 0);
        initial[1] = new ComplexNumber(inv, 0);

        var result = sim.RunBitFlipCorrection(code, initial, null, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Fidelity should be ~1.0 but was {result.Fidelity}");
    }

    [TestMethod]
    public void BitFlip_SingleErrorOnQubit0_Corrected()
    {
        var code = new BitFlipCode3();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(1, 0);
        initial[1] = new ComplexNumber(0, 0);

        var errors = new List<(QuantumGate gate, int qubit)>
        {
            (new PauliXGate(), 0)
        };

        var result = sim.RunBitFlipCorrection(code, initial, errors, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Fidelity after correction should be ~1.0 but was {result.Fidelity}");
        Assert.AreEqual(0b01, result.Syndrome, "Syndrome should indicate error on qubit 0");
    }

    [TestMethod]
    public void BitFlip_SingleErrorOnQubit1_Corrected()
    {
        var code = new BitFlipCode3();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(1, 0);
        initial[1] = new ComplexNumber(0, 0);

        var errors = new List<(QuantumGate gate, int qubit)>
        {
            (new PauliXGate(), 1)
        };

        var result = sim.RunBitFlipCorrection(code, initial, errors, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Fidelity should be ~1.0 but was {result.Fidelity}");
        Assert.AreEqual(0b11, result.Syndrome, "Syndrome should indicate error on qubit 1");
    }

    [TestMethod]
    public void BitFlip_SingleErrorOnQubit2_Corrected()
    {
        var code = new BitFlipCode3();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(1, 0);
        initial[1] = new ComplexNumber(0, 0);

        var errors = new List<(QuantumGate gate, int qubit)>
        {
            (new PauliXGate(), 2)
        };

        var result = sim.RunBitFlipCorrection(code, initial, errors, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Fidelity should be ~1.0 but was {result.Fidelity}");
        Assert.AreEqual(0b10, result.Syndrome, "Syndrome should indicate error on qubit 2");
    }

    [TestMethod]
    public void BitFlip_SuperpositionState_SingleError_Corrected()
    {
        var code = new BitFlipCode3();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        // |ψ⟩ = (3|0⟩ + 4i|1⟩)/5
        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(0.6, 0);
        initial[1] = new ComplexNumber(0, 0.8);

        var errors = new List<(QuantumGate gate, int qubit)>
        {
            (new PauliXGate(), 1)
        };

        var result = sim.RunBitFlipCorrection(code, initial, errors, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Fidelity should be ~1.0 but was {result.Fidelity}");
    }

    // ── ErrorCorrectionSimulator — PhaseFlip ────────────────────

    [TestMethod]
    public void PhaseFlip_NoError_RecoversWith100PercentFidelity()
    {
        var code = new PhaseFlipCode3();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(1, 0);
        initial[1] = new ComplexNumber(0, 0);

        var result = sim.RunPhaseFlipCorrection(code, initial, null, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Fidelity should be ~1.0 but was {result.Fidelity}");
        Assert.AreEqual(0, result.Syndrome);
    }

    [TestMethod]
    public void PhaseFlip_SingleZErrorOnQubit0_Corrected()
    {
        var code = new PhaseFlipCode3();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(1, 0);
        initial[1] = new ComplexNumber(0, 0);

        var errors = new List<(QuantumGate gate, int qubit)>
        {
            (new PauliZGate(), 0)
        };

        var result = sim.RunPhaseFlipCorrection(code, initial, errors, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Fidelity after correction should be ~1.0 but was {result.Fidelity}");
    }

    [TestMethod]
    public void PhaseFlip_SingleZErrorOnQubit1_Corrected()
    {
        var code = new PhaseFlipCode3();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(1, 0);
        initial[1] = new ComplexNumber(0, 0);

        var errors = new List<(QuantumGate gate, int qubit)>
        {
            (new PauliZGate(), 1)
        };

        var result = sim.RunPhaseFlipCorrection(code, initial, errors, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Fidelity should be ~1.0 but was {result.Fidelity}");
    }

    [TestMethod]
    public void PhaseFlip_SingleZErrorOnQubit2_Corrected()
    {
        var code = new PhaseFlipCode3();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(1, 0);
        initial[1] = new ComplexNumber(0, 0);

        var errors = new List<(QuantumGate gate, int qubit)>
        {
            (new PauliZGate(), 2)
        };

        var result = sim.RunPhaseFlipCorrection(code, initial, errors, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Fidelity should be ~1.0 but was {result.Fidelity}");
    }

    [TestMethod]
    public void PhaseFlip_SuperpositionState_SingleZError_Corrected()
    {
        var code = new PhaseFlipCode3();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        // |ψ⟩ = (|0⟩ + |1⟩)/√2
        double inv = 1.0 / Math.Sqrt(2.0);
        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(inv, 0);
        initial[1] = new ComplexNumber(inv, 0);

        var errors = new List<(QuantumGate gate, int qubit)>
        {
            (new PauliZGate(), 1)
        };

        var result = sim.RunPhaseFlipCorrection(code, initial, errors, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Fidelity should be ~1.0 but was {result.Fidelity}");
    }

    // ── Monte Carlo comparison ──────────────────────────────────

    [TestMethod]
    public void MonteCarlo_ProtectedBeatUnprotected_AtModerateErrorRate()
    {
        var code = new BitFlipCode3();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        // |+⟩ state
        double inv = 1.0 / Math.Sqrt(2.0);
        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(inv, 0);
        initial[1] = new ComplexNumber(inv, 0);

        // Error rate 10% per qubit — QEC should help
        var (protectedFid, unprotectedFid) = sim.RunMonteCarloComparison(
            code, initial, errorRate: 0.10, rounds: 500, random: rng);

        Assert.IsTrue(protectedFid > unprotectedFid,
            $"Protected fidelity ({protectedFid:F4}) should exceed unprotected ({unprotectedFid:F4})");
    }

    // ── Logical operators ───────────────────────────────────────

    [TestMethod]
    public void BitFlipCode3_LogicalX_IsXXX()
    {
        var code = new BitFlipCode3();
        var lx = code.GetLogicalX();
        Assert.AreEqual(3, lx.Count);
        Assert.AreEqual((0, 'X'), lx[0]);
        Assert.AreEqual((1, 'X'), lx[1]);
        Assert.AreEqual((2, 'X'), lx[2]);
    }

    [TestMethod]
    public void BitFlipCode3_LogicalZ_IsZ0()
    {
        var code = new BitFlipCode3();
        var lz = code.GetLogicalZ();
        Assert.AreEqual(1, lz.Count);
        Assert.AreEqual((0, 'Z'), lz[0]);
    }

    [TestMethod]
    public void PhaseFlipCode3_LogicalX_IsX0()
    {
        var code = new PhaseFlipCode3();
        var lx = code.GetLogicalX();
        Assert.AreEqual(1, lx.Count);
        Assert.AreEqual((0, 'X'), lx[0]);
    }

    [TestMethod]
    public void PhaseFlipCode3_LogicalZ_IsZZZ()
    {
        var code = new PhaseFlipCode3();
        var lz = code.GetLogicalZ();
        Assert.AreEqual(3, lz.Count);
        Assert.AreEqual((0, 'Z'), lz[0]);
        Assert.AreEqual((1, 'Z'), lz[1]);
        Assert.AreEqual((2, 'Z'), lz[2]);
    }

    // ── Edge cases ──────────────────────────────────────────────

    [TestMethod]
    public void BitFlip_State1_NoError_RecoversPerfectly()
    {
        var code = new BitFlipCode3();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        // |1⟩ state
        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(0, 0);
        initial[1] = new ComplexNumber(1, 0);

        var result = sim.RunBitFlipCorrection(code, initial, null, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Fidelity should be ~1.0 but was {result.Fidelity}");
    }

    [TestMethod]
    public void BitFlip_State1_SingleError_Corrected()
    {
        var code = new BitFlipCode3();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        // |1⟩ state
        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(0, 0);
        initial[1] = new ComplexNumber(1, 0);

        var errors = new List<(QuantumGate gate, int qubit)>
        {
            (new PauliXGate(), 2)
        };

        var result = sim.RunBitFlipCorrection(code, initial, errors, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Fidelity should be ~1.0 but was {result.Fidelity}");
    }

    [TestMethod]
    public void PhaseFlip_State1_NoError_RecoversPerfectly()
    {
        var code = new PhaseFlipCode3();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        // |1⟩ state
        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(0, 0);
        initial[1] = new ComplexNumber(1, 0);

        var result = sim.RunPhaseFlipCorrection(code, initial, null, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Fidelity should be ~1.0 but was {result.Fidelity}");
    }

    // ══════════════════════════════════════════════════════════════
    //  ShorCode9 [[9,1,3]] — Phase 2
    // ══════════════════════════════════════════════════════════════

    // ── ShorCode9 — Code properties ─────────────────────────────

    [TestMethod]
    public void ShorCode9_HasCorrectParameters()
    {
        var code = new ShorCode9();
        Assert.AreEqual(9, code.PhysicalQubits);
        Assert.AreEqual(1, code.LogicalQubits);
        Assert.AreEqual(3, code.Distance);
        Assert.AreEqual(8, code.SyndromeQubits);
    }

    [TestMethod]
    public void ShorCode9_Has8Stabilizers()
    {
        var code = new ShorCode9();
        var stabs = code.GetStabilizers();
        Assert.AreEqual(8, stabs.Count);
        // First 6 are ZZ pairs
        Assert.AreEqual((0, 'Z'), stabs[0][0]);
        Assert.AreEqual((1, 'Z'), stabs[0][1]);
        Assert.AreEqual((7, 'Z'), stabs[5][0]);
        Assert.AreEqual((8, 'Z'), stabs[5][1]);
        // Last 2 are 6-qubit X strings
        Assert.AreEqual(6, stabs[6].Count);
        Assert.AreEqual(6, stabs[7].Count);
    }

    [TestMethod]
    public void ShorCode9_CorrectionMap_Has256Entries()
    {
        var code = new ShorCode9();
        var map = code.GetCorrectionMap();
        Assert.AreEqual(256, map.Count);
        // No error
        Assert.AreEqual(0, map[0].Count);
    }

    [TestMethod]
    public void ShorCode9_LogicalX_IsXOnAll9()
    {
        var code = new ShorCode9();
        var lx = code.GetLogicalX();
        Assert.AreEqual(9, lx.Count);
        for (int i = 0; i < 9; i++)
            Assert.AreEqual((i, 'X'), lx[i]);
    }

    [TestMethod]
    public void ShorCode9_LogicalZ_IsZOnBlockLeaders()
    {
        var code = new ShorCode9();
        var lz = code.GetLogicalZ();
        Assert.AreEqual(3, lz.Count);
        Assert.AreEqual((0, 'Z'), lz[0]);
        Assert.AreEqual((3, 'Z'), lz[1]);
        Assert.AreEqual((6, 'Z'), lz[2]);
    }

    // ── Shor code — No error ────────────────────────────────────

    [TestMethod]
    public void Shor9_NoError_State0_RecoversPerfectly()
    {
        var code = new ShorCode9();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(1, 0);
        initial[1] = new ComplexNumber(0, 0);

        var result = sim.RunShorCorrection(code, initial, null, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Shor |0⟩ no error: fidelity should be ~1.0 but was {result.Fidelity}");
        Assert.AreEqual(0, result.Syndrome);
    }

    [TestMethod]
    public void Shor9_NoError_State1_RecoversPerfectly()
    {
        var code = new ShorCode9();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(0, 0);
        initial[1] = new ComplexNumber(1, 0);

        var result = sim.RunShorCorrection(code, initial, null, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Shor |1⟩ no error: fidelity should be ~1.0 but was {result.Fidelity}");
    }

    [TestMethod]
    public void Shor9_NoError_Superposition_RecoversPerfectly()
    {
        var code = new ShorCode9();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        double inv = 1.0 / Math.Sqrt(2.0);
        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(inv, 0);
        initial[1] = new ComplexNumber(inv, 0);

        var result = sim.RunShorCorrection(code, initial, null, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Shor |+⟩ no error: fidelity should be ~1.0 but was {result.Fidelity}");
    }

    // ── Shor code — Single X (bit-flip) errors ──────────────────

    [TestMethod]
    public void Shor9_BitFlipOnQubit0_Corrected()
    {
        var code = new ShorCode9();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(0.6, 0);
        initial[1] = new ComplexNumber(0, 0.8);

        var errors = new List<(QuantumGate gate, int qubit)> { (new PauliXGate(), 0) };
        var result = sim.RunShorCorrection(code, initial, errors, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Shor X₀ corrected: fidelity={result.Fidelity}");
    }

    [TestMethod]
    public void Shor9_BitFlipOnQubit4_Corrected()
    {
        var code = new ShorCode9();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(0.6, 0);
        initial[1] = new ComplexNumber(0, 0.8);

        var errors = new List<(QuantumGate gate, int qubit)> { (new PauliXGate(), 4) };
        var result = sim.RunShorCorrection(code, initial, errors, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Shor X₄ corrected: fidelity={result.Fidelity}");
    }

    [TestMethod]
    public void Shor9_BitFlipOnQubit8_Corrected()
    {
        var code = new ShorCode9();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(0.6, 0);
        initial[1] = new ComplexNumber(0, 0.8);

        var errors = new List<(QuantumGate gate, int qubit)> { (new PauliXGate(), 8) };
        var result = sim.RunShorCorrection(code, initial, errors, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Shor X₈ corrected: fidelity={result.Fidelity}");
    }

    // ── Shor code — Single Z (phase-flip) errors ────────────────

    [TestMethod]
    public void Shor9_PhaseFlipOnQubit0_Corrected()
    {
        var code = new ShorCode9();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(0.6, 0);
        initial[1] = new ComplexNumber(0, 0.8);

        var errors = new List<(QuantumGate gate, int qubit)> { (new PauliZGate(), 0) };
        var result = sim.RunShorCorrection(code, initial, errors, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Shor Z₀ corrected: fidelity={result.Fidelity}");
    }

    [TestMethod]
    public void Shor9_PhaseFlipOnQubit5_Corrected()
    {
        var code = new ShorCode9();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(0.6, 0);
        initial[1] = new ComplexNumber(0, 0.8);

        var errors = new List<(QuantumGate gate, int qubit)> { (new PauliZGate(), 5) };
        var result = sim.RunShorCorrection(code, initial, errors, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Shor Z₅ corrected: fidelity={result.Fidelity}");
    }

    [TestMethod]
    public void Shor9_PhaseFlipOnQubit7_Corrected()
    {
        var code = new ShorCode9();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(0.6, 0);
        initial[1] = new ComplexNumber(0, 0.8);

        var errors = new List<(QuantumGate gate, int qubit)> { (new PauliZGate(), 7) };
        var result = sim.RunShorCorrection(code, initial, errors, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Shor Z₇ corrected: fidelity={result.Fidelity}");
    }

    // ── Shor code — Single Y (combined XZ) errors ───────────────

    [TestMethod]
    public void Shor9_YErrorOnQubit0_Corrected()
    {
        var code = new ShorCode9();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(0.6, 0);
        initial[1] = new ComplexNumber(0, 0.8);

        var errors = new List<(QuantumGate gate, int qubit)> { (new PauliYGate(), 0) };
        var result = sim.RunShorCorrection(code, initial, errors, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Shor Y₀ corrected: fidelity={result.Fidelity}");
    }

    [TestMethod]
    public void Shor9_YErrorOnQubit4_Corrected()
    {
        var code = new ShorCode9();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(0.6, 0);
        initial[1] = new ComplexNumber(0, 0.8);

        var errors = new List<(QuantumGate gate, int qubit)> { (new PauliYGate(), 4) };
        var result = sim.RunShorCorrection(code, initial, errors, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Shor Y₄ corrected: fidelity={result.Fidelity}");
    }

    // ══════════════════════════════════════════════════════════════
    //  SteaneCode7 [[7,1,3]] — Phase 3
    // ══════════════════════════════════════════════════════════════

    // ── SteaneCode7 — Code properties ───────────────────────────

    [TestMethod]
    public void SteaneCode7_HasCorrectParameters()
    {
        var code = new SteaneCode7();
        Assert.AreEqual(7, code.PhysicalQubits);
        Assert.AreEqual(1, code.LogicalQubits);
        Assert.AreEqual(3, code.Distance);
        Assert.AreEqual(6, code.SyndromeQubits);
    }

    [TestMethod]
    public void SteaneCode7_Has6Stabilizers()
    {
        var code = new SteaneCode7();
        var stabs = code.GetStabilizers();
        Assert.AreEqual(6, stabs.Count);
        // First 3 are Z-type (4 qubits each)
        for (int i = 0; i < 3; i++)
        {
            Assert.AreEqual(4, stabs[i].Count);
            foreach (var (_, p) in stabs[i])
                Assert.AreEqual('Z', p);
        }
        // Last 3 are X-type (4 qubits each)
        for (int i = 3; i < 6; i++)
        {
            Assert.AreEqual(4, stabs[i].Count);
            foreach (var (_, p) in stabs[i])
                Assert.AreEqual('X', p);
        }
    }

    [TestMethod]
    public void SteaneCode7_CorrectionMap_Has64Entries()
    {
        var code = new SteaneCode7();
        var map = code.GetCorrectionMap();
        Assert.AreEqual(64, map.Count);
        // No error
        Assert.AreEqual(0, map[0].Count);
    }

    [TestMethod]
    public void SteaneCode7_CorrectionMap_SingleXErrors()
    {
        var code = new SteaneCode7();
        var map = code.GetCorrectionMap();
        // Z-syndrome for X_j: {1→q3, 2→q1, 3→q5, 4→q0, 5→q4, 6→q2, 7→q6}
        int[] expectedQubit = { -1, 3, 1, 5, 0, 4, 2, 6 };
        for (int s = 1; s <= 7; s++)
        {
            Assert.AreEqual(1, map[s].Count, $"syndrome {s} should have one correction");
            Assert.AreEqual(expectedQubit[s], map[s][0].qubit);
            Assert.AreEqual('X', map[s][0].pauli);
        }
    }

    [TestMethod]
    public void SteaneCode7_CorrectionMap_SingleZErrors()
    {
        var code = new SteaneCode7();
        var map = code.GetCorrectionMap();
        // X-syndrome for Z_j in bits 3-5: same mapping shifted
        int[] expectedQubit = { -1, 3, 1, 5, 0, 4, 2, 6 };
        for (int s = 1; s <= 7; s++)
        {
            int fullSyndrome = s << 3; // bits 3-5
            Assert.AreEqual(1, map[fullSyndrome].Count, $"syndrome {fullSyndrome} should have one correction");
            Assert.AreEqual(expectedQubit[s], map[fullSyndrome][0].qubit);
            Assert.AreEqual('Z', map[fullSyndrome][0].pauli);
        }
    }

    [TestMethod]
    public void SteaneCode7_LogicalX_IsXOnAll7()
    {
        var code = new SteaneCode7();
        var lx = code.GetLogicalX();
        Assert.AreEqual(7, lx.Count);
        for (int i = 0; i < 7; i++)
            Assert.AreEqual((i, 'X'), lx[i]);
    }

    [TestMethod]
    public void SteaneCode7_LogicalZ_IsZOnAll7()
    {
        var code = new SteaneCode7();
        var lz = code.GetLogicalZ();
        Assert.AreEqual(7, lz.Count);
        for (int i = 0; i < 7; i++)
            Assert.AreEqual((i, 'Z'), lz[i]);
    }

    // ── Steane code — No error ──────────────────────────────────

    [TestMethod]
    public void Steane7_NoError_State0_RecoversPerfectly()
    {
        var code = new SteaneCode7();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(1, 0);
        initial[1] = new ComplexNumber(0, 0);

        var result = sim.RunSteaneCorrection(code, initial, null, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Steane |0⟩ no error: fidelity should be ~1.0 but was {result.Fidelity}");
        Assert.AreEqual(0, result.Syndrome);
    }

    [TestMethod]
    public void Steane7_NoError_State1_RecoversPerfectly()
    {
        var code = new SteaneCode7();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(0, 0);
        initial[1] = new ComplexNumber(1, 0);

        var result = sim.RunSteaneCorrection(code, initial, null, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Steane |1⟩ no error: fidelity should be ~1.0 but was {result.Fidelity}");
    }

    [TestMethod]
    public void Steane7_NoError_Superposition_RecoversPerfectly()
    {
        var code = new SteaneCode7();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        double inv = 1.0 / Math.Sqrt(2.0);
        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(inv, 0);
        initial[1] = new ComplexNumber(inv, 0);

        var result = sim.RunSteaneCorrection(code, initial, null, rng);

        Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
            $"Steane |+⟩ no error: fidelity should be ~1.0 but was {result.Fidelity}");
    }

    // ── Steane code — Single X (bit-flip) errors ────────────────

    [TestMethod]
    public void Steane7_BitFlipOnEachQubit_Corrected()
    {
        var code = new SteaneCode7();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(0.6, 0);
        initial[1] = new ComplexNumber(0, 0.8);

        for (int q = 0; q < 7; q++)
        {
            var errors = new List<(QuantumGate gate, int qubit)> { (new PauliXGate(), q) };
            var result = sim.RunSteaneCorrection(code, initial, errors, new Random(42));

            Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
                $"Steane X_{q} corrected: fidelity={result.Fidelity}");
            Assert.IsTrue(result.Syndrome > 0, $"X_{q} should produce non-zero syndrome");
        }
    }

    // ── Steane code — Single Z (phase-flip) errors ──────────────

    [TestMethod]
    public void Steane7_PhaseFlipOnEachQubit_Corrected()
    {
        var code = new SteaneCode7();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(0.6, 0);
        initial[1] = new ComplexNumber(0, 0.8);

        for (int q = 0; q < 7; q++)
        {
            var errors = new List<(QuantumGate gate, int qubit)> { (new PauliZGate(), q) };
            var result = sim.RunSteaneCorrection(code, initial, errors, new Random(42));

            Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
                $"Steane Z_{q} corrected: fidelity={result.Fidelity}");
            Assert.IsTrue(result.Syndrome > 0, $"Z_{q} should produce non-zero syndrome");
        }
    }

    // ── Steane code — Single Y (combined XZ) errors ─────────────

    [TestMethod]
    public void Steane7_YErrorOnEachQubit_Corrected()
    {
        var code = new SteaneCode7();
        var sim = new ErrorCorrectionSimulator();
        var rng = new Random(42);

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(0.6, 0);
        initial[1] = new ComplexNumber(0, 0.8);

        for (int q = 0; q < 7; q++)
        {
            var errors = new List<(QuantumGate gate, int qubit)> { (new PauliYGate(), q) };
            var result = sim.RunSteaneCorrection(code, initial, errors, new Random(42));

            Assert.IsTrue(result.Fidelity > 1.0 - Tolerance,
                $"Steane Y_{q} corrected: fidelity={result.Fidelity}");
            Assert.IsTrue(result.Syndrome > 0, $"Y_{q} should produce non-zero syndrome");
        }
    }

    // ── Steane code — Syndrome value verification ───────────────

    [TestMethod]
    public void Steane7_XError_ProducesCorrectZSyndrome()
    {
        // Z-syndrome for X_j follows the Hamming decoding:
        // qubit 0→4, 1→2, 2→6, 3→1, 4→5, 5→3, 6→7
        int[] expectedZSynd = { 4, 2, 6, 1, 5, 3, 7 };

        var code = new SteaneCode7();
        var sim = new ErrorCorrectionSimulator();

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(1, 0);
        initial[1] = new ComplexNumber(0, 0);

        for (int q = 0; q < 7; q++)
        {
            var errors = new List<(QuantumGate gate, int qubit)> { (new PauliXGate(), q) };
            var result = sim.RunSteaneCorrection(code, initial, errors, new Random(42));

            int zPart = result.Syndrome & 7;
            Assert.AreEqual(expectedZSynd[q], zPart,
                $"X_{q}: expected Z-syndrome={expectedZSynd[q]}, got {zPart}");
        }
    }

    [TestMethod]
    public void Steane7_ZError_ProducesCorrectXSyndrome()
    {
        // X-syndrome for Z_j (bits 3-5) follows the same Hamming map
        int[] expectedXSynd = { 4, 2, 6, 1, 5, 3, 7 };

        var code = new SteaneCode7();
        var sim = new ErrorCorrectionSimulator();

        var initial = new ComplexVectorN(2);
        initial[0] = new ComplexNumber(1, 0);
        initial[1] = new ComplexNumber(0, 0);

        for (int q = 0; q < 7; q++)
        {
            var errors = new List<(QuantumGate gate, int qubit)> { (new PauliZGate(), q) };
            var result = sim.RunSteaneCorrection(code, initial, errors, new Random(42));

            int xPart = (result.Syndrome >> 3) & 7;
            Assert.AreEqual(expectedXSynd[q], xPart,
                $"Z_{q}: expected X-syndrome={expectedXSynd[q]}, got {xPart}");
        }
    }
}
