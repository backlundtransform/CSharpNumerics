using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Quantum;
using CSharpNumerics.Physics.Quantum.ErrorCorrection;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Quantum.ErrorCorrection;

/// <summary>
/// Simulates quantum error correction: encode → noise → syndrome extraction →
/// classical decode → correction → decode. Reports fidelity of the recovered
/// state against the ideal (unencoded, noiseless) state.
/// </summary>
public class ErrorCorrectionSimulator
{
    private readonly QuantumSimulator _simulator = new QuantumSimulator();

    /// <summary>
    /// Result of a single error-correction simulation run.
    /// </summary>
    public class Result
    {
        /// <summary>Fidelity of the recovered state vs the ideal state.</summary>
        public double Fidelity { get; set; }

        /// <summary>The syndrome measured during the run.</summary>
        public int Syndrome { get; set; }

        /// <summary>The corrective operations applied (empty if no error detected).</summary>
        public List<(int qubit, char pauli)> Corrections { get; set; }

        /// <summary>The recovered quantum state (data qubits only).</summary>
        public QuantumState RecoveredState { get; set; }
    }

    /// <summary>
    /// Runs a full error-correction cycle for the bit-flip code:
    ///   1. Prepare the logical state by encoding <paramref name="initialState"/> using CNOT gates
    ///   2. Inject an error (if specified) on data qubits
    ///   3. Extract syndrome using ancilla qubits
    ///   4. Decode syndrome and apply corrective gates
    ///   5. Decode (undo encoding) and return fidelity vs. original state
    /// </summary>
    /// <param name="code">The QEC code to use.</param>
    /// <param name="initialState">Initial 1-qubit state to protect (2 amplitudes).</param>
    /// <param name="errorGates">
    /// Optional list of (gate, qubit index) errors to inject on data qubits before syndrome extraction.
    /// Qubit indices are relative to the data qubits (0-based).
    /// </param>
    /// <param name="random">Random number generator for measurements.</param>
    public Result RunBitFlipCorrection(
        BitFlipCode3 code,
        ComplexVectorN initialState,
        List<(QuantumGate gate, int qubit)> errorGates,
        Random random)
    {
        // Total qubits: 3 data + 2 ancilla = 5
        int totalQubits = code.PhysicalQubits + code.SyndromeQubits;
        var circuit = new QuantumCircuit(totalQubits);

        // --- Prepare initial state on qubit 0 ---
        // We'll build the state manually by running a prep circuit
        // Actually, we set up the state vector directly after encoding

        // First: build encoding circuit on a clean 5-qubit register
        // We need to set the initial state on qubit 0 first, then encode.
        // We'll do this by building the full state vector.

        // Step 1: Create the initial 5-qubit state with qubit 0 in |ψ⟩
        int stateSize = 1 << totalQubits;
        var amplitudes = new ComplexVectorN(stateSize);
        // |ψ⟩ ⊗ |0000⟩: qubit 0 has the logical state, rest are |0⟩
        // In our bit ordering, qubit 0 is the LSB.
        // So amplitude of |00000⟩ = α, amplitude of |00001⟩ = β (qubit 0 = 1)
        amplitudes[0] = initialState[0]; // α coefficient for qubit 0 = |0⟩
        amplitudes[1] = initialState[1]; // β coefficient for qubit 0 = |1⟩

        // Step 2: Encode — CNOT(0,1), CNOT(0,2) to get α|000⟩+β|111⟩ on data qubits
        var cnot = new CNOTGate();
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 1 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 2 }, totalQubits);

        // Step 3: Inject errors (if any) on data qubits
        if (errorGates != null)
        {
            foreach (var (gate, qubit) in errorGates)
            {
                amplitudes = gate.Apply(amplitudes, new[] { qubit }, totalQubits);
            }
        }

        // Step 4: Syndrome extraction using ancilla qubits 3 and 4
        // Stabilizer Z₀Z₁ → ancilla 3: CNOT(0,3), CNOT(1,3)
        // Stabilizer Z₁Z₂ → ancilla 4: CNOT(1,4), CNOT(2,4)
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 3 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 1, 3 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 1, 4 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 2, 4 }, totalQubits);

        // Step 5: Measure ancilla qubits to get syndrome
        var state = new QuantumState(amplitudes);
        int s0 = state.MeasureQubit(3, random); // stabilizer Z₀Z₁
        int s1 = state.MeasureQubit(4, random); // stabilizer Z₁Z₂
        int syndrome = s0 | (s1 << 1);

        // Step 6: Decode syndrome and apply correction to data qubits
        var decoder = new SyndromeDecoder(code);
        var corrections = decoder.Decode(syndrome);

        amplitudes = state.Amplitudes;
        foreach (var (qubit, pauli) in corrections)
        {
            var corrGate = GetPauliGate(pauli);
            amplitudes = corrGate.Apply(amplitudes, new[] { qubit }, totalQubits);
        }

        // Step 7: Decode (undo encoding) — reverse of CNOT(0,2), CNOT(0,1)
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 2 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 1 }, totalQubits);

        // Step 8: Extract logical qubit (qubit 0) state
        // After decoding, qubits 1,2 should be |0⟩ and 3,4 are the measured syndrome.
        // Extract amplitudes where other qubits match expected values.
        int ancillaBits = (s0 << 3) | (s1 << 4);
        var recoveredAmplitudes = ExtractLogicalQubit(amplitudes, 0, ancillaBits, totalQubits);
        var recoveredState = new QuantumState(recoveredAmplitudes);

        // Step 9: Compute fidelity vs ideal state
        double fidelity = QuantumFidelity.Fidelity(initialState, recoveredAmplitudes);

        return new Result
        {
            Fidelity = fidelity,
            Syndrome = syndrome,
            Corrections = corrections,
            RecoveredState = recoveredState
        };
    }

    /// <summary>
    /// Runs a full error-correction cycle for the phase-flip code:
    ///   1. Prepare + encode using CNOT + Hadamard
    ///   2. Inject error
    ///   3. Syndrome extraction in X-basis (H → CNOT → H on ancillas)
    ///   4. Correct (Z gate)
    ///   5. Decode and return fidelity
    /// </summary>
    public Result RunPhaseFlipCorrection(
        PhaseFlipCode3 code,
        ComplexVectorN initialState,
        List<(QuantumGate gate, int qubit)> errorGates,
        Random random)
    {
        int totalQubits = code.PhysicalQubits + code.SyndromeQubits;
        int stateSize = 1 << totalQubits;
        var amplitudes = new ComplexVectorN(stateSize);
        amplitudes[0] = initialState[0];
        amplitudes[1] = initialState[1];

        var cnot = new CNOTGate();
        var hadamard = new HadamardGate();

        // Encode: CNOT(0,1), CNOT(0,2) → α|000⟩+β|111⟩
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 1 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 2 }, totalQubits);

        // Then Hadamard each data qubit → α|+++⟩+β|---⟩
        amplitudes = hadamard.Apply(amplitudes, new[] { 0 }, totalQubits);
        amplitudes = hadamard.Apply(amplitudes, new[] { 1 }, totalQubits);
        amplitudes = hadamard.Apply(amplitudes, new[] { 2 }, totalQubits);

        // Inject errors
        if (errorGates != null)
        {
            foreach (var (gate, qubit) in errorGates)
            {
                amplitudes = gate.Apply(amplitudes, new[] { qubit }, totalQubits);
            }
        }

        // Syndrome extraction for X-stabilizers:
        // Hadamard data qubits to convert X-basis to Z-basis
        amplitudes = hadamard.Apply(amplitudes, new[] { 0 }, totalQubits);
        amplitudes = hadamard.Apply(amplitudes, new[] { 1 }, totalQubits);
        amplitudes = hadamard.Apply(amplitudes, new[] { 2 }, totalQubits);

        // Z checks (same as bit-flip code)
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 3 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 1, 3 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 1, 4 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 2, 4 }, totalQubits);

        // Hadamard back
        amplitudes = hadamard.Apply(amplitudes, new[] { 0 }, totalQubits);
        amplitudes = hadamard.Apply(amplitudes, new[] { 1 }, totalQubits);
        amplitudes = hadamard.Apply(amplitudes, new[] { 2 }, totalQubits);

        // Measure ancillas
        var state = new QuantumState(amplitudes);
        int s0 = state.MeasureQubit(3, random);
        int s1 = state.MeasureQubit(4, random);
        int syndrome = s0 | (s1 << 1);

        // Decode syndrome and apply Z corrections
        var decoder = new SyndromeDecoder(code);
        var corrections = decoder.Decode(syndrome);

        amplitudes = state.Amplitudes;
        foreach (var (qubit, pauli) in corrections)
        {
            var corrGate = GetPauliGate(pauli);
            amplitudes = corrGate.Apply(amplitudes, new[] { qubit }, totalQubits);
        }

        // Decode: H on each data qubit, then reverse CNOTs
        amplitudes = hadamard.Apply(amplitudes, new[] { 0 }, totalQubits);
        amplitudes = hadamard.Apply(amplitudes, new[] { 1 }, totalQubits);
        amplitudes = hadamard.Apply(amplitudes, new[] { 2 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 2 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 1 }, totalQubits);

        var recoveredAmplitudes = ExtractLogicalQubit(amplitudes, 0, (s0 << 3) | (s1 << 4), totalQubits);
        var recoveredState = new QuantumState(recoveredAmplitudes);
        double fidelity = QuantumFidelity.Fidelity(initialState, recoveredAmplitudes);

        return new Result
        {
            Fidelity = fidelity,
            Syndrome = syndrome,
            Corrections = corrections,
            RecoveredState = recoveredState
        };
    }

    /// <summary>
    /// Runs many rounds with random single-qubit errors at given error rate and
    /// returns the average fidelity. Compares protected (with QEC) vs unprotected scenarios.
    /// </summary>
    /// <param name="code">The bit-flip QEC code.</param>
    /// <param name="initialState">The logical state to protect.</param>
    /// <param name="errorRate">Probability of a bit-flip error per qubit per round.</param>
    /// <param name="rounds">Number of Monte Carlo rounds.</param>
    /// <param name="random">Random number generator.</param>
    /// <returns>(averageFidelityProtected, averageFidelityUnprotected)</returns>
    public (double protectedFidelity, double unprotectedFidelity) RunMonteCarloComparison(
        BitFlipCode3 code,
        ComplexVectorN initialState,
        double errorRate,
        int rounds,
        Random random)
    {
        double sumProtected = 0;
        double sumUnprotected = 0;
        var pauliX = new PauliXGate();

        for (int r = 0; r < rounds; r++)
        {
            // Generate random errors for this round
            var errors = new List<(QuantumGate gate, int qubit)>();
            for (int q = 0; q < code.PhysicalQubits; q++)
            {
                if (random.NextDouble() < errorRate)
                    errors.Add((pauliX, q));
            }

            // Protected: encode → error → correct → decode
            var result = RunBitFlipCorrection(code, initialState, errors, random);
            sumProtected += result.Fidelity;

            // Unprotected: apply same error(s) directly to the single-qubit state
            var unprotectedAmps = new ComplexVectorN(2);
            unprotectedAmps[0] = new ComplexNumber(initialState[0].realPart, initialState[0].imaginaryPart);
            unprotectedAmps[1] = new ComplexNumber(initialState[1].realPart, initialState[1].imaginaryPart);

            if (errors.Count > 0)
            {
                // For unprotected: a single X flip destroys the state
                unprotectedAmps = pauliX.Apply(unprotectedAmps, new[] { 0 }, 1);
            }

            double unprotectedFidelity = QuantumFidelity.Fidelity(initialState, unprotectedAmps);
            sumUnprotected += unprotectedFidelity;
        }

        return (sumProtected / rounds, sumUnprotected / rounds);
    }

    /// <summary>
    /// Runs a full error-correction cycle for Shor's 9-qubit code [[9,1,3]]:
    ///   1. Encode: outer phase-flip repetition + inner bit-flip repetition per block
    ///   2. Inject errors (any single-qubit X, Z, or Y)
    ///   3. Syndrome extraction: 8 ancilla qubits for 6 ZZ + 2 XXXXXX stabilizers
    ///   4. Decode + correct + uncompute encoding → fidelity
    ///
    /// Total qubits: 9 data + 8 ancilla = 17.
    /// </summary>
    public Result RunShorCorrection(
        ShorCode9 code,
        ComplexVectorN initialState,
        List<(QuantumGate gate, int qubit)> errorGates,
        Random random)
    {
        int totalQubits = code.PhysicalQubits + code.SyndromeQubits; // 17
        int stateSize = 1 << totalQubits;
        var amplitudes = new ComplexVectorN(stateSize);

        // Step 1: Prepare logical qubit on qubit 0
        amplitudes[0] = initialState[0]; // α for qubit 0 = |0⟩
        amplitudes[1] = initialState[1]; // β for qubit 0 = |1⟩

        var cnot = new CNOTGate();
        var hadamard = new HadamardGate();

        // Step 2: Encode — Shor code encoding
        // Phase-flip level: CNOT(0,3), CNOT(0,6)  → α|000⟩+β|111⟩ across block leaders
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 3 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 6 }, totalQubits);

        // Hadamard on block leaders: H(0), H(3), H(6)
        amplitudes = hadamard.Apply(amplitudes, new[] { 0 }, totalQubits);
        amplitudes = hadamard.Apply(amplitudes, new[] { 3 }, totalQubits);
        amplitudes = hadamard.Apply(amplitudes, new[] { 6 }, totalQubits);

        // Bit-flip level within each block:
        // Block 0: CNOT(0,1), CNOT(0,2)
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 1 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 2 }, totalQubits);
        // Block 1: CNOT(3,4), CNOT(3,5)
        amplitudes = cnot.Apply(amplitudes, new[] { 3, 4 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 3, 5 }, totalQubits);
        // Block 2: CNOT(6,7), CNOT(6,8)
        amplitudes = cnot.Apply(amplitudes, new[] { 6, 7 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 6, 8 }, totalQubits);

        // Step 3: Inject errors
        if (errorGates != null)
        {
            foreach (var (gate, qubit) in errorGates)
            {
                int[] targets;
                if (gate.QubitCount == 1)
                    targets = new[] { qubit };
                else
                    throw new ArgumentException("Only single-qubit error gates are supported.");
                amplitudes = gate.Apply(amplitudes, targets, totalQubits);
            }
        }

        // Step 4: Syndrome extraction
        // Ancilla qubits 9-14 for bit-flip stabilizers (ZZ within blocks)
        // Ancilla qubits 15-16 for phase-flip stabilizers (XXXXXX between blocks)

        // g₁: Z₀Z₁ → ancilla 9
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 9 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 1, 9 }, totalQubits);
        // g₂: Z₁Z₂ → ancilla 10
        amplitudes = cnot.Apply(amplitudes, new[] { 1, 10 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 2, 10 }, totalQubits);
        // g₃: Z₃Z₄ → ancilla 11
        amplitudes = cnot.Apply(amplitudes, new[] { 3, 11 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 4, 11 }, totalQubits);
        // g₄: Z₄Z₅ → ancilla 12
        amplitudes = cnot.Apply(amplitudes, new[] { 4, 12 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 5, 12 }, totalQubits);
        // g₅: Z₆Z₇ → ancilla 13
        amplitudes = cnot.Apply(amplitudes, new[] { 6, 13 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 7, 13 }, totalQubits);
        // g₆: Z₇Z₈ → ancilla 14
        amplitudes = cnot.Apply(amplitudes, new[] { 7, 14 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 8, 14 }, totalQubits);

        // g₇: X₀X₁X₂X₃X₄X₅ → ancilla 15
        // For X-stabilizer measurement: H on data qubits → CNOT → H back
        // Or equivalently: use controlled-X from data to ancilla in X-basis
        // Simpler: H(ancilla) → CNOT(ancilla, data) for each data qubit → H(ancilla)
        // Actually: to measure X_i, we can do CNOT(data_i, ancilla) in X-basis.
        // Standard approach: H on ancilla, then CX(ancilla, data_i) for each i, then H on ancilla, then measure ancilla.
        // But that's a multi-control. Simpler: H on all relevant data qubits, CNOT them to ancilla, H back.
        // Equivalent: for X-type stabilizer measurement of X₀X₁X₂X₃X₄X₅:
        //   Apply H to data qubits 0-5, CNOT each to ancilla 15, apply H back to data qubits 0-5

        // g₇: X₀X₁X₂X₃X₄X₅ → ancilla 15
        for (int q = 0; q <= 5; q++)
            amplitudes = hadamard.Apply(amplitudes, new[] { q }, totalQubits);
        for (int q = 0; q <= 5; q++)
            amplitudes = cnot.Apply(amplitudes, new[] { q, 15 }, totalQubits);
        for (int q = 0; q <= 5; q++)
            amplitudes = hadamard.Apply(amplitudes, new[] { q }, totalQubits);

        // g₈: X₃X₄X₅X₆X₇X₈ → ancilla 16
        for (int q = 3; q <= 8; q++)
            amplitudes = hadamard.Apply(amplitudes, new[] { q }, totalQubits);
        for (int q = 3; q <= 8; q++)
            amplitudes = cnot.Apply(amplitudes, new[] { q, 16 }, totalQubits);
        for (int q = 3; q <= 8; q++)
            amplitudes = hadamard.Apply(amplitudes, new[] { q }, totalQubits);

        // Step 5: Measure all 8 ancilla qubits
        var state = new QuantumState(amplitudes);
        int[] syndromes = new int[8];
        for (int i = 0; i < 8; i++)
            syndromes[i] = state.MeasureQubit(9 + i, random);

        int syndrome = 0;
        for (int i = 0; i < 8; i++)
            if (syndromes[i] == 1)
                syndrome |= (1 << i);

        // Step 6: Decode syndrome and apply corrections
        var decoder = new SyndromeDecoder(code);
        var corrections = decoder.Decode(syndrome);

        amplitudes = state.Amplitudes;
        foreach (var (qubit, pauli) in corrections)
        {
            var corrGate = GetPauliGate(pauli);
            amplitudes = corrGate.Apply(amplitudes, new[] { qubit }, totalQubits);
        }

        // Step 7: Decode (undo encoding) — reverse order of encoding
        // Undo bit-flip within blocks
        amplitudes = cnot.Apply(amplitudes, new[] { 6, 8 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 6, 7 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 3, 5 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 3, 4 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 2 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 1 }, totalQubits);

        // Undo Hadamard on block leaders
        amplitudes = hadamard.Apply(amplitudes, new[] { 6 }, totalQubits);
        amplitudes = hadamard.Apply(amplitudes, new[] { 3 }, totalQubits);
        amplitudes = hadamard.Apply(amplitudes, new[] { 0 }, totalQubits);

        // Undo phase-flip CNOTs
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 6 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 3 }, totalQubits);

        // Step 8: Extract logical qubit
        int ancillaBits = 0;
        for (int i = 0; i < 8; i++)
            if (syndromes[i] == 1)
                ancillaBits |= (1 << (9 + i));

        var recoveredAmplitudes = ExtractLogicalQubit(amplitudes, 0, ancillaBits, totalQubits);
        var recoveredState = new QuantumState(recoveredAmplitudes);
        double fidelity = QuantumFidelity.Fidelity(initialState, recoveredAmplitudes);

        return new Result
        {
            Fidelity = fidelity,
            Syndrome = syndrome,
            Corrections = corrections,
            RecoveredState = recoveredState
        };
    }

    /// <summary>
    /// Runs a full error-correction cycle for the Steane [[7,1,3]] code:
    ///   1. Encode: CSS encoding based on the [7,4,3] Hamming code
    ///   2. Inject errors (any single-qubit X, Z, or Y)
    ///   3. Syndrome extraction: 3 Z-stabilizer + 3 X-stabilizer measurements
    ///   4. Decode + correct + uncompute encoding → fidelity
    ///
    /// Total qubits: 7 data + 6 ancilla = 13.
    ///
    /// Encoding circuit (qubit 0 = |ψ⟩, qubits 1–6 = |0⟩):
    ///   Step 1: CNOT(0,1), CNOT(0,2) — spread logical value via coset representative
    ///   Step 2: H(4), H(5), H(6) — superposition on stabilizer info qubits
    ///   Step 3: CNOTs based on the [7,3,4] dual Hamming generator matrix G₂:
    ///           CNOT(4,1), CNOT(4,2), CNOT(4,3)
    ///           CNOT(5,0), CNOT(5,2), CNOT(5,3)
    ///           CNOT(6,0), CNOT(6,1), CNOT(6,3)
    ///
    /// Ancilla assignment:
    ///   7 → Z₃Z₄Z₅Z₆ (syndrome bit 0)
    ///   8 → Z₁Z₂Z₅Z₆ (syndrome bit 1)
    ///   9 → Z₀Z₂Z₄Z₆ (syndrome bit 2)
    ///  10 → X₃X₄X₅X₆ (syndrome bit 3)
    ///  11 → X₁X₂X₅X₆ (syndrome bit 4)
    ///  12 → X₀X₂X₄X₆ (syndrome bit 5)
    /// </summary>
    public Result RunSteaneCorrection(
        SteaneCode7 code,
        ComplexVectorN initialState,
        List<(QuantumGate gate, int qubit)> errorGates,
        Random random)
    {
        int totalQubits = code.PhysicalQubits + code.SyndromeQubits; // 13
        int stateSize = 1 << totalQubits;
        var amplitudes = new ComplexVectorN(stateSize);

        // Prepare logical qubit on qubit 0
        amplitudes[0] = initialState[0]; // α for |0⟩
        amplitudes[1] = initialState[1]; // β for |1⟩

        var cnot = new CNOTGate();
        var hadamard = new HadamardGate();

        // ── Encode ──────────────────────────────────────────────────────
        // Step 1: Spread logical value via coset representative v = [1,1,1,0,0,0,0]
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 1 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 2 }, totalQubits);

        // Step 2: Hadamard on stabilizer info qubits (4, 5, 6)
        amplitudes = hadamard.Apply(amplitudes, new[] { 4 }, totalQubits);
        amplitudes = hadamard.Apply(amplitudes, new[] { 5 }, totalQubits);
        amplitudes = hadamard.Apply(amplitudes, new[] { 6 }, totalQubits);

        // Step 3: Entangle via [7,3,4] dual Hamming generator matrix G₂
        // G₂ row for q₄: adds to positions 1,2,3
        amplitudes = cnot.Apply(amplitudes, new[] { 4, 1 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 4, 2 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 4, 3 }, totalQubits);
        // G₂ row for q₅: adds to positions 0,2,3
        amplitudes = cnot.Apply(amplitudes, new[] { 5, 0 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 5, 2 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 5, 3 }, totalQubits);
        // G₂ row for q₆: adds to positions 0,1,3
        amplitudes = cnot.Apply(amplitudes, new[] { 6, 0 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 6, 1 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 6, 3 }, totalQubits);

        // ── Inject errors ───────────────────────────────────────────────
        if (errorGates != null)
        {
            foreach (var (gate, qubit) in errorGates)
            {
                amplitudes = gate.Apply(amplitudes, new[] { qubit }, totalQubits);
            }
        }

        // ── Syndrome extraction ─────────────────────────────────────────
        // Z-stabilizer measurement (detect X errors): CNOT from data to ancilla
        // bit 0: Z₃Z₄Z₅Z₆ → ancilla 7
        amplitudes = cnot.Apply(amplitudes, new[] { 3, 7 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 4, 7 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 5, 7 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 6, 7 }, totalQubits);
        // bit 1: Z₁Z₂Z₅Z₆ → ancilla 8
        amplitudes = cnot.Apply(amplitudes, new[] { 1, 8 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 2, 8 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 5, 8 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 6, 8 }, totalQubits);
        // bit 2: Z₀Z₂Z₄Z₆ → ancilla 9
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 9 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 2, 9 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 4, 9 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 6, 9 }, totalQubits);

        // X-stabilizer measurement (detect Z errors): H → CNOT → H
        // bit 3: X₃X₄X₅X₆ → ancilla 10
        int[] xStab0 = { 3, 4, 5, 6 };
        foreach (int q in xStab0)
            amplitudes = hadamard.Apply(amplitudes, new[] { q }, totalQubits);
        foreach (int q in xStab0)
            amplitudes = cnot.Apply(amplitudes, new[] { q, 10 }, totalQubits);
        foreach (int q in xStab0)
            amplitudes = hadamard.Apply(amplitudes, new[] { q }, totalQubits);

        // bit 4: X₁X₂X₅X₆ → ancilla 11
        int[] xStab1 = { 1, 2, 5, 6 };
        foreach (int q in xStab1)
            amplitudes = hadamard.Apply(amplitudes, new[] { q }, totalQubits);
        foreach (int q in xStab1)
            amplitudes = cnot.Apply(amplitudes, new[] { q, 11 }, totalQubits);
        foreach (int q in xStab1)
            amplitudes = hadamard.Apply(amplitudes, new[] { q }, totalQubits);

        // bit 5: X₀X₂X₄X₆ → ancilla 12
        int[] xStab2 = { 0, 2, 4, 6 };
        foreach (int q in xStab2)
            amplitudes = hadamard.Apply(amplitudes, new[] { q }, totalQubits);
        foreach (int q in xStab2)
            amplitudes = cnot.Apply(amplitudes, new[] { q, 12 }, totalQubits);
        foreach (int q in xStab2)
            amplitudes = hadamard.Apply(amplitudes, new[] { q }, totalQubits);

        // ── Measure ancilla qubits ──────────────────────────────────────
        var state = new QuantumState(amplitudes);
        int[] syndromes = new int[6];
        for (int i = 0; i < 6; i++)
            syndromes[i] = state.MeasureQubit(7 + i, random);

        int syndrome = 0;
        for (int i = 0; i < 6; i++)
            if (syndromes[i] == 1)
                syndrome |= (1 << i);

        // ── Decode syndrome and apply corrections ───────────────────────
        var decoder = new SyndromeDecoder(code);
        var corrections = decoder.Decode(syndrome);

        amplitudes = state.Amplitudes;
        foreach (var (qubit, pauli) in corrections)
        {
            var corrGate = GetPauliGate(pauli);
            amplitudes = corrGate.Apply(amplitudes, new[] { qubit }, totalQubits);
        }

        // ── Decode (undo encoding in reverse order) ─────────────────────
        // Undo Step 3: reverse G₂ CNOTs
        amplitudes = cnot.Apply(amplitudes, new[] { 6, 3 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 6, 1 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 6, 0 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 5, 3 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 5, 2 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 5, 0 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 4, 3 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 4, 2 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 4, 1 }, totalQubits);

        // Undo Step 2: H on qubits 4, 5, 6
        amplitudes = hadamard.Apply(amplitudes, new[] { 6 }, totalQubits);
        amplitudes = hadamard.Apply(amplitudes, new[] { 5 }, totalQubits);
        amplitudes = hadamard.Apply(amplitudes, new[] { 4 }, totalQubits);

        // Undo Step 1: CNOT(0,2), CNOT(0,1)
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 2 }, totalQubits);
        amplitudes = cnot.Apply(amplitudes, new[] { 0, 1 }, totalQubits);

        // ── Extract logical qubit ───────────────────────────────────────
        int ancillaBits = 0;
        for (int i = 0; i < 6; i++)
            if (syndromes[i] == 1)
                ancillaBits |= (1 << (7 + i));

        var recoveredAmplitudes = ExtractLogicalQubit(amplitudes, 0, ancillaBits, totalQubits);
        var recoveredState = new QuantumState(recoveredAmplitudes);
        double fidelity = QuantumFidelity.Fidelity(initialState, recoveredAmplitudes);

        return new Result
        {
            Fidelity = fidelity,
            Syndrome = syndrome,
            Corrections = corrections,
            RecoveredState = recoveredState
        };
    }

    /// <summary>
    /// Extracts the logical qubit amplitudes after decoding. After a successful QEC
    /// cycle, data qubits 1,2 should be |0⟩ and ancilla qubits are collapsed to their
    /// measured values. We extract the two amplitudes for the logical qubit (qubit 0)
    /// from the full state vector at the known ancilla+idle qubit configuration.
    /// </summary>
    /// <param name="amplitudes">Full state vector.</param>
    /// <param name="logicalQubit">Index of the logical qubit (0).</param>
    /// <param name="knownBits">Bit pattern for the known qubit values (ancilla bits).</param>
    /// <param name="totalQubits">Total number of qubits.</param>
    private static ComplexVectorN ExtractLogicalQubit(
        ComplexVectorN amplitudes, int logicalQubit, int knownBits, int totalQubits)
    {
        // After decoding, the state should be |ψ⟩_logical ⊗ |0⟩₁ ⊗ |0⟩₂ ⊗ |measured ancilla⟩
        // The two relevant indices are:
        //   idx0 = knownBits (logical qubit = 0, data 1,2 = 0, ancilla = measured)
        //   idx1 = knownBits | (1 << logicalQubit) (logical qubit = 1, rest same)
        int idx0 = knownBits;
        int idx1 = knownBits | (1 << logicalQubit);

        var amp0 = amplitudes[idx0];
        var amp1 = amplitudes[idx1];

        // Renormalize
        double norm = amp0.realPart * amp0.realPart + amp0.imaginaryPart * amp0.imaginaryPart
                    + amp1.realPart * amp1.realPart + amp1.imaginaryPart * amp1.imaginaryPart;

        var result = new ComplexVectorN(2);
        if (norm > 1e-15)
        {
            double scale = 1.0 / Math.Sqrt(norm);
            result[0] = new ComplexNumber(amp0.realPart * scale, amp0.imaginaryPart * scale);
            result[1] = new ComplexNumber(amp1.realPart * scale, amp1.imaginaryPart * scale);
        }
        return result;
    }

    /// <summary>
    /// Returns the appropriate Pauli gate for a given Pauli character.
    /// </summary>
    private static QuantumGate GetPauliGate(char pauli)
    {
        switch (pauli)
        {
            case 'X': return new PauliXGate();
            case 'Y': return new PauliYGate();
            case 'Z': return new PauliZGate();
            default: throw new ArgumentException($"Unknown Pauli operator: {pauli}");
        }
    }
}
