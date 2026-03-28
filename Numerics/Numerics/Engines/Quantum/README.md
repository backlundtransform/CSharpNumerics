## ⚛️ Quantum Circuits Engine

A lightweight quantum circuit simulator that builds, composes, and executes quantum circuits using the gate primitives from `CSharpNumerics.Physics.Quantum`. State vectors use `ComplexVectorN` for amplitudes and `VectorN` for probability distributions.

**Namespace:** `CSharpNumerics.Engines.Quantum`

---

### Quick Start

```csharp
using CSharpNumerics.Physics.Quantum;
using CSharpNumerics.Engines.Quantum;

// Fluent API — build a Bell state in one chain
var circuit = QuantumCircuitBuilder.New(2)
    .H(0)
    .CNOT(0, 1)
    .Build();

// Simulate
var state = new QuantumSimulator().Run(circuit);

// Inspect probabilities
double p00 = state.GetProbability(0);   // 0.5  — |00⟩
double p11 = state.GetProbability(3);   // 0.5  — |11⟩

VectorN probs = state.GetProbabilities(); // [0.5, 0, 0, 0.5]

// Measurement — collapse to a single outcome
int result = state.Measure(new Random());

// Sampling — statistical distribution without collapse
var counts = state.Sample(shots: 10000, new Random());
```

---

### QuantumCircuit

Holds a fixed qubit count and an ordered list of instructions.

```csharp
var circuit = new QuantumCircuit(3);          // 3-qubit register
circuit.AddInstruction(instruction);          // append gate + target qubits
int n = circuit.QubitCount;                   // 3
var list = circuit.Instructions;              // List<QuantumInstruction>
```

| Property / Method | Type | Description |
|---|---|---|
| `QubitCount` | `int` | Number of qubits in the register |
| `Instructions` | `List<QuantumInstruction>` | Ordered gate sequence |
| `AddInstruction(instruction)` | `void` | Appends an instruction (validates qubit indices) |

---

### QuantumInstruction

Pairs a `QuantumGate` with the physical qubit indices it targets.

```csharp
var instr = new QuantumInstruction(new HadamardGate(), new List<int> { 0 });
var cnot  = new QuantumInstruction(new CNOTGate(), new List<int> { 0, 1 }); // control=0, target=1
```

| Property | Type | Description |
|---|---|---|
| `Gate` | `QuantumGate` | The gate to apply |
| `QubitIndices` | `List<int>` | Physical qubit indices (count must match `Gate.QubitCount`) |

---

### QuantumCircuitBuilder

Fluent API for building circuits without manual instruction wiring.

```csharp
var circuit = QuantumCircuitBuilder.New(3)
    .H(0)                   // Hadamard on qubit 0
    .X(1)                   // Pauli-X on qubit 1
    .Ry(2, Math.PI / 4)     // Y-rotation on qubit 2
    .CNOT(0, 1)             // CNOT: control=0, target=1
    .CZ(1, 2)               // Controlled-Z
    .SWAP(0, 2)             // Swap qubits 0 and 2
    .Gate(new TGate(), 1)   // Arbitrary gate
    .Build();
```

| Method | Description |
|---|---|
| `New(int qubits)` | Create a builder for an n-qubit circuit |
| `H(q)`, `X(q)`, `Y(q)`, `Z(q)`, `S(q)`, `T(q)` | Single-qubit standard gates |
| `Phase(q, θ)` | General phase gate $P(\theta)$ |
| `Rx(q, θ)`, `Ry(q, θ)`, `Rz(q, θ)` | Parameterised rotation gates |
| `CNOT(c, t)`, `CZ(c, t)`, `CPhase(c, t, θ)`, `SWAP(a, b)` | Two-qubit gates |
| `Toffoli(c1, c2, t)` | Three-qubit CCNOT gate |
| `Fredkin(c, t1, t2)` | Three-qubit CSWAP gate |
| `ApplyQFT(params qubits)` | Quantum Fourier Transform on specified qubits |
| `ApplyInverseQFT(params qubits)` | Inverse QFT (QFT†) on specified qubits |
| `ApplyGrover(marked, params qubits)` | Grover search with optimal iterations |
| `ApplyGrover(marked, iters, params qubits)` | Grover search with explicit iteration count |
| `ApplyQPE(counting[], target[], gate)` | Quantum Phase Estimation |
| `Controlled(gate, control, params targets)` | Controlled version of any gate |
| `Gate(gate, params qubits)` | Arbitrary gate on specified qubits |
| `Build()` | Returns the constructed `QuantumCircuit` |

---

### QuantumState

Result of a simulation — a complex amplitude vector with probability queries, measurement, and sampling.

```csharp
double p = state.GetProbability(0);          // |α₀|² for basis state |0…0⟩
VectorN probs = state.GetProbabilities();    // all |αᵢ|² as VectorN
int qubits = state.QubitCount;               // number of qubits
ComplexVectorN amps = state.Amplitudes;      // raw amplitude vector
```

| Property / Method | Type | Description |
|---|---|---|
| `Amplitudes` | `ComplexVectorN` | Complex amplitudes (length 2ⁿ) |
| `QubitCount` | `int` | Number of qubits |
| `GetProbability(int basisState)` | `double` | Measurement probability of a specific basis state |
| `GetProbabilities()` | `VectorN` | All basis-state probabilities |
| `GetBlochVector()` | `BlochVector` | Bloch sphere coordinates (single-qubit only) |
| `Measure(Random)` | `int` | Full measurement — collapses state, returns basis index |
| `MeasureQubit(int, Random)` | `int` | Measure one qubit (0 or 1), renormalising the rest |
| `Sample(int shots, Random)` | `Dictionary<int,int>` | Sample distribution without collapse |

**Measurement (collapse)**

```csharp
var state = sim.Run(QuantumCircuitBuilder.New(1).H(0).Build());
var rng = new Random(42);

// Full measurement — collapses to |0⟩ or |1⟩
int result = state.Measure(rng);  // 0 or 1
// State is now a basis state: P(result) = 1.0
```

**Single-qubit measurement**

```csharp
// Bell state: (|00⟩+|11⟩)/√2
var state = sim.Run(QuantumCircuitBuilder.New(2).H(0).CNOT(0, 1).Build());

// Measure only qubit 0 — collapses both qubits due to entanglement
int q0 = state.MeasureQubit(0, rng);  // 0 or 1
```

**Sampling (no collapse)**

```csharp
// Run 10000 shots — state is NOT modified
var counts = state.Sample(10000, rng);
// counts = { 0: ~5000, 1: ~5000 } for H|0⟩
```

---

### QuantumSimulator

Stateless simulator: initialises |0…0⟩ and applies each instruction in sequence.

```csharp
var sim = new QuantumSimulator();
QuantumState result = sim.Run(circuit);
```

---

### NoisyQuantumSimulator

Simulator with quantum noise channels applied after each gate using **Monte Carlo trajectory selection** (stochastic Kraus operator sampling). Noise channels come from `CSharpNumerics.Physics.Quantum.NoiseModels`.

```csharp
using CSharpNumerics.Physics.Quantum.NoiseModels;

var noisy = new NoisyQuantumSimulator(new Random(42))
    .WithNoise(new DepolarizingNoise(0.01))      // 1% depolarizing
    .WithNoise(new DephasingNoise(0.02))          // 2% dephasing
    .Run(circuit);

double fidelity = QuantumFidelity.Fidelity(idealState, noisy);
```

| Method | Description |
|---|---|
| `NoisyQuantumSimulator(Random)` | Constructor with RNG for trajectory sampling |
| `WithNoise(INoiseChannel)` | Adds a noise channel (fluent, stackable) |
| `Run(QuantumCircuit)` | Executes with noise after each gate |

**How it works:** After each gate, every registered noise channel is applied independently to each qubit the gate acted on. For each channel the simulator computes $p_k = \langle\psi|E_k^\dagger E_k|\psi\rangle$ for all Kraus operators $\{E_k\}$, samples one operator, and renormalises.

**Noise channels** (from `CSharpNumerics.Physics.Quantum.NoiseModels`):

| Channel | Parameter | Effect |
|---|---|---|
| `DepolarizingNoise(p)` | $p \in [0,1]$ | Replaces qubit with maximally mixed state with probability $p$ |
| `DephasingNoise(p)` | $p \in [0,1]$ | Randomises relative phase (T₂ decoherence) |
| `AmplitudeDampingNoise(γ)` | $\gamma \in [0,1]$ | Energy dissipation: \|1⟩ decays to \|0⟩ (T₁ decay) |

All channels implement `INoiseChannel` and provide Kraus operators satisfying $\sum E_k^\dagger E_k = I$.

---

### Available Gates (from `CSharpNumerics.Physics.Quantum`)

#### Single-Qubit Gates

| Gate | Qubits | Effect |
|---|---|---|
| `HadamardGate` | 1 | Equal superposition: H\|0⟩ = (\|0⟩ + \|1⟩)/√2 |
| `PauliXGate` | 1 | Bit flip: X\|0⟩ = \|1⟩ |
| `PauliYGate` | 1 | Bit + phase flip: $Y^2 = I$ |
| `PauliZGate` | 1 | Phase flip: Z\|1⟩ = −\|1⟩ |
| `SGate` | 1 | π/2 phase: S\|1⟩ = i\|1⟩, $S^2 = Z$ |
| `TGate` | 1 | π/4 phase: T\|1⟩ = $e^{i\pi/4}$\|1⟩, $T^2 = S$ |
| `PhaseGate(θ)` | 1 | General phase: $P(\theta)\|1\rangle = e^{i\theta}\|1\rangle$, $P(\pi) = Z$ |

#### Rotation Gates

| Gate | Qubits | Effect |
|---|---|---|
| `RxGate(θ)` | 1 | Rotation about X-axis by θ |
| `RyGate(θ)` | 1 | Rotation about Y-axis by θ |
| `RzGate(θ)` | 1 | Rotation about Z-axis by θ |

#### Two-Qubit Gates

| Gate | Qubits | Effect |
|---|---|---|
| `CNOTGate` | 2 | Controlled-NOT: flips target when control = \|1⟩ |
| `CZGate` | 2 | Controlled-Z: phase flip on \|11⟩ |
| `CPhaseGate(θ)` | 2 | Controlled phase: $e^{i\theta}$ on \|11⟩, $CP(\pi) = CZ$ |
| `SWAPGate` | 2 | Swaps the states of two qubits |

#### Three-Qubit Gates

| Gate | Qubits | Effect |
|---|---|---|
| `ToffoliGate` | 3 | CCNOT: flips target when both controls = \|1⟩ |
| `FredkinGate` | 3 | CSWAP: swaps targets when control = \|1⟩ |

#### N-Qubit Gates

| Gate | Qubits | Effect |
|---|---|---|
| `PhaseOracle(n, states)` | n | Flips phase of marked basis states |
| `ControlledGate(U)` | n+1 | Applies U when control = \|1⟩ (wraps any gate) |
| `ModularMultiplyGate(a,N,n)` | n | Permutation: \|y⟩ → \|ay mod N⟩ |

---

### Examples

**Rotation circuit**

```csharp
var circuit = new QuantumCircuit(1);
circuit.AddInstruction(new QuantumInstruction(new RyGate(Math.PI / 2), new List<int> { 0 }));
// Equivalent to Hadamard on |0⟩ — equal superposition

var state = new QuantumSimulator().Run(circuit);
// P(|0⟩) ≈ 0.5, P(|1⟩) ≈ 0.5
```

**3-qubit uniform superposition**

```csharp
var circuit = new QuantumCircuit(3);
circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));
circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 1 }));
circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 2 }));

var state = new QuantumSimulator().Run(circuit);
// Each of the 8 basis states has probability 1/8
```

**Bit-flip verification (X² = I)**

```csharp
var circuit = new QuantumCircuit(1);
circuit.AddInstruction(new QuantumInstruction(new PauliXGate(), new List<int> { 0 }));
circuit.AddInstruction(new QuantumInstruction(new PauliXGate(), new List<int> { 0 }));

var state = new QuantumSimulator().Run(circuit);
// Back to |0⟩ with probability 1.0
```

---

### QuantumEnvironment (RL Integration)

An `IEnvironment` implementation for reinforcement learning — the agent learns a gate sequence to prepare a target quantum state from |0…0⟩.

| Property | Description |
|---|---|
| **Observation** | State probabilities (2ⁿ values) + normalised step progress |
| **Actions** | Discrete — each index maps to a (gate, qubit) pair |
| **Reward** | Fidelity delta + bonus when threshold reached |
| **Done** | Max gates reached or fidelity ≥ threshold |

**Quick start**

```csharp
using CSharpNumerics.Engines.Quantum;
using CSharpNumerics.Physics.Quantum;

// Target: Bell state (|00⟩+|11⟩)/√2
var bellCircuit = QuantumCircuitBuilder.New(2).H(0).CNOT(0, 1).Build();

var env = QuantumEnvironment.Create(2)
    .WithTargetCircuit(bellCircuit)
    .WithMaxGates(20)
    .WithFidelityThreshold(0.99)
    .Build();

// Use with the RL framework
var result = RLExperiment
    .For(env)
    .WithAgent(agent)
    .WithPolicy(new EpsilonGreedy())
    .WithEpisodes(1000)
    .Run();
```

**Builder methods**

| Method | Description |
|---|---|
| `Create(int qubits)` | Start building an environment |
| `.WithTargetState(QuantumState)` | Set the target state directly |
| `.WithTargetCircuit(QuantumCircuit)` | Set the target by simulating a circuit |
| `.WithMaxGates(int)` | Maximum gates per episode (default 20) |
| `.WithFidelityThreshold(double)` | Early termination threshold (default 0.99) |
| `.WithActions(List<QuantumInstruction>)` | Custom action set |
| `.Build()` | Construct the environment |

**Default action set** (auto-generated when no custom actions provided):
- Single-qubit gates: H, X, Z, S, T on each qubit (5 × n)
- Two-qubit gates: CNOT on all ordered qubit pairs (n × (n−1))

**Step info dictionary**

| Key | Type | Description |
|---|---|---|
| `"fidelity"` | `double` | Current fidelity to target state |
| `"gates"` | `int` | Number of gates placed so far |

---

### Algorithms

**Namespace:** `CSharpNumerics.Engines.Quantum.Algorithms`

#### Quantum Fourier Transform (QFT)

Generates a circuit implementing the standard QFT: Hadamard gates, controlled-phase rotations $CP(\pi/2^k)$, and a SWAP cascade to reverse qubit order.

```csharp
using CSharpNumerics.Engines.Quantum.Algorithms;

// Standalone — returns a QuantumCircuit
var qftCircuit = QFT.CreateCircuit(totalQubits: 3, 0, 1, 2);
var state = new QuantumSimulator().Run(qftCircuit);

// Via builder — inline into a larger circuit
var circuit = QuantumCircuitBuilder.New(4)
    .X(0).X(2)                   // Prepare |0101⟩
    .ApplyQFT(0, 1, 2, 3)       // QFT on all 4 qubits
    .Build();
```

QFT on $|0\rangle^{\otimes n}$ produces a uniform superposition: each basis state has probability $1/2^n$.

#### Inverse QFT (QFT†)


Reverses the QFT: SWAP cascade, then reversed controlled-phase rotations with negated angles and Hadamards.

```csharp
// Round-trip: InverseQFT(QFT(|ψ⟩)) = |ψ⟩
var circuit = QuantumCircuitBuilder.New(3)
    .X(0).X(2)                        // |101⟩
    .ApplyQFT(0, 1, 2)
    .ApplyInverseQFT(0, 1, 2)
    .Build();

var state = new QuantumSimulator().Run(circuit);
// state ≡ |101⟩ with fidelity ≈ 1.0
```

#### Grover's Search Algorithm

Amplifies the amplitude of marked basis states in an unstructured search space. Given M marked states in N = 2ⁿ, the algorithm finds a marked state with high probability after ≈ (π/4)√(N/M) iterations.

```csharp
using CSharpNumerics.Engines.Quantum.Algorithms;

// Standalone — search for |101⟩ in a 3-qubit space
var circuit = GroverSearch.CreateCircuit(
    totalQubits: 3,
    searchQubits: new[] { 0, 1, 2 },
    markedStates: new[] { 5 });         // |101⟩ = index 5

var state = new QuantumSimulator().Run(circuit);
double pTarget = state.GetProbability(5); // > 0.94

// Via builder
var circuit2 = QuantumCircuitBuilder.New(3)
    .ApplyGrover(new[] { 5 }, new[] { 0, 1, 2 })
    .Build();

// Multiple targets
var circuit3 = GroverSearch.CreateCircuit(3,
    new[] { 0, 1, 2 }, new[] { 3, 5 }); // find |011⟩ or |101⟩

// Query optimal iteration count
int iters = GroverSearch.OptimalIterations(
    searchSpaceQubits: 3, markedCount: 1); // 2
```

The algorithm uses a `PhaseOracle` gate (from `CSharpNumerics.Physics.Quantum`) — an n-qubit diagonal unitary that flips the phase of specified basis states.

#### Quantum Phase Estimation (QPE)

Estimates the eigenvalue phase φ of a unitary gate U, where $U|\psi\rangle = e^{2\pi i\varphi}|\psi\rangle$. Uses t counting qubits for t bits of precision.

```csharp
using CSharpNumerics.Engines.Quantum.Algorithms;

// Standalone — estimate phase of PhaseGate(π/2), eigenstate |1⟩
// φ = 1/4, with 2 counting qubits → measures 1 (= φ·2²)
int totalQubits = 3;
var circuit = new QuantumCircuit(totalQubits);
circuit.AddInstruction(new QuantumInstruction(new PauliXGate(), new List<int> { 2 })); // |1⟩ eigenstate
var qpe = QPE.CreateCircuit(totalQubits,
    countingQubits: new[] { 0, 1 },
    targetQubits: new[] { 2 },
    unitaryGate: new PhaseGate(Math.PI / 2));
foreach (var instr in qpe.Instructions)
    circuit.AddInstruction(instr);

var state = new QuantumSimulator().Run(circuit);
// Counting register encodes φ·2^t = 1

// Via builder
var circuit2 = QuantumCircuitBuilder.New(3)
    .X(2)
    .ApplyQPE(new[] { 0, 1 }, new[] { 2 }, new PhaseGate(Math.PI / 2))
    .Build();
```

QPE uses `ControlledGate` (from `CSharpNumerics.Physics.Quantum`) to create controlled-U operations automatically. The `countingQubits[0]` is the MSB of the result.

#### Shor's Factoring Algorithm

Factors composite integers using quantum order-finding (QPE + modular multiplication).

```csharp
using CSharpNumerics.Engines.Quantum.Algorithms;

// Factor 15 = 3 × 5
var result = ShorAlgorithm.Factor(15, new Random(42));
Console.WriteLine($"{result.Factor1} × {result.Factor2}"); // 3 × 5
Console.WriteLine($"Success: {result.Success}");           // True
Console.WriteLine($"Base: {result.Base}, Order: {result.Order}");

// Access the quantum circuit directly
var circuit = ShorAlgorithm.CreateOrderFindingCircuit(
    a: 2, N: 15, countingQubits: 8, targetQubits: 4);

// Helper utilities
int iters = GroverSearch.OptimalIterations(3);              // Grover
long order = ShorAlgorithm.ExtractOrder(4, 4, 15, 2);      // continued fractions
long g = ShorAlgorithm.GCD(12, 8);                         // 4
long p = ShorAlgorithm.ModPow(2, 10, 1000);                // 1024 mod 1000 = 24
```

The quantum core uses `ModularMultiplyGate` (from `CSharpNumerics.Physics.Quantum`) — a permutation gate implementing $|y\rangle \to |ay \bmod N\rangle$.

| Class | Role |
|---|---|
| `ShorAlgorithm` | Full factor algorithm (classical + quantum) |
| `ShorResult` | Result with factors, base, order, attempt count |
| `ModularMultiplyGate` | Permutation oracle $U_a\|y\rangle = \|ay \bmod N\rangle$ |
| `ControlledGate` | Wraps any gate with a control qubit |

---

### Module Structure

```
Engines/Quantum/
├── Algorithms/
│   ├── QFT.cs                 Quantum Fourier Transform circuit generator
│   ├── InverseQFT.cs          Inverse QFT circuit generator
│   ├── GroverSearch.cs        Grover's search algorithm
│   ├── QPE.cs                 Quantum Phase Estimation
│   └── ShorAlgorithm.cs       Shor's factoring algorithm
├── ErrorCorrection/
│   ├── SyndromeDecoder.cs     Classical syndrome → correction lookup
│   └── ErrorCorrectionSimulator.cs  Full QEC simulation orchestrator
├── QuantumCircuit.cs          Circuit definition (qubits + instruction list)
├── QuantumCircuitBuilder.cs   Fluent builder API
├── QuantumInstruction.cs      Gate + target qubit binding
├── QuantumSimulator.cs        Stateless ideal circuit executor
├── NoisyQuantumSimulator.cs   Noisy simulator (Monte Carlo trajectories)
├── QuantumEnvironment.cs      RL environment for circuit synthesis
├── QuantumState.cs            Result: amplitudes, measurement, sampling
└── README.md
```

### Quantum Error Correction — Simulation

**Namespace:** `CSharpNumerics.Engines.Quantum.ErrorCorrection`

The `ErrorCorrection` sub-namespace provides the simulation layer for quantum error-correcting codes defined in `Physics.Quantum.ErrorCorrection`.

#### SyndromeDecoder

Classical lookup table that maps measured syndrome bits to corrective operations:

```csharp
using CSharpNumerics.Physics.Quantum.ErrorCorrection;
using CSharpNumerics.Engines.Quantum.ErrorCorrection;

var decoder = new SyndromeDecoder(new BitFlipCode3());

// From integer syndrome
var ops = decoder.Decode(0b01);     // → [(0, 'X')]  — apply X to qubit 0

// From measurement bit array
var ops2 = decoder.Decode(new[] { 1, 1 });  // → [(1, 'X')]  — apply X to qubit 1
```

#### ErrorCorrectionSimulator

Orchestrates the full QEC cycle: **encode → inject error → extract syndrome → decode → correct → decode → measure fidelity**.

**Bit-flip correction:**

```csharp
using CSharpNumerics.Physics.Quantum;
using CSharpNumerics.Physics.Quantum.ErrorCorrection;
using CSharpNumerics.Engines.Quantum.ErrorCorrection;
using CSharpNumerics.Numerics.Objects;

var code = new BitFlipCode3();
var sim  = new ErrorCorrectionSimulator();
var rng  = new Random(42);

// Prepare a superposition state to protect
var initial = new ComplexVectorN(2);
initial[0] = new ComplexNumber(0.6, 0);   // α
initial[1] = new ComplexNumber(0, 0.8);   // β

// Inject a single bit-flip error on qubit 1
var errors = new List<(QuantumGate, int)> { (new PauliXGate(), 1) };

var result = sim.RunBitFlipCorrection(code, initial, errors, rng);

Console.WriteLine(result.Fidelity);   // ≈ 1.0 — error fully corrected
Console.WriteLine(result.Syndrome);   // 3 (0b11) — both stabilizers triggered
```

**Phase-flip correction:**

```csharp
var pfCode = new PhaseFlipCode3();
var zError = new List<(QuantumGate, int)> { (new PauliZGate(), 0) };

var pfResult = sim.RunPhaseFlipCorrection(pfCode, initial, zError, rng);
// pfResult.Fidelity ≈ 1.0
```

**Monte Carlo comparison — protected vs. unprotected:**

```csharp
var (protectedFidelity, unprotectedFidelity) = sim.RunMonteCarloComparison(
    code, initial,
    errorRate: 0.10,  // 10% bit-flip probability per qubit
    rounds: 1000,
    random: rng);

// protectedFidelity > unprotectedFidelity — QEC provides net benefit
```

**Shor 9-qubit code — corrects any single-qubit error (X, Z, or Y):**

```csharp
using CSharpNumerics.Physics.Quantum;
using CSharpNumerics.Physics.Quantum.ErrorCorrection;
using CSharpNumerics.Engines.Quantum.ErrorCorrection;
using CSharpNumerics.Numerics.Objects;

var shor = new ShorCode9();
var sim  = new ErrorCorrectionSimulator();
var rng  = new Random(42);

var initial = new ComplexVectorN(2);
initial[0] = new ComplexNumber(0.6, 0);
initial[1] = new ComplexNumber(0, 0.8);

// Y error = combined bit-flip + phase-flip
var errors = new List<(QuantumGate, int)> { (new PauliYGate(), 4) };
var result = sim.RunShorCorrection(shor, initial, errors, rng);
// result.Fidelity ≈ 1.0 — Y error fully corrected
```

**Steane 7-qubit code — the smallest CSS code correcting any single-qubit error:**

```csharp
using CSharpNumerics.Physics.Quantum;
using CSharpNumerics.Physics.Quantum.ErrorCorrection;
using CSharpNumerics.Engines.Quantum.ErrorCorrection;
using CSharpNumerics.Numerics.Objects;

var steane = new SteaneCode7();
var sim    = new ErrorCorrectionSimulator();
var rng    = new Random(42);

var initial = new ComplexVectorN(2);
initial[0] = new ComplexNumber(0.6, 0);
initial[1] = new ComplexNumber(0, 0.8);

// Z (phase-flip) error on qubit 3
var errors = new List<(QuantumGate, int)> { (new PauliZGate(), 3) };
var result = sim.RunSteaneCorrection(steane, initial, errors, rng);
// result.Fidelity ≈ 1.0 — CSS structure detects and corrects Z error
// result.Syndrome bits 3-5 identify the Z error via Hamming decoding
```

| Return Type | Fields |
|---|---|
| `Result` | `Fidelity` (double), `Syndrome` (int), `Corrections` (list), `RecoveredState` (QuantumState) |
| Monte Carlo | `(protectedFidelity, unprotectedFidelity)` tuple |
