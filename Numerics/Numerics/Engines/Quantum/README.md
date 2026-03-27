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
| `H(q)`, `X(q)`, `Z(q)`, `S(q)`, `T(q)` | Single-qubit standard gates |
| `Rx(q, θ)`, `Ry(q, θ)`, `Rz(q, θ)` | Parameterised rotation gates |
| `CNOT(c, t)`, `CZ(c, t)`, `SWAP(a, b)` | Two-qubit gates |
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
| `PauliZGate` | 1 | Phase flip: Z\|1⟩ = −\|1⟩ |
| `SGate` | 1 | π/2 phase: S\|1⟩ = i\|1⟩, $S^2 = Z$ |
| `TGate` | 1 | π/4 phase: T\|1⟩ = $e^{i\pi/4}$\|1⟩, $T^2 = S$ |

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
| `SWAPGate` | 2 | Swaps the states of two qubits |

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

### Module Structure

```
Engines/Quantum/
├── QuantumCircuit.cs          Circuit definition (qubits + instruction list)
├── QuantumCircuitBuilder.cs   Fluent builder API
├── QuantumInstruction.cs      Gate + target qubit binding
├── QuantumSimulator.cs        Stateless ideal circuit executor
├── NoisyQuantumSimulator.cs   Noisy simulator (Monte Carlo trajectories)
├── QuantumEnvironment.cs      RL environment for circuit synthesis
├── QuantumState.cs            Result: amplitudes, measurement, sampling
└── README.md
```
