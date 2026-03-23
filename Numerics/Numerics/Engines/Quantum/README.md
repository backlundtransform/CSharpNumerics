## ⚛️ Quantum Circuits Engine

A lightweight quantum circuit simulator that builds, composes, and executes quantum circuits using the gate primitives from `CSharpNumerics.Physics.Quantum`. State vectors use `ComplexVectorN` for amplitudes and `VectorN` for probability distributions.

**Namespace:** `CSharpNumerics.Engines.Quantum`

---

### Quick Start

```csharp
using CSharpNumerics.Physics.Quantum;
using CSharpNumerics.Engines.Quantum;

// Build a 2-qubit circuit that creates a Bell state
var circuit = new QuantumCircuit(2);
circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));
circuit.AddInstruction(new QuantumInstruction(new CNOTGate(), new List<int> { 0, 1 }));

// Simulate
var state = new QuantumSimulator().Run(circuit);

// Inspect probabilities
double p00 = state.GetProbability(0);   // 0.5  — |00⟩
double p11 = state.GetProbability(3);   // 0.5  — |11⟩

VectorN probs = state.GetProbabilities(); // [0.5, 0, 0, 0.5]
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

### QuantumState

Immutable result of a simulation — a complex amplitude vector with convenience methods.

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

---

### QuantumSimulator

Stateless simulator: initialises |0…0⟩ and applies each instruction in sequence.

```csharp
var sim = new QuantumSimulator();
QuantumState result = sim.Run(circuit);
```

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

### Module Structure

```
Engines/Quantum/
├── QuantumCircuit.cs        Circuit definition (qubits + instruction list)
├── QuantumInstruction.cs    Gate + target qubit binding
├── QuantumSimulator.cs      Stateless circuit executor
├── QuantumState.cs          Result: amplitudes + probability queries
└── README.md
```
