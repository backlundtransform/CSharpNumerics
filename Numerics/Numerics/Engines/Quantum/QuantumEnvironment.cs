using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Quantum;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Quantum;

/// <summary>
/// Reinforcement learning environment for quantum circuit synthesis.
/// The agent learns a sequence of quantum gates to prepare a target quantum state
/// starting from |0…0⟩. Each action applies one gate to the circuit.
///
/// <b>Observation:</b> current state probabilities + normalised step progress
/// <b>Actions:</b> discrete — each index maps to a (gate, qubit(s)) pair
/// <b>Reward:</b> fidelity improvement (delta) toward the target state
/// <b>Done:</b> max gates reached or fidelity ≥ threshold
/// </summary>
public class QuantumEnvironment : IEnvironment
{
    private readonly int _qubitCount;
    private readonly QuantumState _targetState;
    private readonly int _maxGates;
    private readonly double _fidelityThreshold;
    private readonly List<QuantumInstruction> _actionMap;

    private QuantumCircuit _circuit;
    private QuantumState _currentState;
    private double _currentFidelity;
    private int _step;
    private Random _rng;

    /// <summary>
    /// Probabilities (2^n) + step progress (1).
    /// </summary>
    public int ObservationSize { get; }

    /// <summary>Number of discrete gate actions available.</summary>
    public int ActionSize => _actionMap.Count;

    public bool IsDiscrete => true;

    private QuantumEnvironment(
        int qubitCount,
        QuantumState targetState,
        List<QuantumInstruction> actionMap,
        int maxGates,
        double fidelityThreshold)
    {
        _qubitCount = qubitCount;
        _targetState = targetState;
        _actionMap = actionMap;
        _maxGates = maxGates;
        _fidelityThreshold = fidelityThreshold;
        ObservationSize = (1 << qubitCount) + 1;
        _rng = new Random();
    }

    public (VectorN state, Dictionary<string, object> info) Reset(int? seed = null)
    {
        _rng = seed.HasValue ? new Random(seed.Value) : new Random();
        _circuit = new QuantumCircuit(_qubitCount);
        _currentState = new QuantumSimulator().Run(_circuit);
        _currentFidelity = QuantumFidelity.Fidelity(_currentState, _targetState);
        _step = 0;
        return (GetObservation(), new Dictionary<string, object>());
    }

    public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(int action)
    {
        if (action < 0 || action >= _actionMap.Count)
            throw new ArgumentOutOfRangeException(nameof(action));

        // Apply the selected gate to the circuit
        var instruction = _actionMap[action];
        _circuit.AddInstruction(instruction);

        // Re-simulate from scratch (statevector simulation)
        _currentState = new QuantumSimulator().Run(_circuit);
        _step++;

        double newFidelity = QuantumFidelity.Fidelity(_currentState, _targetState);
        double reward = newFidelity - _currentFidelity;
        _currentFidelity = newFidelity;

        bool done = _step >= _maxGates || _currentFidelity >= _fidelityThreshold;

        // Bonus reward for reaching the target
        if (_currentFidelity >= _fidelityThreshold)
            reward += 1.0;

        var info = new Dictionary<string, object>
        {
            { "fidelity", _currentFidelity },
            { "gates", _step }
        };

        return (GetObservation(), reward, done, info);
    }

    public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(VectorN action)
        => Step((int)action[0]);

    private VectorN GetObservation()
    {
        var probs = _currentState.GetProbabilities();
        int stateSize = 1 << _qubitCount;
        var obs = new double[stateSize + 1];
        for (int i = 0; i < stateSize; i++)
            obs[i] = probs[i];
        obs[stateSize] = (double)_step / _maxGates; // normalised progress
        return new VectorN(obs);
    }

    // ── Builder ────────────────────────────────────────────────

    /// <summary>
    /// Creates a builder for a quantum RL environment with the specified qubit count.
    /// </summary>
    public static QuantumEnvironmentBuilder Create(int qubitCount)
        => new QuantumEnvironmentBuilder(qubitCount);

    /// <summary>
    /// Fluent builder for configuring a <see cref="QuantumEnvironment"/>.
    /// </summary>
    public class QuantumEnvironmentBuilder
    {
        private readonly int _qubitCount;
        private QuantumState _targetState;
        private int _maxGates = 20;
        private double _fidelityThreshold = 0.99;
        private List<QuantumInstruction> _actionMap;

        internal QuantumEnvironmentBuilder(int qubitCount)
        {
            if (qubitCount <= 0)
                throw new ArgumentException("Qubit count must be positive.", nameof(qubitCount));
            _qubitCount = qubitCount;
        }

        /// <summary>Sets the target state the agent should learn to prepare.</summary>
        public QuantumEnvironmentBuilder WithTargetState(QuantumState target)
        {
            _targetState = target ?? throw new ArgumentNullException(nameof(target));
            return this;
        }

        /// <summary>
        /// Sets the target state by providing a circuit that produces it.
        /// The circuit is simulated once to obtain the target amplitudes.
        /// </summary>
        public QuantumEnvironmentBuilder WithTargetCircuit(QuantumCircuit circuit)
        {
            if (circuit == null) throw new ArgumentNullException(nameof(circuit));
            _targetState = new QuantumSimulator().Run(circuit);
            return this;
        }

        /// <summary>Maximum number of gates the agent can place per episode.</summary>
        public QuantumEnvironmentBuilder WithMaxGates(int maxGates)
        {
            if (maxGates <= 0)
                throw new ArgumentException("Max gates must be positive.", nameof(maxGates));
            _maxGates = maxGates;
            return this;
        }

        /// <summary>Fidelity threshold for early termination (default 0.99).</summary>
        public QuantumEnvironmentBuilder WithFidelityThreshold(double threshold)
        {
            if (threshold < 0 || threshold > 1)
                throw new ArgumentOutOfRangeException(nameof(threshold));
            _fidelityThreshold = threshold;
            return this;
        }

        /// <summary>
        /// Provides a custom action map. Each entry maps an action index to a gate instruction.
        /// If not set, a default action set is generated for single-qubit gates on all qubits
        /// plus CNOT on all qubit pairs.
        /// </summary>
        public QuantumEnvironmentBuilder WithActions(List<QuantumInstruction> actions)
        {
            _actionMap = actions ?? throw new ArgumentNullException(nameof(actions));
            return this;
        }

        /// <summary>Builds the environment.</summary>
        public QuantumEnvironment Build()
        {
            if (_targetState == null)
                throw new InvalidOperationException(
                    "Target state is required. Call WithTargetState() or WithTargetCircuit().");

            var actions = _actionMap ?? BuildDefaultActions();

            return new QuantumEnvironment(
                _qubitCount,
                _targetState,
                actions,
                _maxGates,
                _fidelityThreshold);
        }

        private List<QuantumInstruction> BuildDefaultActions()
        {
            var actions = new List<QuantumInstruction>();

            // Single-qubit gates on each qubit: H, X, Z, S, T
            QuantumGate[] singleGates = new QuantumGate[]
            {
                new HadamardGate(),
                new PauliXGate(),
                new PauliZGate(),
                new SGate(),
                new TGate()
            };

            for (int q = 0; q < _qubitCount; q++)
            {
                foreach (var gate in singleGates)
                    actions.Add(new QuantumInstruction(gate, new List<int> { q }));
            }

            // CNOT on all ordered qubit pairs
            for (int c = 0; c < _qubitCount; c++)
                for (int t = 0; t < _qubitCount; t++)
                    if (c != t)
                        actions.Add(new QuantumInstruction(new CNOTGate(), new List<int> { c, t }));

            return actions;
        }
    }
}
