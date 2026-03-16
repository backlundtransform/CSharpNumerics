# Reinforcement Learning — Roadmap

> Extends the CSharpNumerics ML framework with **Reinforcement Learning** (RL).
> Follows the same fluent, pipeline-based API patterns established by `SupervisedExperiment` and `ClusteringExperiment`.

---

## Feasibility Summary

| Question | Answer |
|----------|--------|
| Can RL use the existing fluent API style? | **Yes.** A `ReinforcementLearningExperiment` builder follows the identical `.For().With*().Run()` chain. |
| Can the existing MLP be reused? | **Yes.** `MLPClassifier`/`MLPRegressor` forward pass, backprop, and weight-update logic transfer directly to policy and value networks. A thin `NeuralNetwork` wrapper (extract from MLP) enables both actor and critic nets without copy-paste. |
| What about Matrix/VectorN? | Already sufficient for Q-tables, state/action tensors, and batch experience storage. |
| Scalers? | `StandardScaler` / `MinMaxScaler` map directly to **state normalisation** — critical for stable RL training. |
| Hyperparameter search? | `PipelineGrid` cartesian expansion works as-is; RL-specific hyperparams (γ, ε-schedule, τ) slot into the existing `Add("Gamma", 0.95, 0.99)` pattern. |
| Cross-validation? | RL doesn't use train/test splits. Instead, evaluation is **episodic**: average return over N evaluation episodes. `MonteCarloCrossValidator` concept maps to bootstrap confidence intervals over episode returns. |

**Conclusion:** the library already has ~80 % of the infrastructure. The main new work is the RL *loop* (environment ↔ agent ↔ replay buffer) and the RL-specific algorithms.

---

## Architecture Overview

```
ML/ReinforcementLearning/
├── Interfaces/
│   ├── IEnvironment.cs          // Environment contract
│   ├── IAgent.cs                // Agent contract
│   ├── IPolicy.cs               // Action-selection policy
│   └── IReplayBuffer.cs         // Experience storage
├── Core/
│   ├── State.cs                 // VectorN wrapper + metadata
│   ├── Action.cs                // Discrete / Continuous action
│   ├── Transition.cs            // (s, a, r, s', done)
│   ├── Episode.cs               // List<Transition> + total return
│   └── TrainingResult.cs        // Return curves, loss, diagnostics
├── Environments/
│   ├── DiscreteEnvironment.cs   // Base for discrete-action envs
│   ├── ContinuousEnvironment.cs // Base for continuous-action envs
│   ├── GridWorld.cs             // Classic grid navigation
│   ├── CartPole.cs              // Pole balancing (continuous state, discrete action)
│   ├── MountainCar.cs           // Energy-based control
│   └── Pendulum.cs              // Continuous control (torque)
├── Policies/
│   ├── EpsilonGreedy.cs         // ε-greedy with decay schedule
│   ├── Softmax.cs               // Boltzmann exploration
│   ├── GaussianNoise.cs         // Continuous exploration noise
│   └── OrnsteinUhlenbeck.cs     // Correlated noise for DDPG
├── Buffers/
│   ├── ReplayBuffer.cs          // Uniform random replay
│   └── PrioritizedReplayBuffer.cs // Proportional priority
├── Algorithms/
│   ├── Tabular/
│   │   ├── QLearning.cs         // Off-policy TD(0)
│   │   ├── SARSA.cs             // On-policy TD(0)
│   │   └── MonteCarloControl.cs // First-visit MC
│   ├── ValueBased/
│   │   ├── DQN.cs               // Deep Q-Network
│   │   ├── DoubleDQN.cs         // Double DQN (bias reduction)
│   │   └── DuelingDQN.cs        // Dueling architecture
│   ├── PolicyGradient/
│   │   ├── REINFORCE.cs         // Vanilla policy gradient
│   │   ├── ActorCritic.cs       // A2C (advantage actor-critic)
│   │   └── PPO.cs               // Proximal Policy Optimization
│   └── ContinuousControl/
│       └── DDPG.cs              // Deep Deterministic Policy Gradient
├── Experiment/
│   ├── RLExperiment.cs          // Fluent builder (entry point)
│   ├── RLExperimentResult.cs    // Results container
│   ├── RLPipelineGrid.cs        // Hyperparameter grid for RL
│   └── EpisodeEvaluator.cs      // Evaluation over N episodes
└── README.md
```

---

## Interfaces

### IEnvironment

```csharp
public interface IEnvironment
{
    int ObservationSize { get; }      // dim(state)
    int ActionSize { get; }           // number of discrete actions, or dim(continuous action)
    bool IsDiscrete { get; }

    (VectorN state, Dictionary<string, object> info) Reset(int? seed = null);
    (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(int action);        // discrete
    (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(VectorN action);    // continuous
}
```

### IAgent

```csharp
public interface IAgent
{
    string Name { get; }

    int SelectAction(VectorN state);              // discrete
    VectorN SelectContinuousAction(VectorN state); // continuous

    void Train(Transition transition);             // online single-step update
    void TrainBatch(List<Transition> batch);       // batch/mini-batch update
    void EndEpisode(Episode episode);              // post-episode hook (e.g. MC returns)

    IAgent Clone();
    Dictionary<string, object> GetHyperParameters();
    void SetHyperParameters(Dictionary<string, object> parameters);
}
```

### IPolicy

```csharp
public interface IPolicy
{
    int SelectAction(VectorN qValues);             // discrete: Q-values → action index
    VectorN SelectAction(VectorN mean, VectorN std); // continuous: mean+std → sampled action
    void Decay();                                   // per-episode decay (e.g. ε decay)
    IPolicy Clone();
}
```

### IReplayBuffer

```csharp
public interface IReplayBuffer
{
    void Add(Transition transition);
    List<Transition> Sample(int batchSize);
    int Count { get; }
    int Capacity { get; }
}
```

---

## Fluent API Design

Matches the existing `SupervisedExperiment` / `ClusteringExperiment` patterns exactly:

```csharp
// ── Tabular Q-Learning ──────────────────────────────────────
var result = RLExperiment
    .For(new GridWorld(rows: 5, cols: 5))
    .WithAgent(new QLearning
    {
        LearningRate = 0.1,
        Gamma = 0.99
    })
    .WithPolicy(new EpsilonGreedy
    {
        Epsilon = 1.0,
        EpsilonMin = 0.01,
        EpsilonDecay = 0.995
    })
    .WithEpisodes(maxEpisodes: 1000, maxStepsPerEpisode: 200)
    .Run();

double avgReturn = result.AverageReturn;
List<Serie> curve = result.ReturnCurve;        // episode vs return
List<Serie> loss  = result.LossCurve;          // episode vs loss
```

```csharp
// ── Deep Q-Network ──────────────────────────────────────────
var result = RLExperiment
    .For(new CartPole())
    .WithAgent(new DQN
    {
        HiddenLayers = new[] { 64, 64 },
        Activation = ActivationType.ReLU,
        LearningRate = 0.001,
        Gamma = 0.99,
        TargetUpdateFrequency = 100,
        BatchSize = 32
    })
    .WithPolicy(new EpsilonGreedy
    {
        Epsilon = 1.0,
        EpsilonMin = 0.05,
        EpsilonDecay = 0.999
    })
    .WithReplayBuffer(capacity: 10000)
    .WithEpisodes(maxEpisodes: 500, maxStepsPerEpisode: 500)
    .WithEvaluation(evalEpisodes: 20, evalInterval: 50)
    .Run();
```

```csharp
// ── Policy Gradient (REINFORCE) ─────────────────────────────
var result = RLExperiment
    .For(new CartPole())
    .WithAgent(new REINFORCE
    {
        HiddenLayers = new[] { 32 },
        LearningRate = 0.01,
        Gamma = 0.99
    })
    .WithEpisodes(maxEpisodes: 1000)
    .Run();
```

```csharp
// ── Actor-Critic (A2C) ─────────────────────────────────────
var result = RLExperiment
    .For(new CartPole())
    .WithAgent(new ActorCritic
    {
        ActorHiddenLayers = new[] { 32 },
        CriticHiddenLayers = new[] { 32 },
        ActorLearningRate = 0.001,
        CriticLearningRate = 0.002,
        Gamma = 0.99
    })
    .WithEpisodes(maxEpisodes: 1000)
    .Run();
```

```csharp
// ── Continuous Control (DDPG) ───────────────────────────────
var result = RLExperiment
    .For(new Pendulum())
    .WithAgent(new DDPG
    {
        ActorHiddenLayers = new[] { 64, 64 },
        CriticHiddenLayers = new[] { 64, 64 },
        ActorLearningRate = 1e-4,
        CriticLearningRate = 1e-3,
        Gamma = 0.99,
        Tau = 0.005,
        BatchSize = 64
    })
    .WithPolicy(new OrnsteinUhlenbeck
    {
        Sigma = 0.2,
        Theta = 0.15
    })
    .WithReplayBuffer(capacity: 100000)
    .WithEpisodes(maxEpisodes: 500, maxStepsPerEpisode: 200)
    .Run();
```

```csharp
// ── Hyperparameter Grid Search ──────────────────────────────
var result = RLExperiment
    .For(new CartPole())
    .WithGrid(new RLPipelineGrid()
        .AddAgent<DQN>(g => g
            .Add("LearningRate", 0.001, 0.0005)
            .Add("Gamma", 0.95, 0.99)
            .Add("HiddenLayers", new[] { 32 }, new[] { 64, 64 }))
        .AddAgent<REINFORCE>(g => g
            .Add("LearningRate", 0.01, 0.005)))
    .WithPolicy(new EpsilonGreedy { Epsilon = 1.0, EpsilonDecay = 0.995 })
    .WithReplayBuffer(capacity: 10000)
    .WithEpisodes(maxEpisodes: 300)
    .WithEvaluation(evalEpisodes: 20)
    .Run();

// Rankings sorted by average evaluation return
foreach (var r in result.Rankings)
    Console.WriteLine($"{r.AgentName} [{r.Parameters}]: {r.AverageReturn:F1} ± {r.StdDev:F1}");
```

```csharp
// ── Monte Carlo Evaluation (confidence intervals) ───────────
var result = RLExperiment
    .For(new CartPole())
    .WithAgent(new DQN { /* ... */ })
    .WithPolicy(new EpsilonGreedy { /* ... */ })
    .WithReplayBuffer(capacity: 10000)
    .WithEpisodes(maxEpisodes: 500)
    .WithMonteCarloEvaluation(runs: 50, evalEpisodes: 30, seed: 42)
    .Run();

double mean = result.MonteCarlo.MeanReturn;
var (lo, hi) = result.MonteCarlo.ConfidenceInterval(0.95);
```

---

## Phased Implementation Plan

### Phase 0 — Neural Network Refactor (prerequisite)

> Extract the forward/backward pass from `MLPClassifier`/`MLPRegressor` into a shared `NeuralNetwork` class that all three paradigms (supervised, unsupervised-autoencoders, RL) can reuse.

| Task | Description |
|------|-------------|
| 0.1 | Create `ML/NeuralNetwork/NeuralNetwork.cs` — forward pass, backprop, weight update, activation functions, Xavier init |
| 0.2 | Add `Predict(VectorN input) → VectorN` (single-sample forward, no batch) |
| 0.3 | Add `Backward(VectorN loss gradient) → gradient update` (for RL policy gradient) |
| 0.4 | Add `CopyWeightsFrom(NeuralNetwork source)` (for DQN target-network sync) |
| 0.5 | Add `SoftUpdate(NeuralNetwork source, double tau)` (for DDPG/SAC Polyak averaging) |
| 0.6 | Refactor `MLPClassifier` and `MLPRegressor` to delegate to `NeuralNetwork` internally |
| 0.7 | Verify all existing ML tests still pass |

**Impact:** zero breaking changes — `MLPClassifier` and `MLPRegressor` keep their public API.

---

### Phase 1 — Core Abstractions & Tabular RL

> Establish the foundational interfaces, data types, and the simplest RL algorithms.

| Task | Description |
|------|-------------|
| 1.1 | `Interfaces/` — `IEnvironment`, `IAgent`, `IPolicy`, `IReplayBuffer` |
| 1.2 | `Core/` — `State`, `Action`, `Transition`, `Episode`, `TrainingResult` |
| 1.3 | `Policies/EpsilonGreedy.cs` — ε-greedy with linear/exponential decay |
| 1.4 | `Environments/GridWorld.cs` — simple grid navigation (discrete states/actions) |
| 1.5 | `Algorithms/Tabular/QLearning.cs` — off-policy TD(0) with Q-table as `Matrix` |
| 1.6 | `Algorithms/Tabular/SARSA.cs` — on-policy TD(0) |
| 1.7 | `Algorithms/Tabular/MonteCarloControl.cs` — first-visit MC with ε-greedy |
| 1.8 | `Experiment/RLExperiment.cs` — fluent builder (`.For().WithAgent().WithPolicy().WithEpisodes().Run()`) |
| 1.9 | `Experiment/RLExperimentResult.cs` — return curves, loss curves, diagnostics |
| 1.10 | Unit tests: GridWorld convergence for Q-Learning, SARSA, MC |

**Deliverable:** fully working tabular RL with the fluent API.

---

### Phase 2 — Deep RL (Value-Based)

> Introduce function approximation using the `NeuralNetwork` from Phase 0.

| Task | Description |
|------|-------------|
| 2.1 | `Buffers/ReplayBuffer.cs` — uniform random sampling, circular buffer |
| 2.2 | `Environments/CartPole.cs` — continuous state, discrete action |
| 2.3 | `Environments/MountainCar.cs` — sparse reward, energy-based |
| 2.4 | `Algorithms/ValueBased/DQN.cs` — online + target network, Huber loss |
| 2.5 | `Algorithms/ValueBased/DoubleDQN.cs` — decoupled action selection / evaluation |
| 2.6 | `Algorithms/ValueBased/DuelingDQN.cs` — V + A stream architecture |
| 2.7 | `Buffers/PrioritizedReplayBuffer.cs` — proportional priority via sum-tree |
| 2.8 | Extend `RLExperiment` with `.WithReplayBuffer()` and `.WithEvaluation()` |
| 2.9 | Unit tests: CartPole solved (avg return ≥ 195) with DQN |

**Deliverable:** deep Q-learning family with replay buffers and target networks.

---

### Phase 3 — Policy Gradient Methods

> Direct policy optimisation — no Q-values, learn π(a|s) directly.

| Task | Description |
|------|-------------|
| 3.1 | `Algorithms/PolicyGradient/REINFORCE.cs` — Monte Carlo policy gradient with baseline |
| 3.2 | `Algorithms/PolicyGradient/ActorCritic.cs` — A2C (advantage actor-critic) with separate value baseline |
| 3.3 | `Algorithms/PolicyGradient/PPO.cs` — clipped surrogate objective, mini-batch updates |
| 3.4 | `Policies/Softmax.cs` — Boltzmann exploration for discrete policy gradient |
| 3.5 | Unit tests: CartPole with REINFORCE, A2C, PPO |

**Deliverable:** policy gradient suite — the most versatile RL methods.

---

### Phase 4 — Continuous Control

> Extend to continuous action spaces.

| Task | Description |
|------|-------------|
| 4.1 | `Environments/Pendulum.cs` — continuous torque control |
| 4.2 | `Policies/GaussianNoise.cs` — additive Gaussian exploration |
| 4.3 | `Policies/OrnsteinUhlenbeck.cs` — temporally correlated noise |
| 4.4 | `Algorithms/ContinuousControl/DDPG.cs` — deterministic policy gradient + soft target updates |
| 4.5 | Extend `IAgent` continuous action path |
| 4.6 | Unit tests: Pendulum swing-up with DDPG |

**Deliverable:** continuous control with DDPG.

---

### Phase 5 — Experiment Grid & Evaluation ✅

> Full hyperparameter search and statistical evaluation, mirroring `SupervisedExperiment`.

| Task | Description | Status |
|------|-------------|--------|
| 5.1 | `Experiment/RLPipelineGrid.cs` — cartesian expansion over agent hyperparameters | ✅ |
| 5.2 | `Experiment/EpisodeEvaluator.cs` — greedy evaluation with summary statistics (mean, std, CI, median) | ✅ |
| 5.3 | `Experiment/RLGridSearchResult.cs` — grid search result, entry, and Monte Carlo result types | ✅ |
| 5.4 | Extend `RLExperiment` with `.WithGrid()`, `.WithMonteCarloEvaluation()`, `.RunGrid()`, `.RunMonteCarlo()` | ✅ |
| 5.5 | Unit tests: 16 tests covering grid expansion, evaluator, grid search, Monte Carlo evaluation | ✅ |

**Deliverable:** production-quality experiment framework matching the supervised/clustering API. **112 total RL tests pass.**

---

### Phase 6 — Diagnostics & Visualisation Hooks ✅

> Tools for understanding agent behaviour.

| Task | Description | Status |
|------|-------------|--------|
| 6.1 | `TrainingResult.ReturnCurve` — episode return as `List<Serie>` | ✅ (Phase 1) |
| 6.2 | `TrainingResult.LossCurve` — training loss as `List<Serie>` | ✅ (Phase 1) |
| 6.3 | `TrainingResult.ExplorationCurve` — ε or noise σ over time | ✅ (Phase 1) |
| 6.4 | `Diagnostics/QValueHeatmap.cs` — Q-value heatmap, greedy policy, Q-table matrix | ✅ |
| 6.5 | `Diagnostics/PolicyVisualizer.cs` — action probabilities, softmax probs, entropy, dominant action | ✅ |
| 6.6 | `Diagnostics/ValueFunctionSurface.cs` — 1D/2D value surface, DQN/A2C/PPO extractors, ToMatrix | ✅ |
| 6.7 | `PPO.GetValue(VectorN)` — public value accessor for diagnostics | ✅ |
| 6.8 | Unit tests: 16 tests covering all diagnostics | ✅ |

**Deliverable:** `List<Serie>`-based diagnostics ready for the existing export pipeline. **128 total RL tests pass.**

---

## Mapping to Existing Infrastructure

| Existing Component | RL Usage |
|--------------------|----------|
| `Matrix` | Q-table storage, batch transition matrices, weight matrices |
| `VectorN` | State vectors, action vectors, reward vectors |
| `MLPClassifier` backprop | Policy network gradient computation |
| `MLPRegressor` forward pass | Value / Q-network forward evaluation |
| `ActivationType` enum | Shared across policy and value networks |
| `StandardScaler` | State normalisation (running mean/std) |
| `IHasHyperparameters` | Agent hyperparameter injection for grid search |
| `PipelineGrid` pattern | `RLPipelineGrid` follows identical cartesian expansion |
| `MonteCarloCrossValidator` concept | Bootstrap confidence intervals over episode returns |
| `Serie` / `List<Serie>` | Return curves, loss curves, exploration schedules |
| `MonteCarloSimulator` | Stochastic evaluation of trained policies |

---

## Algorithms Summary

| Algorithm | Type | State | Action | Network(s) | Replay Buffer | Key Idea |
|-----------|------|-------|--------|------------|--------------|----------|
| Q-Learning | Tabular | Discrete | Discrete | — | — | Off-policy TD: $Q(s,a) \leftarrow Q + \alpha[r + \gamma \max_{a'} Q(s',a') - Q(s,a)]$ |
| SARSA | Tabular | Discrete | Discrete | — | — | On-policy TD: uses $Q(s',a')$ where $a' \sim \pi$ |
| MC Control | Tabular | Discrete | Discrete | — | — | Averages full episodic returns |
| DQN | Deep | Continuous | Discrete | Q-net + target | Yes | Experience replay + target network stabilisation |
| Double DQN | Deep | Continuous | Discrete | Q-net + target | Yes | Decouples action selection from evaluation |
| Dueling DQN | Deep | Continuous | Discrete | V + A streams | Yes | $Q(s,a) = V(s) + A(s,a) - \text{mean}(A)$ |
| REINFORCE | PG | Any | Discrete | Policy net | — | $\nabla J \propto \sum_t \nabla \log \pi(a_t|s_t) \cdot G_t$ |
| A2C | PG | Any | Discrete | Actor + Critic | — | Advantage reduces variance: $A(s,a) = Q(s,a) - V(s)$ |
| PPO | PG | Any | Discrete | Actor + Critic | — | Clipped ratio: $\min(r_t A_t, \text{clip}(r_t) A_t)$ |
| DDPG | Continuous | Continuous | Continuous | Actor + Critic × 2 | Yes | Deterministic policy + soft target updates |

---

## Why This Fits CSharpNumerics

1. **Educational transparency.** Every algorithm is implemented from first principles on `Matrix`/`VectorN` — no opaque tensor libraries. Users can step through Q-updates, policy gradients, and replay sampling.

2. **Fluent consistency.** `RLExperiment.For(env).WithAgent(agent).Run()` is immediately familiar to anyone who has used `SupervisedExperiment.For(X,y).WithGrid(grid).Run()`.

3. **Composable.** Mix and match environments, agents, policies, buffers, and evaluators — the same pipeline philosophy.

4. **Numerical.** RL is fundamentally numerical: matrix multiplications for Q-networks, gradient descent for policy updates, stochastic sampling for exploration. Everything maps to existing primitives.

---

## Estimated Scale

| Phase | New files | Lines of code (est.) |
|-------|-----------|---------------------|
| 0 — NN refactor | 1 new + 2 refactored | ~300 |
| 1 — Core + Tabular | ~12 | ~1200 |
| 2 — Deep RL (DQN) | ~8 | ~1000 |
| 3 — Policy Gradient | ~5 | ~800 |
| 4 — Continuous | ~5 | ~700 |
| 5 — Experiment Grid | ~4 | ~600 |
| 6 — Diagnostics | ~3 | ~400 |
| **Total** | **~38 files** | **~5000 lines** |

Comparable in scope to the Physics/Oscillations + Waves namespaces combined.
