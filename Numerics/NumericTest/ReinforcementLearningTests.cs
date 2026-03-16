using CSharpNumerics.ML.ReinforcementLearning.Algorithms.ContinuousControl;
using CSharpNumerics.ML.ReinforcementLearning.Algorithms.PolicyGradient;
using CSharpNumerics.ML.ReinforcementLearning.Algorithms.Tabular;
using CSharpNumerics.ML.ReinforcementLearning.Algorithms.ValueBased;
using CSharpNumerics.ML.ReinforcementLearning.Buffers;
using CSharpNumerics.ML.ReinforcementLearning.Core;
using CSharpNumerics.ML.ReinforcementLearning.Diagnostics;
using CSharpNumerics.ML.ReinforcementLearning.Environments;
using CSharpNumerics.ML.ReinforcementLearning.Experiment;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.ML.ReinforcementLearning.Policies;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Linq;

namespace NumericTest;

[TestClass]
public class ReinforcementLearningTests
{
    // ── GridWorld Environment Tests ─────────────────────────────

    [TestMethod]
    public void GridWorld_Reset_ShouldReturnOriginState()
    {
        var env = new GridWorld(5, 5);
        var (state, _) = env.Reset(42);
        Assert.AreEqual(0.0, state[0]); // row = 0
        Assert.AreEqual(0.0, state[1]); // col = 0
    }

    [TestMethod]
    public void GridWorld_StepRight_ShouldMoveAgent()
    {
        var env = new GridWorld(5, 5);
        env.Reset(42);
        var (state, reward, done, _) = env.Step(1); // Right
        Assert.AreEqual(0.0, state[0]);
        Assert.AreEqual(1.0, state[1]);
        Assert.IsFalse(done);
        Assert.AreEqual(-0.01, reward, 1e-10);
    }

    [TestMethod]
    public void GridWorld_ReachGoal_ShouldReturnDone()
    {
        var env = new GridWorld(2, 2); // goal at (1,1)
        env.Reset(42);
        env.Step(1); // (0,1)
        var (state, reward, done, _) = env.Step(2); // (1,1) = goal
        Assert.AreEqual(1.0, state[0]);
        Assert.AreEqual(1.0, state[1]);
        Assert.IsTrue(done);
        Assert.AreEqual(1.0, reward);
    }

    [TestMethod]
    public void GridWorld_WallBlocks_ShouldStayInPlace()
    {
        var walls = new HashSet<(int, int)> { (0, 1) };
        var env = new GridWorld(3, 3, walls);
        env.Reset(42);
        var (state, _, _, _) = env.Step(1); // try to go right into wall
        Assert.AreEqual(0.0, state[0]); // stayed at (0,0)
        Assert.AreEqual(0.0, state[1]);
    }

    [TestMethod]
    public void GridWorld_BoundaryBlocks_ShouldStayInPlace()
    {
        var env = new GridWorld(3, 3);
        env.Reset(42);
        var (state, _, _, _) = env.Step(0); // try to go up from (0,0)
        Assert.AreEqual(0.0, state[0]);
        Assert.AreEqual(0.0, state[1]);
    }

    [TestMethod]
    public void GridWorld_StateToIndex_ShouldBeCorrect()
    {
        var env = new GridWorld(5, 5);
        var state = new VectorN(new double[] { 2, 3 });
        Assert.AreEqual(13, env.StateToIndex(state)); // 2*5 + 3
    }

    [TestMethod]
    public void GridWorld_Properties_ShouldBeCorrect()
    {
        var env = new GridWorld(4, 6);
        Assert.AreEqual(2, env.ObservationSize);
        Assert.AreEqual(4, env.ActionSize);
        Assert.IsTrue(env.IsDiscrete);
        Assert.AreEqual(24, env.StateCount);
    }

    // ── EpsilonGreedy Policy Tests ──────────────────────────────

    [TestMethod]
    public void EpsilonGreedy_FullExploration_ShouldBeRandom()
    {
        var policy = new EpsilonGreedy(seed: 42) { Epsilon = 1.0 };
        var qValues = new VectorN(new double[] { 10, 1, 1, 1 });

        // With ε=1, not all actions should be 0 (greedy)
        int nonGreedy = 0;
        for (int i = 0; i < 100; i++)
            if (policy.SelectAction(qValues) != 0) nonGreedy++;

        Assert.IsTrue(nonGreedy > 10, "Should explore frequently when ε=1");
    }

    [TestMethod]
    public void EpsilonGreedy_NoExploration_ShouldBeGreedy()
    {
        var policy = new EpsilonGreedy(seed: 42) { Epsilon = 0.0 };
        var qValues = new VectorN(new double[] { 1, 5, 2, 3 });

        for (int i = 0; i < 50; i++)
            Assert.AreEqual(1, policy.SelectAction(qValues), "Should always pick argmax when ε=0");
    }

    [TestMethod]
    public void EpsilonGreedy_Decay_ShouldReduceEpsilon()
    {
        var policy = new EpsilonGreedy { Epsilon = 1.0, EpsilonMin = 0.01, EpsilonDecay = 0.9 };
        policy.Decay();
        Assert.AreEqual(0.9, policy.Epsilon, 1e-10);
        policy.Decay();
        Assert.AreEqual(0.81, policy.Epsilon, 1e-10);
    }

    [TestMethod]
    public void EpsilonGreedy_Decay_ShouldRespectMinimum()
    {
        var policy = new EpsilonGreedy { Epsilon = 0.02, EpsilonMin = 0.01, EpsilonDecay = 0.1 };
        policy.Decay();
        Assert.AreEqual(0.01, policy.Epsilon, 1e-10); // clamped to min
    }

    // ── Core Data Types ─────────────────────────────────────────

    [TestMethod]
    public void Episode_DiscountedReturns_ShouldComputeCorrectly()
    {
        var ep = new Episode();
        var s = new VectorN(new double[] { 0, 0 });
        ep.Transitions.Add(new Transition(s, 0, 1.0, s, false));
        ep.Transitions.Add(new Transition(s, 0, 2.0, s, false));
        ep.Transitions.Add(new Transition(s, 0, 3.0, s, true));

        double gamma = 0.9;
        var returns = ep.DiscountedReturns(gamma);

        // G_2 = 3
        // G_1 = 2 + 0.9*3 = 4.7
        // G_0 = 1 + 0.9*4.7 = 5.23
        Assert.AreEqual(3.0, returns[2], 1e-10);
        Assert.AreEqual(4.7, returns[1], 1e-10);
        Assert.AreEqual(5.23, returns[0], 1e-10);
    }

    [TestMethod]
    public void Episode_TotalReturn_ShouldBeSumOfRewards()
    {
        var ep = new Episode();
        var s = new VectorN(new double[] { 0, 0 });
        ep.Transitions.Add(new Transition(s, 0, 1.0, s, false));
        ep.Transitions.Add(new Transition(s, 0, -0.5, s, false));
        ep.Transitions.Add(new Transition(s, 0, 2.0, s, true));

        Assert.AreEqual(2.5, ep.TotalReturn, 1e-10);
        Assert.AreEqual(3, ep.Length);
    }

    // ── Q-Learning Tests ────────────────────────────────────────

    [TestMethod]
    public void QLearning_SingleUpdate_ShouldModifyQValue()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, s => (int)s[0] * 3 + (int)s[1])
        {
            LearningRate = 0.1,
            Gamma = 0.9
        };

        var s = new VectorN(new double[] { 0, 0 });
        var sNext = new VectorN(new double[] { 0, 1 });

        // Initial Q should be 0
        var q0 = agent.GetQValues(s);
        Assert.AreEqual(0.0, q0[1]); // action 1 (right)

        // After one update: Q(0,0, right) += 0.1 * (-0.01 + 0.9*0 - 0) = -0.001
        agent.Train(new Transition(s, 1, -0.01, sNext, false));
        var q1 = agent.GetQValues(s);
        Assert.AreEqual(-0.001, q1[1], 1e-10);
    }

    [TestMethod]
    public void QLearning_GridWorld_ShouldConverge()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex)
        {
            LearningRate = 0.1,
            Gamma = 0.99
        };

        var result = RLExperiment
            .For(env)
            .WithAgent(agent)
            .WithPolicy(new EpsilonGreedy(seed: 42)
            {
                Epsilon = 1.0,
                EpsilonMin = 0.01,
                EpsilonDecay = 0.995
            })
            .WithEpisodes(maxEpisodes: 2000, maxStepsPerEpisode: 100)
            .WithSeed(42)
            .Run();

        // Should learn to reach the goal
        double avgLast100 = result.AverageReturnLastN(100);
        Assert.IsTrue(avgLast100 > 0.5,
            $"Q-Learning should converge on 3×3 grid, but avg last 100 = {avgLast100:F3}");
        Assert.IsTrue(result.TotalEpisodes == 2000);
        Assert.IsTrue(result.ReturnCurve.Count == 2000);
    }

    // ── SARSA Tests ─────────────────────────────────────────────

    [TestMethod]
    public void SARSA_GridWorld_ShouldConverge()
    {
        var env = new GridWorld(3, 3);
        var agent = new SARSA(9, 4, env.StateToIndex)
        {
            LearningRate = 0.1,
            Gamma = 0.99
        };

        var result = RLExperiment
            .For(env)
            .WithAgent(agent)
            .WithPolicy(new EpsilonGreedy(seed: 42)
            {
                Epsilon = 1.0,
                EpsilonMin = 0.01,
                EpsilonDecay = 0.995
            })
            .WithEpisodes(maxEpisodes: 2000, maxStepsPerEpisode: 100)
            .WithSeed(42)
            .Run();

        double avgLast100 = result.AverageReturnLastN(100);
        Assert.IsTrue(avgLast100 > 0.5,
            $"SARSA should converge on 3×3 grid, but avg last 100 = {avgLast100:F3}");
    }

    // ── Monte Carlo Control Tests ───────────────────────────────

    [TestMethod]
    public void MonteCarloControl_GridWorld_ShouldConverge()
    {
        var env = new GridWorld(3, 3);
        var agent = new MonteCarloControl(9, 4, env.StateToIndex)
        {
            LearningRate = 0.01,
            Gamma = 0.99
        };

        var result = RLExperiment
            .For(env)
            .WithAgent(agent)
            .WithPolicy(new EpsilonGreedy(seed: 42)
            {
                Epsilon = 1.0,
                EpsilonMin = 0.05,
                EpsilonDecay = 0.999
            })
            .WithEpisodes(maxEpisodes: 5000, maxStepsPerEpisode: 100)
            .WithSeed(42)
            .Run();

        double avgLast100 = result.AverageReturnLastN(100);
        Assert.IsTrue(avgLast100 > 0.3,
            $"MC Control should converge on 3×3 grid, but avg last 100 = {avgLast100:F3}");
    }

    // ── RLExperiment Fluent API Tests ───────────────────────────

    [TestMethod]
    public void RLExperiment_WithoutAgent_ShouldThrow()
    {
        var env = new GridWorld(3, 3);
        Assert.ThrowsException<InvalidOperationException>(() =>
            RLExperiment.For(env).WithEpisodes(10).Run());
    }

    [TestMethod]
    public void RLExperiment_DefaultPolicy_ShouldWork()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex);

        // No explicit policy — should default to EpsilonGreedy
        var result = RLExperiment
            .For(env)
            .WithAgent(agent)
            .WithEpisodes(maxEpisodes: 10, maxStepsPerEpisode: 50)
            .Run();

        Assert.AreEqual(10, result.TotalEpisodes);
        Assert.AreEqual("QLearning", result.AgentName);
    }

    [TestMethod]
    public void RLExperiment_Result_ShouldContainDiagnostics()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex);

        var result = RLExperiment
            .For(env)
            .WithAgent(agent)
            .WithPolicy(new EpsilonGreedy(seed: 42))
            .WithEpisodes(maxEpisodes: 50, maxStepsPerEpisode: 30)
            .WithSeed(0)
            .Run();

        Assert.AreEqual(50, result.TotalEpisodes);
        Assert.AreEqual(50, result.ReturnCurve.Count);
        Assert.IsTrue(result.ExplorationCurve.Count > 0);
        Assert.IsTrue(result.Duration.TotalMilliseconds >= 0);
        Assert.IsNotNull(result.Parameters);
        Assert.IsTrue(result.Parameters.ContainsKey("LearningRate"));
        Assert.IsTrue(result.Parameters.ContainsKey("Gamma"));
    }

    [TestMethod]
    public void RLExperiment_ExplorationDecays_ShouldBeDecreasing()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex);

        var result = RLExperiment
            .For(env)
            .WithAgent(agent)
            .WithPolicy(new EpsilonGreedy(seed: 42)
            {
                Epsilon = 1.0,
                EpsilonDecay = 0.9
            })
            .WithEpisodes(maxEpisodes: 20, maxStepsPerEpisode: 10)
            .Run();

        var curve = result.ExplorationCurve;
        Assert.IsTrue(curve.Count == 20);
        Assert.IsTrue(curve[0].Value > curve[^1].Value, "Exploration should decrease over time");
    }

    // ── TabularAgent Tests ──────────────────────────────────────

    [TestMethod]
    public void TabularAgent_Clone_ShouldCopyQTable()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex) { LearningRate = 0.5 };

        var s = new VectorN(new double[] { 1, 1 });
        var sNext = new VectorN(new double[] { 1, 2 });
        agent.Train(new Transition(s, 0, 1.0, sNext, true));

        var clone = (QLearning)agent.Clone();
        var qOriginal = agent.GetQValues(s);
        var qClone = clone.GetQValues(s);

        Assert.AreEqual(qOriginal[0], qClone[0], 1e-10);
    }

    [TestMethod]
    public void TabularAgent_GetQTable_ShouldReturnMatrix()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex);

        var qTable = agent.GetQTable();
        Assert.AreEqual(9, qTable.rowLength);
        Assert.AreEqual(4, qTable.columnLength);
    }

    // ── TrainingResult Tests ────────────────────────────────────

    [TestMethod]
    public void TrainingResult_AverageReturnLastN_ShouldWork()
    {
        var tr = new TrainingResult();
        tr.EpisodeReturns.AddRange(new[] { 1.0, 2.0, 3.0, 4.0, 5.0 });

        Assert.AreEqual(3.0, tr.AverageReturn, 1e-10);
        Assert.AreEqual(4.5, tr.AverageReturnLastN(2), 1e-10);
        Assert.AreEqual(5.0, tr.BestReturn, 1e-10);
        Assert.AreEqual(5, tr.TotalEpisodes);
    }

    // ════════════════════════════════════════════════════════════
    //  Phase 2 — Deep RL (Value-Based) Tests
    // ════════════════════════════════════════════════════════════

    // ── CartPole Environment Tests ──────────────────────────────

    [TestMethod]
    public void CartPole_Properties_ShouldBeCorrect()
    {
        var env = new CartPole();
        Assert.AreEqual(4, env.ObservationSize);
        Assert.AreEqual(2, env.ActionSize);
        Assert.IsTrue(env.IsDiscrete);
    }

    [TestMethod]
    public void CartPole_Reset_ShouldReturnSmallState()
    {
        var env = new CartPole();
        var (state, _) = env.Reset(42);
        Assert.AreEqual(4, state.Length);
        for (int i = 0; i < 4; i++)
            Assert.IsTrue(Math.Abs(state[i]) < 0.1, $"state[{i}] should be near zero after reset");
    }

    [TestMethod]
    public void CartPole_Step_ShouldChangeState()
    {
        var env = new CartPole();
        var (s0, _) = env.Reset(42);
        var (s1, reward, done, _) = env.Step(1);
        Assert.AreEqual(4, s1.Length);
        Assert.AreEqual(1.0, reward);
        // State should have changed
        bool changed = false;
        for (int i = 0; i < 4; i++)
            if (Math.Abs(s1[i] - s0[i]) > 1e-10) changed = true;
        Assert.IsTrue(changed, "State should change after a step");
    }

    [TestMethod]
    public void CartPole_ExtremeAngle_ShouldTerminate()
    {
        var env = new CartPole();
        env.Reset(42);
        // Push right repeatedly — eventually pole falls
        bool terminated = false;
        for (int i = 0; i < 500; i++)
        {
            var (_, _, done, _) = env.Step(1);
            if (done) { terminated = true; break; }
        }
        Assert.IsTrue(terminated, "CartPole should terminate when pole falls");
    }

    // ── MountainCar Environment Tests ───────────────────────────

    [TestMethod]
    public void MountainCar_Properties_ShouldBeCorrect()
    {
        var env = new MountainCar();
        Assert.AreEqual(2, env.ObservationSize);
        Assert.AreEqual(3, env.ActionSize);
        Assert.IsTrue(env.IsDiscrete);
    }

    [TestMethod]
    public void MountainCar_Reset_ShouldReturnValidState()
    {
        var env = new MountainCar();
        var (state, _) = env.Reset(42);
        Assert.AreEqual(2, state.Length);
        Assert.IsTrue(state[0] >= -0.6 && state[0] <= -0.4, "Position should be in [-0.6, -0.4]");
        Assert.AreEqual(0.0, state[1], 1e-10, "Velocity should be 0");
    }

    [TestMethod]
    public void MountainCar_Step_NegativeReward()
    {
        var env = new MountainCar();
        env.Reset(42);
        var (_, reward, _, _) = env.Step(1); // no push
        Assert.AreEqual(-1.0, reward);
    }

    [TestMethod]
    public void MountainCar_MaxSteps_ShouldTerminate()
    {
        var env = new MountainCar { MaxSteps = 50 };
        env.Reset(42);
        bool terminated = false;
        for (int i = 0; i < 60; i++)
        {
            var (_, _, done, _) = env.Step(1); // do nothing
            if (done) { terminated = true; break; }
        }
        Assert.IsTrue(terminated, "MountainCar should terminate at MaxSteps");
    }

    // ── ReplayBuffer Tests ──────────────────────────────────────

    [TestMethod]
    public void ReplayBuffer_AddAndCount_ShouldWork()
    {
        var buf = new ReplayBuffer(100, seed: 42);
        Assert.AreEqual(0, buf.Count);
        Assert.AreEqual(100, buf.Capacity);

        var s = new VectorN(new[] { 1.0, 2.0 });
        buf.Add(new Transition(s, 0, 1.0, s, false));
        Assert.AreEqual(1, buf.Count);
    }

    [TestMethod]
    public void ReplayBuffer_Circular_ShouldOverwrite()
    {
        var buf = new ReplayBuffer(3, seed: 42);
        var s = new VectorN(new[] { 0.0 });

        for (int i = 0; i < 5; i++)
            buf.Add(new Transition(new VectorN(new[] { (double)i }), 0, 0, s, false));

        Assert.AreEqual(3, buf.Count); // capped at capacity
    }

    [TestMethod]
    public void ReplayBuffer_Sample_ShouldReturnCorrectSize()
    {
        var buf = new ReplayBuffer(100, seed: 42);
        var s = new VectorN(new[] { 1.0 });
        for (int i = 0; i < 20; i++)
            buf.Add(new Transition(s, 0, 1.0, s, false));

        var batch = buf.Sample(10);
        Assert.AreEqual(10, batch.Count);
    }

    [TestMethod]
    public void ReplayBuffer_SampleTooMany_ShouldThrow()
    {
        var buf = new ReplayBuffer(10, seed: 42);
        var s = new VectorN(new[] { 1.0 });
        buf.Add(new Transition(s, 0, 1.0, s, false));
        Assert.ThrowsException<InvalidOperationException>(() => buf.Sample(5));
    }

    // ── PrioritizedReplayBuffer Tests ───────────────────────────

    [TestMethod]
    public void PrioritizedReplayBuffer_AddAndSample_ShouldWork()
    {
        var buf = new PrioritizedReplayBuffer(100, seed: 42);
        var s = new VectorN(new[] { 1.0, 2.0 });
        for (int i = 0; i < 10; i++)
            buf.Add(new Transition(s, 0, 1.0, s, false));

        Assert.AreEqual(10, buf.Count);
        var batch = buf.Sample(5);
        Assert.AreEqual(5, batch.Count);
    }

    [TestMethod]
    public void PrioritizedReplayBuffer_UpdatePriority_ShouldAffectSampling()
    {
        var buf = new PrioritizedReplayBuffer(10, seed: 42);
        var low = new VectorN(new[] { 0.0 });
        var high = new VectorN(new[] { 1.0 });

        buf.Add(new Transition(low, 0, 0, low, false));   // index 0
        buf.Add(new Transition(high, 1, 1, high, false));  // index 1

        // Give index 1 much higher priority
        buf.UpdatePriority(0, 0.001);
        buf.UpdatePriority(1, 100.0);

        // Sample many times — index 1 (action=1) should dominate
        int highCount = 0;
        for (int i = 0; i < 100; i++)
        {
            var batch = buf.Sample(1);
            if (batch[0].Action == 1) highCount++;
        }
        Assert.IsTrue(highCount > 80, $"High-priority transition should be sampled most ({highCount}/100)");
    }

    [TestMethod]
    public void PrioritizedReplayBuffer_SampleWithIndices_ShouldReturnIndices()
    {
        var buf = new PrioritizedReplayBuffer(50, seed: 42);
        var s = new VectorN(new[] { 1.0 });
        for (int i = 0; i < 20; i++)
            buf.Add(new Transition(s, i % 3, 1.0, s, false));

        var (transitions, indices) = buf.SampleWithIndices(8);
        Assert.AreEqual(8, transitions.Count);
        Assert.AreEqual(8, indices.Count);
        foreach (var idx in indices)
            Assert.IsTrue(idx >= 0 && idx < 20);
    }

    // ── DQN Tests ───────────────────────────────────────────────

    [TestMethod]
    public void DQN_Initialize_ShouldCreateNetworks()
    {
        var dqn = new DQN
        {
            HiddenLayers = new[] { 16 },
            LearningRate = 0.01
        };
        dqn.Initialize(4, 2, seed: 42);

        // Should be able to get Q-values after initialization
        var state = new VectorN(new[] { 0.1, 0.2, 0.3, 0.4 });
        var q = dqn.GetQValues(state);
        Assert.AreEqual(2, q.Length);
    }

    [TestMethod]
    public void DQN_WithoutInitialize_ShouldThrow()
    {
        var dqn = new DQN();
        var state = new VectorN(new[] { 0.1, 0.2 });
        Assert.ThrowsException<InvalidOperationException>(() => dqn.SelectAction(state));
    }

    [TestMethod]
    public void DQN_Train_ShouldNotCrash()
    {
        var dqn = new DQN
        {
            HiddenLayers = new[] { 16 },
            BatchSize = 4,
            MinBufferSize = 4
        };
        dqn.Initialize(4, 2, seed: 42);

        var rng = new Random(42);
        for (int i = 0; i < 10; i++)
        {
            var s = new VectorN(new[] { rng.NextDouble(), rng.NextDouble(), rng.NextDouble(), rng.NextDouble() });
            var sn = new VectorN(new[] { rng.NextDouble(), rng.NextDouble(), rng.NextDouble(), rng.NextDouble() });
            dqn.Train(new Transition(s, rng.Next(2), rng.NextDouble(), sn, false));
        }
        // Should complete without exception
    }

    [TestMethod]
    public void DQN_Clone_ShouldCopyHyperparameters()
    {
        var dqn = new DQN
        {
            HiddenLayers = new[] { 32, 32 },
            LearningRate = 0.005,
            Gamma = 0.95
        };

        var clone = (DQN)dqn.Clone();
        Assert.AreEqual(0.005, clone.LearningRate);
        Assert.AreEqual(0.95, clone.Gamma);
        Assert.AreEqual(2, clone.HiddenLayers.Length);
    }

    // ── DoubleDQN Tests ─────────────────────────────────────────

    [TestMethod]
    public void DoubleDQN_Initialize_ShouldWork()
    {
        var agent = new DoubleDQN
        {
            HiddenLayers = new[] { 16 },
            LearningRate = 0.01
        };
        agent.Initialize(4, 2, seed: 42);

        var state = new VectorN(new[] { 0.1, 0.2, 0.3, 0.4 });
        var q = agent.GetQValues(state);
        Assert.AreEqual(2, q.Length);
    }

    [TestMethod]
    public void DoubleDQN_Train_ShouldNotCrash()
    {
        var agent = new DoubleDQN
        {
            HiddenLayers = new[] { 16 },
            BatchSize = 4,
            MinBufferSize = 4
        };
        agent.Initialize(4, 2, seed: 42);

        var rng = new Random(42);
        for (int i = 0; i < 10; i++)
        {
            var s = new VectorN(new[] { rng.NextDouble(), rng.NextDouble(), rng.NextDouble(), rng.NextDouble() });
            var sn = new VectorN(new[] { rng.NextDouble(), rng.NextDouble(), rng.NextDouble(), rng.NextDouble() });
            agent.Train(new Transition(s, rng.Next(2), rng.NextDouble(), sn, false));
        }
    }

    // ── DuelingDQN Tests ────────────────────────────────────────

    [TestMethod]
    public void DuelingDQN_Initialize_ShouldWork()
    {
        var agent = new DuelingDQN
        {
            SharedLayers = new[] { 16 },
            ValueLayers = new[] { 8 },
            AdvantageLayers = new[] { 8 }
        };
        agent.Initialize(4, 2, seed: 42);

        var state = new VectorN(new[] { 0.1, 0.2, 0.3, 0.4 });
        var q = agent.GetQValues(state);
        Assert.AreEqual(2, q.Length);
    }

    [TestMethod]
    public void DuelingDQN_QDecomposition_ShouldSatisfyIdentity()
    {
        // Q(s,a) = V(s) + A(s,a) - mean(A) means Σ_a A(s,a) should center around 0
        var agent = new DuelingDQN
        {
            SharedLayers = new[] { 16 },
            ValueLayers = new[] { 8 },
            AdvantageLayers = new[] { 8 }
        };
        agent.Initialize(4, 3, seed: 42);

        var state = new VectorN(new[] { 0.5, 0.5, 0.5, 0.5 });
        var q = agent.GetQValues(state);
        Assert.AreEqual(3, q.Length);
        // All Q-values should be finite
        for (int i = 0; i < 3; i++)
            Assert.IsFalse(double.IsNaN(q[i]), $"Q[{i}] should not be NaN");
    }

    [TestMethod]
    public void DuelingDQN_Train_ShouldNotCrash()
    {
        var agent = new DuelingDQN
        {
            SharedLayers = new[] { 16 },
            ValueLayers = new[] { 8 },
            AdvantageLayers = new[] { 8 },
            BatchSize = 4,
            MinBufferSize = 4
        };
        agent.Initialize(4, 2, seed: 42);

        var rng = new Random(42);
        for (int i = 0; i < 10; i++)
        {
            var s = new VectorN(new[] { rng.NextDouble(), rng.NextDouble(), rng.NextDouble(), rng.NextDouble() });
            var sn = new VectorN(new[] { rng.NextDouble(), rng.NextDouble(), rng.NextDouble(), rng.NextDouble() });
            agent.Train(new Transition(s, rng.Next(2), rng.NextDouble(), sn, false));
        }
    }

    // ── RLExperiment with Deep Agents ───────────────────────────

    [TestMethod]
    public void RLExperiment_DQN_CartPole_ShouldRun()
    {
        var result = RLExperiment
            .For(new CartPole())
            .WithAgent(new DQN
            {
                HiddenLayers = new[] { 32 },
                LearningRate = 0.001,
                BatchSize = 16,
                MinBufferSize = 32,
                TargetUpdateFrequency = 50
            })
            .WithPolicy(new EpsilonGreedy(seed: 42)
            {
                Epsilon = 1.0,
                EpsilonMin = 0.05,
                EpsilonDecay = 0.995
            })
            .WithReplayBuffer(5000, seed: 42)
            .WithEpisodes(maxEpisodes: 50, maxStepsPerEpisode: 200)
            .WithSeed(42)
            .Run();

        Assert.AreEqual("DQN", result.AgentName);
        Assert.AreEqual(50, result.TotalEpisodes);
        Assert.IsTrue(result.Training.EpisodeReturns.Count == 50);
    }

    [TestMethod]
    public void RLExperiment_DoubleDQN_CartPole_ShouldRun()
    {
        var result = RLExperiment
            .For(new CartPole())
            .WithAgent(new DoubleDQN
            {
                HiddenLayers = new[] { 32 },
                LearningRate = 0.001,
                BatchSize = 16,
                MinBufferSize = 32
            })
            .WithPolicy(new EpsilonGreedy(seed: 42)
            {
                Epsilon = 1.0,
                EpsilonMin = 0.05,
                EpsilonDecay = 0.995
            })
            .WithReplayBuffer(5000, seed: 42)
            .WithEpisodes(maxEpisodes: 30, maxStepsPerEpisode: 200)
            .WithSeed(42)
            .Run();

        Assert.AreEqual(30, result.TotalEpisodes);
    }

    [TestMethod]
    public void RLExperiment_DuelingDQN_CartPole_ShouldRun()
    {
        var result = RLExperiment
            .For(new CartPole())
            .WithAgent(new DuelingDQN
            {
                SharedLayers = new[] { 32 },
                ValueLayers = new[] { 16 },
                AdvantageLayers = new[] { 16 },
                LearningRate = 0.001,
                BatchSize = 16,
                MinBufferSize = 32
            })
            .WithPolicy(new EpsilonGreedy(seed: 42)
            {
                Epsilon = 1.0,
                EpsilonMin = 0.05,
                EpsilonDecay = 0.995
            })
            .WithReplayBuffer(5000, seed: 42)
            .WithEpisodes(maxEpisodes: 30, maxStepsPerEpisode: 200)
            .WithSeed(42)
            .Run();

        Assert.AreEqual(30, result.TotalEpisodes);
    }

    [TestMethod]
    public void RLExperiment_WithEvaluation_ShouldRecordEvalReturns()
    {
        var result = RLExperiment
            .For(new CartPole())
            .WithAgent(new DQN
            {
                HiddenLayers = new[] { 16 },
                LearningRate = 0.001,
                BatchSize = 8,
                MinBufferSize = 16
            })
            .WithPolicy(new EpsilonGreedy(seed: 42)
            {
                Epsilon = 0.5,
                EpsilonDecay = 0.99
            })
            .WithReplayBuffer(2000, seed: 42)
            .WithEpisodes(maxEpisodes: 100, maxStepsPerEpisode: 100)
            .WithEvaluation(evalEpisodes: 5, evalInterval: 50)
            .WithSeed(42)
            .Run();

        // 100 episodes / 50 interval = 2 evaluation checkpoints
        Assert.AreEqual(2, result.Training.EvalReturns.Count);
    }

    [TestMethod]
    public void RLExperiment_WithExternalReplayBuffer_ShouldUseIt()
    {
        var buffer = new PrioritizedReplayBuffer(3000, seed: 42);

        var result = RLExperiment
            .For(new CartPole())
            .WithAgent(new DQN
            {
                HiddenLayers = new[] { 16 },
                BatchSize = 8,
                MinBufferSize = 16
            })
            .WithPolicy(new EpsilonGreedy(seed: 42))
            .WithReplayBuffer(buffer)
            .WithEpisodes(maxEpisodes: 20, maxStepsPerEpisode: 50)
            .WithSeed(42)
            .Run();

        Assert.AreEqual(20, result.TotalEpisodes);
        Assert.IsTrue(buffer.Count > 0, "External buffer should have been filled");
    }

    // ════════════════════════════════════════════════════════════
    //  Phase 3 — Policy Gradient Methods Tests
    // ════════════════════════════════════════════════════════════

    // ── SoftmaxPolicy Tests ─────────────────────────────────────

    [TestMethod]
    public void SoftmaxPolicy_HighTemp_ShouldBeUniform()
    {
        var policy = new SoftmaxPolicy(seed: 42) { Temperature = 100.0 };
        var qValues = new VectorN(new double[] { 10.0, 1.0, 1.0, 1.0 });

        var counts = new int[4];
        for (int i = 0; i < 1000; i++)
            counts[policy.SelectAction(qValues)]++;

        // At high temperature all actions should be roughly equally likely
        foreach (var c in counts)
            Assert.IsTrue(c > 100, $"High temperature should spread actions, count={c}");
    }

    [TestMethod]
    public void SoftmaxPolicy_LowTemp_ShouldBeGreedy()
    {
        var policy = new SoftmaxPolicy(seed: 42) { Temperature = 0.001 };
        var qValues = new VectorN(new double[] { 1.0, 5.0, 2.0, 3.0 });

        int greedyCount = 0;
        for (int i = 0; i < 100; i++)
            if (policy.SelectAction(qValues) == 1) greedyCount++;

        Assert.IsTrue(greedyCount > 95, $"Low temperature should be near-greedy, got {greedyCount}/100");
    }

    [TestMethod]
    public void SoftmaxPolicy_Decay_ShouldReduceTemperature()
    {
        var policy = new SoftmaxPolicy { Temperature = 1.0, TemperatureMin = 0.01, TemperatureDecay = 0.9 };
        policy.Decay();
        Assert.AreEqual(0.9, policy.Temperature, 1e-10);
        for (int i = 0; i < 100; i++) policy.Decay();
        Assert.IsTrue(policy.Temperature >= 0.01, "Temperature should respect minimum");
    }

    [TestMethod]
    public void SoftmaxPolicy_Clone_ShouldCopy()
    {
        var policy = new SoftmaxPolicy { Temperature = 0.5, TemperatureMin = 0.1, TemperatureDecay = 0.99 };
        var clone = (SoftmaxPolicy)policy.Clone();
        Assert.AreEqual(0.5, clone.Temperature, 1e-10);
        Assert.AreEqual(0.1, clone.TemperatureMin, 1e-10);
    }

    // ── REINFORCE Tests ─────────────────────────────────────────

    [TestMethod]
    public void REINFORCE_Initialize_ShouldCreateNetwork()
    {
        var agent = new REINFORCE
        {
            HiddenLayers = new[] { 16 },
            LearningRate = 0.01
        };
        agent.Initialize(4, 2, seed: 42);

        var state = new VectorN(new double[] { 0.1, 0.2, 0.3, 0.4 });
        var probs = agent.GetActionProbabilities(state);
        Assert.AreEqual(2, probs.Length);
        // Softmax output sums to 1
        double sum = 0;
        for (int i = 0; i < probs.Length; i++) sum += probs[i];
        Assert.AreEqual(1.0, sum, 1e-6);
    }

    [TestMethod]
    public void REINFORCE_WithoutInitialize_ShouldThrow()
    {
        var agent = new REINFORCE();
        var state = new VectorN(new double[] { 0.1, 0.2 });
        Assert.ThrowsException<InvalidOperationException>(() => agent.SelectAction(state));
    }

    [TestMethod]
    public void REINFORCE_SelectAction_ShouldReturnValidAction()
    {
        var agent = new REINFORCE { HiddenLayers = new[] { 16 } };
        agent.Initialize(4, 3, seed: 42);

        var state = new VectorN(new double[] { 0.1, 0.2, 0.3, 0.4 });
        int action = agent.SelectAction(state);
        Assert.IsTrue(action >= 0 && action < 3);
    }

    [TestMethod]
    public void REINFORCE_Clone_ShouldCopyHyperparameters()
    {
        var agent = new REINFORCE
        {
            HiddenLayers = new[] { 32 },
            LearningRate = 0.005,
            Gamma = 0.95,
            UseBaseline = false
        };

        var clone = (REINFORCE)agent.Clone();
        Assert.AreEqual(0.005, clone.LearningRate);
        Assert.AreEqual(0.95, clone.Gamma);
        Assert.IsFalse(clone.UseBaseline);
    }

    [TestMethod]
    public void REINFORCE_CartPole_ShouldRunViaExperiment()
    {
        var result = RLExperiment
            .For(new CartPole())
            .WithAgent(new REINFORCE
            {
                HiddenLayers = new[] { 16 },
                LearningRate = 0.01,
                Gamma = 0.99
            })
            .WithEpisodes(maxEpisodes: 30, maxStepsPerEpisode: 200)
            .WithSeed(42)
            .Run();

        Assert.AreEqual("REINFORCE", result.AgentName);
        Assert.AreEqual(30, result.TotalEpisodes);
        Assert.IsTrue(result.Parameters.ContainsKey("LearningRate"));
        Assert.IsTrue(result.Parameters.ContainsKey("Gamma"));
    }

    // ── ActorCritic (A2C) Tests ─────────────────────────────────

    [TestMethod]
    public void ActorCritic_Initialize_ShouldCreateNetworks()
    {
        var agent = new ActorCritic
        {
            ActorHiddenLayers = new[] { 16 },
            CriticHiddenLayers = new[] { 16 }
        };
        agent.Initialize(4, 2, seed: 42);

        var state = new VectorN(new double[] { 0.1, 0.2, 0.3, 0.4 });

        var probs = agent.GetActionProbabilities(state);
        Assert.AreEqual(2, probs.Length);
        double sum = 0;
        for (int i = 0; i < probs.Length; i++) sum += probs[i];
        Assert.AreEqual(1.0, sum, 1e-6);

        double value = agent.GetValue(state);
        Assert.IsFalse(double.IsNaN(value));
    }

    [TestMethod]
    public void ActorCritic_WithoutInitialize_ShouldThrow()
    {
        var agent = new ActorCritic();
        var state = new VectorN(new double[] { 0.1, 0.2 });
        Assert.ThrowsException<InvalidOperationException>(() => agent.SelectAction(state));
    }

    [TestMethod]
    public void ActorCritic_Train_ShouldNotCrash()
    {
        var agent = new ActorCritic
        {
            ActorHiddenLayers = new[] { 16 },
            CriticHiddenLayers = new[] { 16 }
        };
        agent.Initialize(4, 2, seed: 42);

        var rng = new Random(42);
        for (int i = 0; i < 10; i++)
        {
            var s = new VectorN(new[] { rng.NextDouble(), rng.NextDouble(), rng.NextDouble(), rng.NextDouble() });
            var sn = new VectorN(new[] { rng.NextDouble(), rng.NextDouble(), rng.NextDouble(), rng.NextDouble() });
            agent.Train(new Transition(s, rng.Next(2), rng.NextDouble(), sn, false));
        }
    }

    [TestMethod]
    public void ActorCritic_CartPole_ShouldRunViaExperiment()
    {
        var result = RLExperiment
            .For(new CartPole())
            .WithAgent(new ActorCritic
            {
                ActorHiddenLayers = new[] { 16 },
                CriticHiddenLayers = new[] { 16 },
                ActorLearningRate = 0.001,
                CriticLearningRate = 0.002,
                Gamma = 0.99
            })
            .WithEpisodes(maxEpisodes: 30, maxStepsPerEpisode: 200)
            .WithSeed(42)
            .Run();

        Assert.AreEqual("A2C", result.AgentName);
        Assert.AreEqual(30, result.TotalEpisodes);
        Assert.IsTrue(result.Parameters.ContainsKey("ActorLearningRate"));
        Assert.IsTrue(result.Parameters.ContainsKey("CriticLearningRate"));
    }

    [TestMethod]
    public void ActorCritic_Clone_ShouldCopyHyperparameters()
    {
        var agent = new ActorCritic
        {
            ActorHiddenLayers = new[] { 64 },
            CriticHiddenLayers = new[] { 32 },
            ActorLearningRate = 0.005,
            CriticLearningRate = 0.01,
            Gamma = 0.95,
            EntropyCoefficient = 0.02
        };

        var clone = (ActorCritic)agent.Clone();
        Assert.AreEqual(0.005, clone.ActorLearningRate);
        Assert.AreEqual(0.01, clone.CriticLearningRate);
        Assert.AreEqual(0.95, clone.Gamma);
        Assert.AreEqual(0.02, clone.EntropyCoefficient);
    }

    // ── PPO Tests ───────────────────────────────────────────────

    [TestMethod]
    public void PPO_Initialize_ShouldCreateNetworks()
    {
        var agent = new PPO
        {
            ActorHiddenLayers = new[] { 16 },
            CriticHiddenLayers = new[] { 16 }
        };
        agent.Initialize(4, 2, seed: 42);

        var state = new VectorN(new double[] { 0.1, 0.2, 0.3, 0.4 });
        var probs = agent.GetActionProbabilities(state);
        Assert.AreEqual(2, probs.Length);
        double sum = 0;
        for (int i = 0; i < probs.Length; i++) sum += probs[i];
        Assert.AreEqual(1.0, sum, 1e-6);
    }

    [TestMethod]
    public void PPO_WithoutInitialize_ShouldThrow()
    {
        var agent = new PPO();
        var state = new VectorN(new double[] { 0.1, 0.2 });
        Assert.ThrowsException<InvalidOperationException>(() => agent.SelectAction(state));
    }

    [TestMethod]
    public void PPO_SelectAction_ShouldReturnValidAction()
    {
        var agent = new PPO
        {
            ActorHiddenLayers = new[] { 16 },
            CriticHiddenLayers = new[] { 16 }
        };
        agent.Initialize(4, 3, seed: 42);

        var state = new VectorN(new double[] { 0.1, 0.2, 0.3, 0.4 });
        int action = agent.SelectAction(state);
        Assert.IsTrue(action >= 0 && action < 3);
    }

    [TestMethod]
    public void PPO_CartPole_ShouldRunViaExperiment()
    {
        var result = RLExperiment
            .For(new CartPole())
            .WithAgent(new PPO
            {
                ActorHiddenLayers = new[] { 16 },
                CriticHiddenLayers = new[] { 16 },
                ActorLearningRate = 0.001,
                CriticLearningRate = 0.002,
                Gamma = 0.99,
                ClipEpsilon = 0.2,
                UpdateEpochs = 2,
                MiniBatchSize = 32
            })
            .WithEpisodes(maxEpisodes: 30, maxStepsPerEpisode: 200)
            .WithSeed(42)
            .Run();

        Assert.AreEqual("PPO", result.AgentName);
        Assert.AreEqual(30, result.TotalEpisodes);
        Assert.IsTrue(result.Parameters.ContainsKey("ClipEpsilon"));
        Assert.IsTrue(result.Parameters.ContainsKey("Lambda"));
    }

    [TestMethod]
    public void PPO_Clone_ShouldCopyHyperparameters()
    {
        var agent = new PPO
        {
            ActorHiddenLayers = new[] { 64, 64 },
            CriticHiddenLayers = new[] { 64, 64 },
            ClipEpsilon = 0.1,
            UpdateEpochs = 10,
            Lambda = 0.9
        };

        var clone = (PPO)agent.Clone();
        Assert.AreEqual(0.1, clone.ClipEpsilon);
        Assert.AreEqual(10, clone.UpdateEpochs);
        Assert.AreEqual(0.9, clone.Lambda);
        Assert.AreEqual(2, clone.ActorHiddenLayers.Length);
    }

    [TestMethod]
    public void PPO_EndEpisode_ShouldUpdateWeights()
    {
        var agent = new PPO
        {
            ActorHiddenLayers = new[] { 8 },
            CriticHiddenLayers = new[] { 8 },
            UpdateEpochs = 1,
            MiniBatchSize = 4
        };
        agent.Initialize(2, 2, seed: 42);

        var state = new VectorN(new double[] { 0.5, 0.5 });
        var probsBefore = agent.GetActionProbabilities(state);

        // Simulate an episode by calling SelectAction + Train + EndEpisode
        var episode = new Episode();
        var rng = new Random(42);
        for (int i = 0; i < 10; i++)
        {
            var s = new VectorN(new[] { rng.NextDouble(), rng.NextDouble() });
            int act = agent.SelectAction(s);
            var sn = new VectorN(new[] { rng.NextDouble(), rng.NextDouble() });
            var t = new Transition(s, act, 1.0, sn, i == 9);
            agent.Train(t);
            episode.Transitions.Add(t);
        }
        agent.EndEpisode(episode);

        var probsAfter = agent.GetActionProbabilities(state);
        // Probabilities should have changed after training
        bool changed = false;
        for (int i = 0; i < probsBefore.Length; i++)
            if (Math.Abs(probsBefore[i] - probsAfter[i]) > 1e-10) changed = true;
        Assert.IsTrue(changed, "PPO should update policy after EndEpisode");
    }

    // ── Policy Gradient with Softmax Policy + Q-Learning ────────

    [TestMethod]
    public void SoftmaxPolicy_WithQLearning_ShouldRun()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex);

        var result = RLExperiment
            .For(env)
            .WithAgent(agent)
            .WithPolicy(new SoftmaxPolicy(seed: 42)
            {
                Temperature = 1.0,
                TemperatureMin = 0.1,
                TemperatureDecay = 0.99
            })
            .WithEpisodes(maxEpisodes: 100, maxStepsPerEpisode: 50)
            .WithSeed(42)
            .Run();

        Assert.AreEqual(100, result.TotalEpisodes);
    }

    // ════════════════════════════════════════════════════════════
    //  Phase 4 — Continuous Control Tests
    // ════════════════════════════════════════════════════════════

    // ── Pendulum Environment Tests ──────────────────────────────

    [TestMethod]
    public void Pendulum_Properties_ShouldBeCorrect()
    {
        var env = new Pendulum();
        Assert.AreEqual(3, env.ObservationSize);
        Assert.AreEqual(1, env.ActionSize);
        Assert.IsFalse(env.IsDiscrete);
    }

    [TestMethod]
    public void Pendulum_Reset_ShouldReturnValidState()
    {
        var env = new Pendulum();
        var (state, _) = env.Reset(42);
        Assert.AreEqual(3, state.Length);
        // cos²(θ) + sin²(θ) = 1
        double cosTheta = state[0];
        double sinTheta = state[1];
        Assert.AreEqual(1.0, cosTheta * cosTheta + sinTheta * sinTheta, 1e-10);
    }

    [TestMethod]
    public void Pendulum_Step_ShouldChangeState()
    {
        var env = new Pendulum();
        var (s0, _) = env.Reset(42);
        var action = new VectorN(new[] { 1.0 });
        var (s1, reward, done, _) = env.Step(action);
        Assert.AreEqual(3, s1.Length);
        Assert.IsTrue(reward <= 0, "Reward should be non-positive");
        Assert.IsFalse(done, "First step should not be terminal");
    }

    [TestMethod]
    public void Pendulum_MaxSteps_ShouldTerminate()
    {
        var env = new Pendulum { MaxSteps = 10 };
        env.Reset(42);
        bool terminated = false;
        for (int i = 0; i < 15; i++)
        {
            var (_, _, done, _) = env.Step(new VectorN(new[] { 0.0 }));
            if (done) { terminated = true; break; }
        }
        Assert.IsTrue(terminated, "Pendulum should terminate at MaxSteps");
    }

    [TestMethod]
    public void Pendulum_ActionClamp_ShouldRespectBounds()
    {
        var env = new Pendulum { MaxTorque = 2.0 };
        env.Reset(42);
        // Apply extreme torque — should be clamped internally
        var (_, reward1, _, _) = env.Step(new VectorN(new[] { 100.0 }));
        Assert.IsTrue(reward1 <= 0);
    }

    [TestMethod]
    public void Pendulum_DiscreteStep_ShouldWork()
    {
        var env = new Pendulum();
        env.Reset(42);
        // Discrete step: 0 = -MaxTorque, 1 = 0, 2 = +MaxTorque
        var (state, reward, done, _) = env.Step(2);
        Assert.AreEqual(3, state.Length);
        Assert.IsTrue(reward <= 0);
    }

    // ── Transition Continuous Action Tests ───────────────────────

    [TestMethod]
    public void Transition_ContinuousAction_ShouldStoreAction()
    {
        var s = new VectorN(new[] { 1.0, 2.0 });
        var a = new VectorN(new[] { 0.5, -0.3 });
        var sn = new VectorN(new[] { 1.1, 2.1 });
        var t = new Transition(s, a, -1.0, sn, false);

        Assert.AreEqual(-1, t.Action); // discrete action is -1 for continuous
        Assert.IsTrue(t.IsContinuous);
        Assert.AreEqual(0.5, t.ContinuousAction[0], 1e-10);
        Assert.AreEqual(-0.3, t.ContinuousAction[1], 1e-10);
        Assert.AreEqual(-1.0, t.Reward);
    }

    [TestMethod]
    public void Transition_DiscreteAction_ShouldNotBeContinuous()
    {
        var s = new VectorN(new[] { 1.0 });
        var t = new Transition(s, 2, 1.0, s, true);
        Assert.AreEqual(2, t.Action);
        Assert.IsFalse(t.IsContinuous);
    }

    // ── GaussianNoise Policy Tests ──────────────────────────────

    [TestMethod]
    public void GaussianNoise_ShouldAddNoise()
    {
        var policy = new GaussianNoise(seed: 42) { Sigma = 0.5 };
        var mean = new VectorN(new[] { 0.0, 0.0 });
        var std = new VectorN(new[] { 1.0, 1.0 });

        int different = 0;
        for (int i = 0; i < 100; i++)
        {
            var noisy = policy.SelectAction(mean, std);
            if (Math.Abs(noisy[0]) > 0.01 || Math.Abs(noisy[1]) > 0.01)
                different++;
        }
        Assert.IsTrue(different > 80, $"Gaussian noise should perturb action, got {different}/100 different");
    }

    [TestMethod]
    public void GaussianNoise_Decay_ShouldReduceSigma()
    {
        var policy = new GaussianNoise { Sigma = 1.0, SigmaMin = 0.01, SigmaDecay = 0.9 };
        policy.Decay();
        Assert.AreEqual(0.9, policy.Sigma, 1e-10);
        for (int i = 0; i < 100; i++) policy.Decay();
        Assert.IsTrue(policy.Sigma >= 0.01);
    }

    [TestMethod]
    public void GaussianNoise_Clone_ShouldCopy()
    {
        var policy = new GaussianNoise { Sigma = 0.3, SigmaMin = 0.05 };
        var clone = (GaussianNoise)policy.Clone();
        Assert.AreEqual(0.3, clone.Sigma, 1e-10);
        Assert.AreEqual(0.05, clone.SigmaMin, 1e-10);
    }

    // ── OrnsteinUhlenbeck Policy Tests ──────────────────────────

    [TestMethod]
    public void OrnsteinUhlenbeck_ShouldAddCorrelatedNoise()
    {
        var ou = new OrnsteinUhlenbeck(seed: 42)
        {
            Theta = 0.15,
            Sigma = 0.3
        };
        ou.Reset(1);

        var mean = new VectorN(new[] { 0.0 });
        var std = new VectorN(new[] { 1.0 });

        var prev = ou.SelectAction(mean, std);
        int correlatedCount = 0;
        for (int i = 0; i < 50; i++)
        {
            var next = ou.SelectAction(mean, std);
            // OU noise should be correlated (same sign consecutive often)
            if (Math.Sign(prev[0]) == Math.Sign(next[0]) && Math.Abs(next[0]) > 0.001)
                correlatedCount++;
            prev = next;
        }
        // OU should show some temporal correlation
        Assert.IsTrue(correlatedCount > 15, $"OU noise should be correlated, got {correlatedCount}/50");
    }

    [TestMethod]
    public void OrnsteinUhlenbeck_Reset_ShouldClearState()
    {
        var ou = new OrnsteinUhlenbeck(seed: 42) { Sigma = 1.0 };
        ou.Reset(2);

        var mean = new VectorN(new[] { 0.0, 0.0 });
        var std = new VectorN(new[] { 1.0, 1.0 });

        // Take several steps to accumulate state
        for (int i = 0; i < 10; i++)
            ou.SelectAction(mean, std);

        // Reset and first value should be close to mean (noise from Mu=0)
        ou.Reset(2);
        var first = ou.SelectAction(mean, std);
        // After reset, noise state is 0 + single step of noise
        Assert.IsTrue(Math.Abs(first[0]) < 2.0, "After reset, noise should be small");
    }

    [TestMethod]
    public void OrnsteinUhlenbeck_Decay_ShouldReduceSigma()
    {
        var ou = new OrnsteinUhlenbeck { Sigma = 0.5, SigmaMin = 0.01, SigmaDecay = 0.9 };
        ou.Decay();
        Assert.AreEqual(0.45, ou.Sigma, 1e-10);
    }

    // ── DDPG Tests ──────────────────────────────────────────────

    [TestMethod]
    public void DDPG_Initialize_ShouldCreateNetworks()
    {
        var agent = new DDPG
        {
            ActorHiddenLayers = new[] { 16 },
            CriticHiddenLayers = new[] { 16 },
            ActionScale = 2.0
        };
        agent.Initialize(3, 1, seed: 42);

        var state = new VectorN(new[] { 0.5, 0.5, 0.1 });
        var action = agent.SelectContinuousAction(state);
        Assert.AreEqual(1, action.Length);
        Assert.IsTrue(Math.Abs(action[0]) <= 2.0 + 1e-10, "Action should be within [-scale, scale]");
    }

    [TestMethod]
    public void DDPG_WithoutInitialize_ShouldThrow()
    {
        var agent = new DDPG();
        var state = new VectorN(new[] { 0.1, 0.2, 0.3 });
        Assert.ThrowsException<InvalidOperationException>(() => agent.SelectContinuousAction(state));
    }

    [TestMethod]
    public void DDPG_ActionScale_ShouldClampOutput()
    {
        var agent = new DDPG
        {
            ActorHiddenLayers = new[] { 16 },
            CriticHiddenLayers = new[] { 16 },
            ActionScale = 2.0
        };
        agent.Initialize(3, 1, seed: 42);

        // Run many states — all actions should be within [-2, 2]
        var rng = new Random(42);
        for (int i = 0; i < 100; i++)
        {
            var s = new VectorN(new[] { rng.NextDouble() * 2 - 1, rng.NextDouble() * 2 - 1, rng.NextDouble() * 2 - 1 });
            var a = agent.SelectContinuousAction(s);
            Assert.IsTrue(Math.Abs(a[0]) <= 2.0 + 1e-10, $"Action {a[0]} exceeds scale");
        }
    }

    [TestMethod]
    public void DDPG_Train_ShouldNotCrash()
    {
        var agent = new DDPG
        {
            ActorHiddenLayers = new[] { 16 },
            CriticHiddenLayers = new[] { 16 },
            BatchSize = 4,
            MinBufferSize = 4,
            ActionScale = 2.0
        };
        agent.Initialize(3, 1, seed: 42);

        var rng = new Random(42);
        for (int i = 0; i < 10; i++)
        {
            var s = new VectorN(new[] { rng.NextDouble(), rng.NextDouble(), rng.NextDouble() });
            var a = new VectorN(new[] { rng.NextDouble() * 4 - 2 });
            var sn = new VectorN(new[] { rng.NextDouble(), rng.NextDouble(), rng.NextDouble() });
            agent.Train(new Transition(s, a, -rng.NextDouble(), sn, false));
        }
    }

    [TestMethod]
    public void DDPG_Clone_ShouldCopyHyperparameters()
    {
        var agent = new DDPG
        {
            ActorHiddenLayers = new[] { 64, 64 },
            CriticHiddenLayers = new[] { 64, 64 },
            ActorLearningRate = 5e-4,
            CriticLearningRate = 1e-3,
            Gamma = 0.95,
            Tau = 0.01,
            ActionScale = 2.0
        };

        var clone = (DDPG)agent.Clone();
        Assert.AreEqual(5e-4, clone.ActorLearningRate);
        Assert.AreEqual(1e-3, clone.CriticLearningRate);
        Assert.AreEqual(0.95, clone.Gamma);
        Assert.AreEqual(0.01, clone.Tau);
        Assert.AreEqual(2.0, clone.ActionScale);
    }

    // ── DDPG via RLExperiment ───────────────────────────────────

    [TestMethod]
    public void RLExperiment_DDPG_Pendulum_ShouldRun()
    {
        var result = RLExperiment
            .For(new Pendulum { MaxSteps = 50 })
            .WithAgent(new DDPG
            {
                ActorHiddenLayers = new[] { 16 },
                CriticHiddenLayers = new[] { 16 },
                ActorLearningRate = 1e-3,
                CriticLearningRate = 1e-3,
                Gamma = 0.99,
                Tau = 0.005,
                BatchSize = 16,
                MinBufferSize = 32,
                ActionScale = 2.0
            })
            .WithPolicy(new GaussianNoise(seed: 42) { Sigma = 0.2 })
            .WithReplayBuffer(5000, seed: 42)
            .WithEpisodes(maxEpisodes: 20, maxStepsPerEpisode: 50)
            .WithSeed(42)
            .Run();

        Assert.AreEqual("DDPG", result.AgentName);
        Assert.AreEqual(20, result.TotalEpisodes);
        Assert.IsTrue(result.Parameters.ContainsKey("Tau"));
        Assert.IsTrue(result.Parameters.ContainsKey("ActionScale"));
    }

    [TestMethod]
    public void RLExperiment_DDPG_WithOU_ShouldRun()
    {
        var result = RLExperiment
            .For(new Pendulum { MaxSteps = 50 })
            .WithAgent(new DDPG
            {
                ActorHiddenLayers = new[] { 16 },
                CriticHiddenLayers = new[] { 16 },
                BatchSize = 16,
                MinBufferSize = 32,
                ActionScale = 2.0
            })
            .WithPolicy(new OrnsteinUhlenbeck(seed: 42)
            {
                Sigma = 0.2,
                Theta = 0.15
            })
            .WithReplayBuffer(5000, seed: 42)
            .WithEpisodes(maxEpisodes: 20, maxStepsPerEpisode: 50)
            .WithSeed(42)
            .Run();

        Assert.AreEqual("DDPG", result.AgentName);
        Assert.AreEqual(20, result.TotalEpisodes);
    }

    [TestMethod]
    public void RLExperiment_DDPG_WithEvaluation_ShouldRecordEvals()
    {
        var result = RLExperiment
            .For(new Pendulum { MaxSteps = 30 })
            .WithAgent(new DDPG
            {
                ActorHiddenLayers = new[] { 8 },
                CriticHiddenLayers = new[] { 8 },
                BatchSize = 8,
                MinBufferSize = 16,
                ActionScale = 2.0
            })
            .WithPolicy(new GaussianNoise(seed: 42) { Sigma = 0.1 })
            .WithReplayBuffer(2000, seed: 42)
            .WithEpisodes(maxEpisodes: 40, maxStepsPerEpisode: 30)
            .WithEvaluation(evalEpisodes: 3, evalInterval: 20)
            .WithSeed(42)
            .Run();

        // 40 episodes / 20 interval = 2 evaluations
        Assert.AreEqual(2, result.Training.EvalReturns.Count);
    }

    // ══════════════════════════════════════════════════════════════
    //  Phase 5 – Grid Search, Evaluation & Monte Carlo
    // ══════════════════════════════════════════════════════════════

    // ── RLPipelineGrid Tests ────────────────────────────────────

    [TestMethod]
    public void RLPipelineGrid_SingleAgent_NoParams_ShouldExpandToOne()
    {
        var env = new GridWorld(3, 3);
        var grid = new RLPipelineGrid()
            .AddAgent<QLearning>(() => new QLearning(9, 4, env.StateToIndex), g => { });

        var configs = grid.Expand();
        Assert.AreEqual(1, configs.Count);
        Assert.AreEqual("QLearning", configs[0].AgentTypeName);
    }

    [TestMethod]
    public void RLPipelineGrid_SingleAgent_WithParams_ShouldExpandCartesian()
    {
        var env = new GridWorld(3, 3);
        var grid = new RLPipelineGrid()
            .AddAgent<QLearning>(() => new QLearning(9, 4, env.StateToIndex), g => g
                .Add("LearningRate", 0.1, 0.5)
                .Add("Gamma", 0.9, 0.99));

        var configs = grid.Expand();
        // 2 x 2 = 4 configurations
        Assert.AreEqual(4, configs.Count);
    }

    [TestMethod]
    public void RLPipelineGrid_MultipleAgents_ShouldExpandAll()
    {
        var env = new GridWorld(3, 3);
        var grid = new RLPipelineGrid()
            .AddAgent<QLearning>(() => new QLearning(9, 4, env.StateToIndex), g => g.Add("LearningRate", 0.1, 0.5))
            .AddAgent<SARSA>(() => new SARSA(9, 4, env.StateToIndex), g => g.Add("LearningRate", 0.1));

        var configs = grid.Expand();
        // QLearning: 2, SARSA: 1 → 3 total
        Assert.AreEqual(3, configs.Count);
        Assert.AreEqual(2, configs.Count(c => c.AgentTypeName == "QLearning"));
        Assert.AreEqual(1, configs.Count(c => c.AgentTypeName == "SARSA"));
    }

    [TestMethod]
    public void RLAgentConfig_CreateAgent_ShouldSetHyperParameters()
    {
        var env = new GridWorld(3, 3);
        var grid = new RLPipelineGrid()
            .AddAgent<QLearning>(() => new QLearning(9, 4, env.StateToIndex), g => g
                .Add("LearningRate", 0.42)
                .Add("Gamma", 0.95));

        var config = grid.Expand().Single();
        var agent = config.CreateAgent();

        var p = agent.GetHyperParameters();
        Assert.AreEqual(0.42, (double)p["LearningRate"], 1e-9);
        Assert.AreEqual(0.95, (double)p["Gamma"], 1e-9);
    }

    [TestMethod]
    public void RLAgentConfig_Description_ShouldContainAgentNameAndParams()
    {
        var env = new GridWorld(3, 3);
        var grid = new RLPipelineGrid()
            .AddAgent<QLearning>(() => new QLearning(9, 4, env.StateToIndex), g => g.Add("LearningRate", 0.1));

        var config = grid.Expand().Single();
        Assert.IsTrue(config.Description.Contains("QLearning"));
        Assert.IsTrue(config.Description.Contains("LearningRate"));
    }

    // ── EpisodeEvaluator Tests ──────────────────────────────────

    [TestMethod]
    public void EpisodeEvaluator_ShouldReturnCorrectEpisodeCount()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex);
        var evaluator = new EpisodeEvaluator(env, maxStepsPerEpisode: 50);

        var result = evaluator.Evaluate(agent, numEpisodes: 5, seed: 42);

        Assert.AreEqual(5, result.NumEpisodes);
        Assert.AreEqual(5, result.Returns.Count);
    }

    [TestMethod]
    public void EpisodeEvaluator_EvalResult_ShouldComputeStatistics()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex);
        var evaluator = new EpisodeEvaluator(env, maxStepsPerEpisode: 50);

        var result = evaluator.Evaluate(agent, numEpisodes: 20, seed: 42);

        Assert.AreEqual(20, result.NumEpisodes);
        Assert.IsTrue(result.StdDev >= 0);
        // Use tolerance for floating-point comparisons
        Assert.IsTrue(result.MinReturn <= result.MeanReturn + 1e-10);
        Assert.IsTrue(result.MeanReturn <= result.MaxReturn + 1e-10);
    }

    [TestMethod]
    public void EvalResult_ConfidenceInterval_ShouldContainMean()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex);
        var evaluator = new EpisodeEvaluator(env, maxStepsPerEpisode: 50);

        var result = evaluator.Evaluate(agent, numEpisodes: 20, seed: 42);
        var (lower, upper) = result.ConfidenceInterval(0.95);

        Assert.IsTrue(lower <= result.MeanReturn);
        Assert.IsTrue(upper >= result.MeanReturn);
    }

    [TestMethod]
    public void EpisodeEvaluator_TrainedAgent_ShouldPerformBetter()
    {
        var env = new GridWorld(5, 5);
        var agent = new QLearning(25, 4, env.StateToIndex)
        {
            LearningRate = 0.2,
            Gamma = 0.99,
        };

        // Evaluate untrained
        var evaluator = new EpisodeEvaluator(env, maxStepsPerEpisode: 100);
        var beforeEval = evaluator.Evaluate(agent, numEpisodes: 10, seed: 42);

        // Train
        RLExperiment
            .For(env)
            .WithAgent(agent)
            .WithPolicy(new EpsilonGreedy(seed: 42) { Epsilon = 1.0, EpsilonMin = 0.01, EpsilonDecay = 0.995 })
            .WithEpisodes(500, 100)
            .WithSeed(42)
            .Run();

        // Evaluate trained
        var afterEval = evaluator.Evaluate(agent, numEpisodes: 10, seed: 42);

        Assert.IsTrue(afterEval.MeanReturn > beforeEval.MeanReturn,
            $"Trained mean {afterEval.MeanReturn} should beat untrained {beforeEval.MeanReturn}");
    }

    // ── Grid Search Tests ───────────────────────────────────────

    [TestMethod]
    public void RunGrid_ShouldRankConfigsByEvalReturn()
    {
        var env = new GridWorld(3, 3);
        var grid = new RLPipelineGrid()
            .AddAgent<QLearning>(() => new QLearning(9, 4, env.StateToIndex), g => g
                .Add("LearningRate", 0.1, 0.5)
                .Add("Gamma", 0.9, 0.99));

        var result = RLExperiment
            .For(env)
            .WithGrid(grid, evalEpisodes: 5)
            .WithPolicyFactory(() => new EpsilonGreedy(seed: 42) { Epsilon = 1.0, EpsilonMin = 0.01, EpsilonDecay = 0.99 })
            .WithEpisodes(200, 50)
            .WithSeed(42)
            .RunGrid();

        Assert.AreEqual(4, result.Rankings.Count);
        Assert.IsNotNull(result.Best);
        Assert.IsTrue(result.TotalDuration.TotalMilliseconds > 0);

        // Rankings should be descending by eval return
        for (int i = 1; i < result.Rankings.Count; i++)
            Assert.IsTrue(result.Rankings[i - 1].AverageEvalReturn >= result.Rankings[i].AverageEvalReturn);
    }

    [TestMethod]
    public void RunGrid_ShouldPopulateEntryFields()
    {
        var env = new GridWorld(3, 3);
        var grid = new RLPipelineGrid()
            .AddAgent<QLearning>(() => new QLearning(9, 4, env.StateToIndex), g => g.Add("LearningRate", 0.3));

        var result = RLExperiment
            .For(env)
            .WithGrid(grid, evalEpisodes: 3)
            .WithPolicyFactory(() => new EpsilonGreedy(seed: 42) { Epsilon = 1.0, EpsilonMin = 0.01, EpsilonDecay = 0.99 })
            .WithEpisodes(100, 50)
            .WithSeed(42)
            .RunGrid();

        var entry = result.Rankings.Single();
        Assert.AreEqual("QLearning", entry.AgentName);
        Assert.IsNotNull(entry.Parameters);
        Assert.IsNotNull(entry.EvalResult);
        Assert.IsNotNull(entry.TrainingResult);
        Assert.IsTrue(entry.Duration.TotalMilliseconds > 0);
        Assert.IsTrue(entry.Description.Contains("QLearning"));
    }

    [TestMethod]
    public void RunGrid_MultipleAgentTypes_ShouldWork()
    {
        var env = new GridWorld(3, 3);
        var grid = new RLPipelineGrid()
            .AddAgent<QLearning>(() => new QLearning(9, 4, env.StateToIndex), g => g.Add("LearningRate", 0.3))
            .AddAgent<SARSA>(() => new SARSA(9, 4, env.StateToIndex), g => g.Add("LearningRate", 0.3));

        var result = RLExperiment
            .For(env)
            .WithGrid(grid, evalEpisodes: 3)
            .WithPolicyFactory(() => new EpsilonGreedy(seed: 42) { Epsilon = 1.0, EpsilonMin = 0.01, EpsilonDecay = 0.99 })
            .WithEpisodes(100, 50)
            .WithSeed(42)
            .RunGrid();

        Assert.AreEqual(2, result.Rankings.Count);
        var agentNames = result.Rankings.Select(r => r.AgentName).Distinct().ToList();
        Assert.IsTrue(agentNames.Contains("QLearning"));
        Assert.IsTrue(agentNames.Contains("SARSA"));
    }

    // ── Monte Carlo Evaluation Tests ────────────────────────────

    [TestMethod]
    public void RunMonteCarlo_ShouldReturnMultipleRuns()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex)
        {
            LearningRate = 0.3,
            Gamma = 0.95
        };

        var mc = RLExperiment
            .For(env)
            .WithAgent(agent)
            .WithPolicy(new EpsilonGreedy(seed: 42) { Epsilon = 1.0, EpsilonMin = 0.01, EpsilonDecay = 0.99 })
            .WithEpisodes(100, 50)
            .WithMonteCarloEvaluation(runs: 3, evalEpisodesPerRun: 5)
            .RunMonteCarlo();

        Assert.AreEqual(3, mc.NumRuns);
        Assert.AreEqual(3, mc.Runs.Count);
        Assert.AreEqual(3, mc.EvalResults.Count);
    }

    [TestMethod]
    public void RunMonteCarlo_ShouldComputeAggregateStatistics()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex)
        {
            LearningRate = 0.3,
            Gamma = 0.95
        };

        var mc = RLExperiment
            .For(env)
            .WithAgent(agent)
            .WithPolicy(new EpsilonGreedy(seed: 42) { Epsilon = 1.0, EpsilonMin = 0.01, EpsilonDecay = 0.99 })
            .WithEpisodes(100, 50)
            .WithMonteCarloEvaluation(runs: 5, evalEpisodesPerRun: 5)
            .RunMonteCarlo();

        Assert.IsFalse(double.IsNaN(mc.MeanReturn));
        Assert.IsTrue(mc.StdDev >= 0);

        var (lower, upper) = mc.ConfidenceInterval(0.95);
        Assert.IsTrue(lower <= mc.MeanReturn);
        Assert.IsTrue(upper >= mc.MeanReturn);
    }

    [TestMethod]
    public void RunMonteCarlo_ShouldUseDifferentSeedsPerRun()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex)
        {
            LearningRate = 0.3,
            Gamma = 0.95
        };

        var mc = RLExperiment
            .For(env)
            .WithAgent(agent)
            .WithPolicy(new EpsilonGreedy(seed: 42) { Epsilon = 1.0, EpsilonMin = 0.01, EpsilonDecay = 0.99 })
            .WithEpisodes(100, 50)
            .WithMonteCarloEvaluation(runs: 3, evalEpisodesPerRun: 5)
            .WithSeed(42)
            .RunMonteCarlo();

        // Each run should potentially produce different results (different seeds)
        var returns = mc.EvalResults.Select(e => e.MeanReturn).ToList();
        Assert.AreEqual(3, returns.Count);
    }

    // ── EvalResult Median Test ──────────────────────────────────

    [TestMethod]
    public void EvalResult_MedianReturn_ShouldBeCorrect()
    {
        // Construct an EvalResult directly with known returns
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex);
        var evaluator = new EpisodeEvaluator(env, maxStepsPerEpisode: 50);

        var result = evaluator.Evaluate(agent, numEpisodes: 10, seed: 42);
        var sorted = result.Returns.OrderBy(r => r).ToList();
        double expectedMedian = (sorted[4] + sorted[5]) / 2.0;

        Assert.AreEqual(expectedMedian, result.MedianReturn, 1e-9);
    }

    // ══════════════════════════════════════════════════════════════
    //  Phase 6 – Diagnostics & Visualisation Hooks
    // ══════════════════════════════════════════════════════════════

    // ── QValueHeatmap Tests ─────────────────────────────────────

    [TestMethod]
    public void QValueHeatmap_MaxQValues_ShouldReturnOnePerState()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex);

        var heatmap = QValueHeatmap.GetMaxQValues(agent, env);

        Assert.AreEqual(9, heatmap.Count);
        // Untrained: all Q=0 → max=0
        foreach (var s in heatmap)
            Assert.AreEqual(0.0, s.Value, 1e-9);
    }

    [TestMethod]
    public void QValueHeatmap_TrainedAgent_ShouldHaveNonZeroValues()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex)
        {
            LearningRate = 0.2,
            Gamma = 0.99
        };

        RLExperiment
            .For(env)
            .WithAgent(agent)
            .WithPolicy(new EpsilonGreedy(seed: 42) { Epsilon = 1.0, EpsilonMin = 0.01, EpsilonDecay = 0.99 })
            .WithEpisodes(300, 50)
            .WithSeed(42)
            .Run();

        var heatmap = QValueHeatmap.GetMaxQValues(agent, env);

        // At least some states should have learned non-zero Q values
        Assert.IsTrue(heatmap.Any(s => Math.Abs(s.Value) > 0.01),
            "Trained agent should have non-zero Q values");
    }

    [TestMethod]
    public void QValueHeatmap_QValuesForAction_ShouldReturnCorrectCount()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex);

        for (int action = 0; action < 4; action++)
        {
            var actionValues = QValueHeatmap.GetQValuesForAction(agent, env, action);
            Assert.AreEqual(9, actionValues.Count);
        }
    }

    [TestMethod]
    public void QValueHeatmap_GetQTableMatrix_ShouldMatchDimensions()
    {
        var env = new GridWorld(4, 4);
        var agent = new QLearning(16, 4, env.StateToIndex);

        var matrix = QValueHeatmap.GetQTableMatrix(agent);

        Assert.AreEqual(16, matrix.rowLength);
        Assert.AreEqual(4, matrix.columnLength);
    }

    [TestMethod]
    public void QValueHeatmap_GreedyPolicy_ShouldReturnValidActions()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex);

        var policy = QValueHeatmap.GetGreedyPolicy(agent, env);

        Assert.AreEqual(9, policy.Count);
        foreach (var s in policy)
            Assert.IsTrue(s.Value >= 0 && s.Value < 4);
    }

    // ── PolicyVisualizer Tests ──────────────────────────────────

    [TestMethod]
    public void PolicyVisualizer_SoftmaxProbabilities_ShouldSumToOne()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex);

        var probs = PolicyVisualizer.GetSoftmaxProbabilities(agent, env, temperature: 1.0);

        Assert.AreEqual(4, probs.Count); // one list per action
        // For each state, probabilities across actions should sum to ~1
        for (int stateIdx = 0; stateIdx < 9; stateIdx++)
        {
            double sum = 0;
            for (int a = 0; a < 4; a++)
                sum += probs[a][stateIdx].Value;
            Assert.AreEqual(1.0, sum, 1e-6);
        }
    }

    [TestMethod]
    public void PolicyVisualizer_Entropy_UntrainedAgent_ShouldBeMaximal()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex);

        var entropy = PolicyVisualizer.GetPolicyEntropy(env, state =>
        {
            var q = agent.GetQValues(state);
            return Softmax(q, 1.0);
        });

        double maxEntropy = Math.Log(4); // uniform over 4 actions
        Assert.AreEqual(9, entropy.Count);
        // Untrained Q=0 → uniform → max entropy
        foreach (var s in entropy)
            Assert.AreEqual(maxEntropy, s.Value, 1e-6);
    }

    [TestMethod]
    public void PolicyVisualizer_DominantAction_ShouldReturnValidActions()
    {
        var env = new GridWorld(3, 3);
        var agent = new QLearning(9, 4, env.StateToIndex);

        var dominant = PolicyVisualizer.GetDominantAction(env, state =>
        {
            var q = agent.GetQValues(state);
            return Softmax(q, 1.0);
        });

        Assert.AreEqual(9, dominant.Count);
        foreach (var s in dominant)
            Assert.IsTrue(s.Value >= 0 && s.Value < 4);
    }

    [TestMethod]
    public void PolicyVisualizer_ActionProbabilities_REINFORCE_ShouldWork()
    {
        var agent = new REINFORCE
        {
            HiddenLayers = new[] { 8 },
            LearningRate = 0.01
        };
        agent.Initialize(2, 4, seed: 42);
        var env = new GridWorld(3, 3);

        var probs = PolicyVisualizer.GetActionProbabilities(env,
            state => agent.GetActionProbabilities(state));

        Assert.AreEqual(4, probs.Count);
        for (int stateIdx = 0; stateIdx < 9; stateIdx++)
        {
            double sum = 0;
            for (int a = 0; a < 4; a++)
                sum += probs[a][stateIdx].Value;
            Assert.AreEqual(1.0, sum, 1e-6);
        }
    }

    // ── ValueFunctionSurface Tests ──────────────────────────────

    [TestMethod]
    public void ValueFunctionSurface_Sample1D_ShouldReturnCorrectCount()
    {
        Func<VectorN, double> linear = s => s[0] * 2;

        var curve = ValueFunctionSurface.Sample1D(linear, min: -1, max: 1, numPoints: 50);

        Assert.AreEqual(50, curve.Count);
        // First point at x=-1 → value=-2, last at x=1 → value=2
        Assert.AreEqual(-2.0, curve[0].Value, 1e-6);
        Assert.AreEqual(2.0, curve[49].Value, 1e-6);
    }

    [TestMethod]
    public void ValueFunctionSurface_Sample2D_ShouldReturnGrid()
    {
        Func<VectorN, double> sum = s => s[0] + s[1];

        var surface = ValueFunctionSurface.Sample2D(sum,
            minX: 0, maxX: 1, numX: 5,
            minY: 0, maxY: 1, numY: 5);

        Assert.AreEqual(25, surface.Points.Count);
        Assert.AreEqual(5, surface.NumX);
        Assert.AreEqual(5, surface.NumY);

        // Corner (0,0) → 0, corner (1,1) → 2
        Assert.AreEqual(0.0, surface.ValueAt(0, 0), 1e-6);
        Assert.AreEqual(2.0, surface.ValueAt(4, 4), 1e-6);
    }

    [TestMethod]
    public void ValueFunctionSurface_ToMatrix_ShouldHaveCorrectDimensions()
    {
        Func<VectorN, double> constant = _ => 1.0;

        var surface = ValueFunctionSurface.Sample2D(constant,
            minX: 0, maxX: 1, numX: 3,
            minY: 0, maxY: 1, numY: 4);

        var matrix = surface.ToMatrix();
        Assert.AreEqual(3, matrix.rowLength);
        Assert.AreEqual(4, matrix.columnLength);
    }

    [TestMethod]
    public void ValueFunctionSurface_DQNMaxQ_ShouldReturnValues()
    {
        var env = new CartPole();
        var dqn = new DQN
        {
            HiddenLayers = new[] { 8 },
            LearningRate = 0.001,
            BatchSize = 8,
            MinBufferSize = 16
        };
        dqn.Initialize(env.ObservationSize, env.ActionSize, seed: 42);

        var maxQFn = ValueFunctionSurface.MaxQFunction(dqn);
        var curve = ValueFunctionSurface.Sample1D(
            s => maxQFn(new VectorN(new[] { s[0], 0, 0, 0 })),
            min: -2.4, max: 2.4, numPoints: 20);

        Assert.AreEqual(20, curve.Count);
        // Values should be finite
        Assert.IsTrue(curve.All(s => !double.IsNaN(s.Value) && !double.IsInfinity(s.Value)));
    }

    [TestMethod]
    public void ValueFunctionSurface_ActorCriticValue_ShouldReturnValues()
    {
        var env = new CartPole();
        var ac = new ActorCritic
        {
            ActorHiddenLayers = new[] { 8 },
            CriticHiddenLayers = new[] { 8 }
        };
        ac.Initialize(env.ObservationSize, env.ActionSize, seed: 42);

        var vFn = ValueFunctionSurface.ValueFunction(ac);
        var curve = ValueFunctionSurface.Sample1D(
            s => vFn(new VectorN(new[] { s[0], 0, 0, 0 })),
            min: -2.4, max: 2.4, numPoints: 20);

        Assert.AreEqual(20, curve.Count);
        Assert.IsTrue(curve.All(s => !double.IsNaN(s.Value) && !double.IsInfinity(s.Value)));
    }

    [TestMethod]
    public void ValueFunctionSurface_PPOValue_ShouldReturnValues()
    {
        var env = new CartPole();
        var ppo = new PPO
        {
            ActorHiddenLayers = new[] { 8 },
            CriticHiddenLayers = new[] { 8 }
        };
        ppo.Initialize(env.ObservationSize, env.ActionSize, seed: 42);

        var vFn = ValueFunctionSurface.ValueFunction(ppo);
        var curve = ValueFunctionSurface.Sample1D(
            s => vFn(new VectorN(new[] { s[0], 0, 0, 0 })),
            min: -2.4, max: 2.4, numPoints: 20);

        Assert.AreEqual(20, curve.Count);
        Assert.IsTrue(curve.All(s => !double.IsNaN(s.Value) && !double.IsInfinity(s.Value)));
    }

    // ── TrainingResult Curves (6.1-6.3 verification) ────────────

    [TestMethod]
    public void TrainingResult_Curves_ShouldMatchEpisodeCount()
    {
        var env = new GridWorld(3, 3);
        var result = RLExperiment
            .For(env)
            .WithAgent(new QLearning(9, 4, env.StateToIndex))
            .WithPolicy(new EpsilonGreedy(seed: 42))
            .WithEpisodes(50, 30)
            .WithSeed(42)
            .Run();

        Assert.AreEqual(50, result.Training.ReturnCurve.Count);
        Assert.AreEqual(50, result.Training.ExplorationCurve.Count);
        // Return curve indices should be 0..49
        Assert.AreEqual(0.0, result.Training.ReturnCurve[0].Index);
        Assert.AreEqual(49.0, result.Training.ReturnCurve[49].Index);
    }

    // ── Helper ───────────────────────────────────────────────────

    private static VectorN Softmax(VectorN values, double temperature)
    {
        double max = values[0];
        for (int i = 1; i < values.Length; i++)
            if (values[i] > max) max = values[i];

        var exps = new double[values.Length];
        double sum = 0;
        for (int i = 0; i < values.Length; i++)
        {
            exps[i] = Math.Exp((values[i] - max) / temperature);
            sum += exps[i];
        }
        for (int i = 0; i < exps.Length; i++)
            exps[i] /= sum;

        return new VectorN(exps);
    }
}
