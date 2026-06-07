using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.RL;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.ML.ReinforcementLearning.Algorithms.Tabular;
using CSharpNumerics.ML.ReinforcementLearning.Algorithms.ValueBased;
using CSharpNumerics.ML.ReinforcementLearning.Experiment;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.ML.ReinforcementLearning.Policies;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Environmental.Enums;
using System;
using System.Linq;

namespace NumericTest;

[TestClass]
public class PlumeEnvironmentTests
{
    // ═══════════════════════════════════════════════════════════════
    //  PlumeEnvironment — Basic Environment Contract
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void PlumeEnvironment_Properties_AreCorrect()
    {
        var env = CreateSmallEnvironment();

        Assert.AreEqual(8, env.ObservationSize);
        Assert.AreEqual(6, env.ActionSize);
        Assert.IsTrue(env.IsDiscrete);
    }

    [TestMethod]
    public void PlumeEnvironment_Reset_ReturnsValidState()
    {
        var env = CreateSmallEnvironment();
        var (state, info) = env.Reset(42);

        Assert.IsNotNull(state);
        Assert.AreEqual(8, state.Length);
        Assert.IsTrue(info.ContainsKey("snapshot"));

        // Normalized time should be 0 at start
        Assert.AreEqual(0.0, state[7], 1e-10);
    }

    [TestMethod]
    public void PlumeEnvironment_Step_None_DoesNotIncurActionCost()
    {
        var env = CreateSmallEnvironment();
        env.Reset(42);

        var (state, reward, done, info) = env.Step((int)MitigationAction.None);

        Assert.IsNotNull(state);
        Assert.AreEqual(8, state.Length);
        // Reward should be >= -1 (exposed fraction only, no action cost)
        Assert.IsTrue(reward >= -1.0);
        Assert.IsFalse(double.IsNaN(reward));
    }

    [TestMethod]
    public void PlumeEnvironment_Step_WithAction_IncursActionCost()
    {
        var env = CreateSmallEnvironment();
        env.ActionCost = 0.1;
        env.Reset(42);

        // Step with no action
        var (_, rewardNone, _, _) = env.Step((int)MitigationAction.None);

        // Reset and step with barrier action
        env.Reset(42);
        var (_, rewardAction, _, _) = env.Step((int)MitigationAction.BarrierNorth);

        // Action reward should be lower by approximately the action cost
        Assert.IsTrue(rewardAction < rewardNone + 1e-10,
            $"Action reward ({rewardAction}) should be less than no-action reward ({rewardNone})");
    }

    [TestMethod]
    public void PlumeEnvironment_RunsFullEpisode_Terminates()
    {
        var env = CreateSmallEnvironment();
        env.Reset(42);

        bool terminated = false;
        int steps = 0;

        while (!terminated && steps < 1000)
        {
            var (_, _, done, _) = env.Step((int)MitigationAction.None);
            terminated = done;
            steps++;
        }

        Assert.IsTrue(terminated, "Episode should terminate when all time steps are exhausted.");
        Assert.AreEqual(env.MaxSteps, steps);
    }

    [TestMethod]
    public void PlumeEnvironment_ActivateFilter_ReducesEmission()
    {
        var env = CreateSmallEnvironment();
        env.Reset(42);

        // Step once with None to get baseline
        env.Step((int)MitigationAction.None);

        // Activate filter
        var (state, _, _, info) = env.Step((int)MitigationAction.ActivateFilter);

        Assert.IsTrue((bool)info["filterActive"]);
        double emissionRate = (double)info["emissionRate"];
        Assert.IsTrue(emissionRate < 5.0, "Emission rate should be reduced after filter activation.");
    }

    [TestMethod]
    public void PlumeEnvironment_Barrier_ReducesConcentration()
    {
        var env = CreateSmallEnvironment();
        env.Reset(42);

        // Step with no barrier
        var (stateNone, _, _, infoNone) = env.Step((int)MitigationAction.None);
        double maxNone = stateNone[0]; // max concentration

        // Reset and apply barrier
        env.Reset(42);
        var (stateBarrier, _, _, _) = env.Step((int)MitigationAction.BarrierNorth);
        double maxBarrier = stateBarrier[0];

        // Barrier should reduce or maintain concentration
        Assert.IsTrue(maxBarrier <= maxNone + 1e-15,
            $"Barrier max ({maxBarrier}) should be <= no-barrier max ({maxNone})");
    }

    [TestMethod]
    public void PlumeEnvironment_StepVectorN_DelegatesToInt()
    {
        var env = CreateSmallEnvironment();
        env.Reset(42);

        var (stateInt, rewardInt, doneInt, _) = env.Step(0);

        env.Reset(42);
        var (stateVec, rewardVec, doneVec, _) = env.Step(new VectorN(new[] { 0.0 }));

        Assert.AreEqual(rewardInt, rewardVec, 1e-10);
        Assert.AreEqual(doneInt, doneVec);
    }

    [TestMethod]
    public void PlumeEnvironment_StateNormalization_TimeProgressesCorrectly()
    {
        var env = CreateSmallEnvironment();
        env.Reset(42);

        // First step: time = 1/maxSteps
        var (state1, _, _, _) = env.Step(0);
        double t1 = state1[7];

        var (state2, _, _, _) = env.Step(0);
        double t2 = state2[7];

        Assert.IsTrue(t2 > t1, "Time should increase with each step.");
        Assert.IsTrue(t2 <= 1.0, "Normalized time should be <= 1.0");
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentOutOfRangeException))]
    public void PlumeEnvironment_InvalidAction_Throws()
    {
        var env = CreateSmallEnvironment();
        env.Reset(42);
        env.Step(99);
    }

    // ═══════════════════════════════════════════════════════════════
    //  ScenarioRLAnalyzer — Fluent API
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void ScenarioRLAnalyzer_SingleAgent_TrainsSuccessfully()
    {
        var env = CreateSmallEnvironment();

        var result = ScenarioRLAnalyzer
            .For(5.0, new Vector(0, 0, 50))
            .WithWind(10, new Vector(1, 0, 0))
            .WithStability(StabilityClass.D)
            .OverGrid(new GeoGrid(-50, 50, -50, 50, 0, 0, 50))
            .OverTime(0, 120, 60)
            .WithThreshold(1e-6)
            .WithAgent(new DQN
            {
                HiddenLayers = new[] { 32, 32 },
                LearningRate = 0.001,
                Gamma = 0.99,
                BatchSize = 16,
                MinBufferSize = 32
            })
            .WithPolicy(new EpsilonGreedy(seed: 42))
            .WithReplayBuffer(1000, seed: 42)
            .WithEpisodes(5, maxStepsPerEpisode: 10)
            .WithSeed(42)
            .Run();

        Assert.IsNotNull(result);
        Assert.IsNotNull(result.ExperimentResult);
        Assert.IsNotNull(result.Environment);
        Assert.AreEqual("DQN", result.AgentName);
        Assert.IsFalse(double.IsNaN(result.AverageReturn));
    }

    [TestMethod]
    public void ScenarioRLAnalyzer_EnvironmentTuning_IsApplied()
    {
        var result = ScenarioRLAnalyzer
            .For(5.0, new Vector(0, 0, 50))
            .WithWind(10, new Vector(1, 0, 0))
            .OverGrid(new GeoGrid(-50, 50, -50, 50, 0, 0, 50))
            .OverTime(0, 120, 60)
            .WithThreshold(1e-4)
            .WithActionCost(0.2)
            .WithMitigationEfficiency(barrierEfficiency: 0.6, filterEfficiency: 0.7)
            .WithAgent(new DQN
            {
                HiddenLayers = new[] { 16 },
                MinBufferSize = 16,
                BatchSize = 8
            })
            .WithPolicy(new EpsilonGreedy(seed: 42))
            .WithReplayBuffer(500, seed: 42)
            .WithEpisodes(3, maxStepsPerEpisode: 5)
            .WithSeed(42)
            .Run();

        Assert.IsInstanceOfType(result.Environment, typeof(PlumeEnvironment));
        var plume = (PlumeEnvironment)result.Environment;
        Assert.AreEqual(1e-4, plume.Threshold);
        Assert.AreEqual(0.2, plume.ActionCost);
        Assert.AreEqual(0.6, plume.BarrierEfficiency);
        Assert.AreEqual(0.7, plume.FilterEfficiency);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Integration with standard RLExperiment
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void PlumeEnvironment_WorksWithRLExperiment_Directly()
    {
        var env = new PlumeEnvironment(
            emissionRate: 5.0,
            windSpeed: 10,
            windDirection: new Vector(1, 0, 0),
            stackHeight: 50,
            sourcePosition: new Vector(0, 0, 50),
            grid: new GeoGrid(-50, 50, -50, 50, 0, 0, 50),
            timeFrame: new TimeFrame(0, 120, 60),
            stability: StabilityClass.D);

        var result = RLExperiment
            .For(env)
            .WithAgent(new DQN
            {
                HiddenLayers = new[] { 32, 32 },
                LearningRate = 0.001,
                BatchSize = 16,
                MinBufferSize = 32
            })
            .WithPolicy(new EpsilonGreedy(seed: 42))
            .WithReplayBuffer(1000, seed: 42)
            .WithEpisodes(5, maxStepsPerEpisode: 10)
            .WithSeed(42)
            .Run();

        Assert.IsNotNull(result);
        Assert.AreEqual("DQN", result.AgentName);
        Assert.IsTrue(result.TotalEpisodes > 0);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Helpers
    // ═══════════════════════════════════════════════════════════════

    private static PlumeEnvironment CreateSmallEnvironment()
    {
        return new PlumeEnvironment(
            emissionRate: 5.0,
            windSpeed: 10,
            windDirection: new Vector(1, 0, 0),
            stackHeight: 50,
            sourcePosition: new Vector(0, 0, 50),
            grid: new GeoGrid(-50, 50, -50, 50, 0, 0, 50),
            timeFrame: new TimeFrame(0, 180, 60),
            stability: StabilityClass.D);
    }
}
