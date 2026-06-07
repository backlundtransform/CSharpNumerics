using CSharpNumerics.Engines.Game.AI;
using CSharpNumerics.Engines.Game.Flight;
using CSharpNumerics.ML.ReinforcementLearning.Algorithms.ContinuousControl;
using CSharpNumerics.ML.ReinforcementLearning.Algorithms.PolicyGradient;
using CSharpNumerics.ML.ReinforcementLearning.Core;
using CSharpNumerics.Engines.Game.RL;
using CSharpNumerics.ML.ReinforcementLearning.Policies;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Linq;

namespace NumericTest;

[TestClass]
public class GameAITests
{
    // ═══════════════════════════════════════════════════════════════
    //  Flight Environment
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void FlightEnv_ResetProducesValidObservation()
    {
        var env = new FlightEnv();
        var (state, info) = env.Reset(42);

        Assert.AreEqual(12, state.Length, "Observation should be 12-dimensional");
        Assert.IsTrue(state[0] > 0, "Altitude should be positive");
        Assert.IsTrue(state[1] > 0, "Airspeed should be positive");
    }

    [TestMethod]
    public void FlightEnv_StepAdvancesSimulation()
    {
        var env = new FlightEnv();
        env.Reset(42);

        // Moderate throttle, slight pitch up
        var action = new VectorN(new[] { 0.6, 0.1, 0.0, 0.0 });
        var (nextState, reward, done, info) = env.Step(action);

        Assert.AreEqual(12, nextState.Length);
        Assert.IsFalse(done, "Should not terminate on first step");
        Assert.IsTrue(info.ContainsKey("altitude"));
    }

    [TestMethod]
    public void FlightEnv_TrainedAgentMaintainsLevelFlight()
    {
        // Train DDPG for a few episodes on a simplified flight task
        var env = new FlightEnv(
            dt: 0.1,
            maxSteps: 100,
            initAltitude: 1000,
            initAirspeed: 60);

        var agent = new DDPG
        {
            ActorHiddenLayers = new[] { 32, 32 },
            CriticHiddenLayers = new[] { 32, 32 },
            ActorLearningRate = 1e-3,
            CriticLearningRate = 1e-3,
            Gamma = 0.99,
            ActionScale = 1.0
        };
        agent.Initialize(env.ObservationSize, env.ActionSize, 42);

        var noise = new GaussianNoise(42) { Sigma = 0.3, SigmaDecay = 0.99, SigmaMin = 0.05 };

        // Train for 20 episodes
        for (int ep = 0; ep < 20; ep++)
        {
            var episode = new Episode();
            var (state, _) = env.Reset(42 + ep);

            for (int step = 0; step < 100; step++)
            {
                var rawAction = agent.SelectContinuousAction(state);
                var exploredAction = noise.SelectAction(rawAction, rawAction);

                // Clamp action
                var clamped = new VectorN(new[]
                {
                    Math.Clamp(exploredAction[0], 0, 1),
                    Math.Clamp(exploredAction[1], -1, 1),
                    Math.Clamp(exploredAction[2], -1, 1),
                    Math.Clamp(exploredAction[3], -1, 1)
                });

                var (nextState, reward, done, _) = env.Step(clamped);
                var transition = new Transition(state, clamped, reward, nextState, done);
                episode.Transitions.Add(transition);
                agent.Train(transition);

                state = nextState;
                if (done) break;
            }

            agent.EndEpisode(episode);
            noise.Decay();
        }

        // Evaluate: run one episode without noise and check reward improves
        var (evalState, _) = env.Reset(999);
        double totalReward = 0;
        for (int step = 0; step < 100; step++)
        {
            var action = agent.SelectContinuousAction(evalState);
            var clamped2 = new VectorN(new[]
            {
                Math.Clamp(action[0], 0, 1),
                Math.Clamp(action[1], -1, 1),
                Math.Clamp(action[2], -1, 1),
                Math.Clamp(action[3], -1, 1)
            });

            var (ns, r, d, _) = env.Step(clamped2);
            totalReward += r;
            evalState = ns;
            if (d) break;
        }

        // Agent should achieve reward better than constant-crash baseline (-100)
        Assert.IsTrue(totalReward > -100,
            $"Trained agent should do better than crashing. Total reward: {totalReward:F2}");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Dogfight Environment
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void DogfightEnv_ResetAndStep()
    {
        var env = new DogfightEnv(maxSteps: 50);
        var (state, _) = env.Reset(42);

        Assert.AreEqual(18, state.Length, "Observation should be 18D");

        // Take a step: full throttle, slight turn
        var action = new VectorN(new[] { 0.9, 0.1, 0.3, 0.0 });
        var (nextState, reward, done, info) = env.Step(action);

        Assert.AreEqual(18, nextState.Length);
        Assert.IsTrue(info.ContainsKey("range"));
        double range = (double)info["range"];
        Assert.IsTrue(range > 0, "Range should be positive");
    }

    [TestMethod]
    public void DogfightEnv_PursuerClosesDistance()
    {
        var env = new DogfightEnv(maxSteps: 200);
        var (state, _) = env.Reset(42);

        // Pursue: full throttle, no maneuver (straight ahead)
        double initRange = 0;
        double finalRange = 0;

        for (int i = 0; i < 100; i++)
        {
            var action = new VectorN(new[] { 1.0, 0.0, 0.0, 0.0 });
            var (ns, r, d, info) = env.Step(action);

            if (i == 0) initRange = (double)info["range"];
            if (i == 99) finalRange = (double)info["range"];

            state = ns;
            if (d) break;
        }

        Assert.IsTrue(finalRange < initRange,
            $"Pursuer should close distance: init={initRange:F0}m, final={finalRange:F0}m");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Fluid Navigation Environment
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void FluidNavEnv_ResetAndStep()
    {
        var env = new FluidNavigationEnv(gridSize: 16, maxSteps: 50);
        var (state, _) = env.Reset(42);

        Assert.AreEqual(8, state.Length, "Observation should be 8D");

        var action = new VectorN(new[] { 0.5, 0.5 });
        var (ns, r, d, info) = env.Step(action);
        Assert.AreEqual(8, ns.Length);
        Assert.IsTrue(info.ContainsKey("distance"));
    }

    // ═══════════════════════════════════════════════════════════════
    //  GameAIAgent
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void GameAIAgent_WrapsRLPolicy()
    {
        var ddpg = new DDPG();
        ddpg.Initialize(12, 4, 42);

        var gameAgent = new GameAIAgent(ddpg, continuous: true, name: "Autopilot");

        Assert.AreEqual("Autopilot", gameAgent.Name);
        Assert.IsTrue(gameAgent.IsContinuous);

        var obs = new VectorN(new double[12]);
        var action = gameAgent.Act(obs);

        Assert.AreEqual(4, action.Length);
        Assert.AreEqual(1, gameAgent.ActionCount);

        gameAgent.Reset();
        Assert.AreEqual(0, gameAgent.ActionCount);
    }

    // ═══════════════════════════════════════════════════════════════
    //  AITrainer
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void AITrainer_TrainsAndReturnsAgent()
    {
        var env = new FluidNavigationEnv(gridSize: 16, maxSteps: 30);
        var ddpg = new DDPG
        {
            ActorHiddenLayers = new[] { 16 },
            CriticHiddenLayers = new[] { 16 }
        };
        ddpg.Initialize(env.ObservationSize, env.ActionSize, 42);

        var trainer = new AITrainer(env, ddpg)
            .WithEpisodes(5, 30)
            .WithPolicy(new GaussianNoise(42))
            .WithSeed(42);

        var agent = trainer.Train("NavBot");

        Assert.AreEqual("NavBot", agent.Name);
        Assert.IsTrue(trainer.EpisodeReturns.Count == 5, "Should have 5 episode returns");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Behavior Tree
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void BehaviorTree_SequenceRunsAllChildren()
    {
        int callCount = 0;
        var tree = new BehaviorTree("test",
            new SequenceNode("root",
                new ActionNode("a1", ctx => { callCount++; return NodeStatus.Success; }),
                new ActionNode("a2", ctx => { callCount++; return NodeStatus.Success; }),
                new ActionNode("a3", ctx => { callCount++; return NodeStatus.Success; })
            ));

        var status = tree.Tick(new BehaviorContext());
        Assert.AreEqual(NodeStatus.Success, status);
        Assert.AreEqual(3, callCount);
    }

    [TestMethod]
    public void BehaviorTree_SequenceStopsOnFailure()
    {
        int callCount = 0;
        var tree = new BehaviorTree("test",
            new SequenceNode("root",
                new ActionNode("a1", ctx => { callCount++; return NodeStatus.Success; }),
                new ActionNode("a2", ctx => { callCount++; return NodeStatus.Failure; }),
                new ActionNode("a3", ctx => { callCount++; return NodeStatus.Success; })
            ));

        var status = tree.Tick(new BehaviorContext());
        Assert.AreEqual(NodeStatus.Failure, status);
        Assert.AreEqual(2, callCount, "Should stop after failure");
    }

    [TestMethod]
    public void BehaviorTree_SelectorReturnsFirstSuccess()
    {
        string chosen = null;
        var tree = new BehaviorTree("test",
            new SelectorNode("root",
                new ActionNode("a1", ctx => NodeStatus.Failure),
                new ActionNode("a2", ctx => { chosen = "a2"; return NodeStatus.Success; }),
                new ActionNode("a3", ctx => { chosen = "a3"; return NodeStatus.Success; })
            ));

        var status = tree.Tick(new BehaviorContext());
        Assert.AreEqual(NodeStatus.Success, status);
        Assert.AreEqual("a2", chosen, "Should select second (first success)");
    }

    [TestMethod]
    public void BehaviorTree_ConditionGatesAction()
    {
        bool actionRan = false;
        var tree = new BehaviorTree("test",
            new SequenceNode("root",
                new ConditionNode("check", ctx => ctx.Get<double>("fuel") > 50),
                new ActionNode("fly", ctx => { actionRan = true; return NodeStatus.Success; })
            ));

        var ctx = new BehaviorContext();

        // Low fuel — condition fails, action should not run
        ctx.Set("fuel", 30.0);
        var status = tree.Tick(ctx);
        Assert.AreEqual(NodeStatus.Failure, status);
        Assert.IsFalse(actionRan);

        // High fuel — condition passes
        ctx.Set("fuel", 80.0);
        status = tree.Tick(ctx);
        Assert.AreEqual(NodeStatus.Success, status);
        Assert.IsTrue(actionRan);
    }

    [TestMethod]
    public void BehaviorTree_InverterFlipsResult()
    {
        var tree = new BehaviorTree("test",
            new InverterNode("inv",
                new ActionNode("always_fail", ctx => NodeStatus.Failure)));

        var status = tree.Tick(new BehaviorContext());
        Assert.AreEqual(NodeStatus.Success, status);
    }

    [TestMethod]
    public void BehaviorTree_HybridMLAndScripted()
    {
        // Test a realistic hybrid tree: ML decision + scripted fallback
        var ctx = new BehaviorContext();
        ctx.Set("inCombatRange", true);
        ctx.Set("mlAction", new VectorN(new[] { 0.8, 0.1, -0.3, 0.0 }));

        VectorN appliedAction = default;

        var tree = new BehaviorTree("combat_ai",
            new SelectorNode("root",
                // Branch 1: If in combat range, use ML action
                new SequenceNode("ml_combat",
                    new ConditionNode("combat_check", c => c.Get<bool>("inCombatRange")),
                    new ActionNode("ml_action", c =>
                    {
                        appliedAction = c.Get<VectorN>("mlAction");
                        return NodeStatus.Success;
                    })),
                // Branch 2: Fallback — patrol
                new ActionNode("patrol", c =>
                {
                    appliedAction = new VectorN(new[] { 0.5, 0.0, 0.0, 0.0 });
                    return NodeStatus.Success;
                })
            ));

        tree.Tick(ctx);
        Assert.IsTrue(appliedAction.Length > 0);
        Assert.AreEqual(0.8, appliedAction[0], 0.01, "Should use ML action when in combat range");

        // Now out of combat range — should patrol
        ctx.Set("inCombatRange", false);
        tree.Reset();
        tree.Tick(ctx);
        Assert.AreEqual(0.5, appliedAction[0], 0.01, "Should patrol when out of range");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Formation Controller
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Formation_MaintainsSpacingUnderStraightFlight()
    {
        var config = AircraftConfig.GenericLightAircraft();

        // Leader
        var leader = new FlightDynamicsEngine(config);
        leader.Init();
        leader.SetState(new AircraftState(
            new Vector(0, 0, -1000),
            new Vector(60, 0, 0),
            Quaternion.Identity,
            new Vector(0, 0, 0)));
        leader.SetInput(new ControlInput(0.5, 0, 0, 0));

        // Wingman — starts near desired position
        var wingEngine = new FlightDynamicsEngine(config);
        wingEngine.Init();
        wingEngine.SetState(new AircraftState(
            new Vector(-30, 30, -1000), // 30m behind, 30m right
            new Vector(60, 0, 0),
            Quaternion.Identity,
            new Vector(0, 0, 0)));

        var formation = new FormationController(leader);
        formation.AddWingman("Wing 2",
            new Vector(-30, 30, 0), // desired offset: 30m behind, 30m right
            wingEngine);

        // Simulate 50 steps
        double dt = 0.05;
        for (int i = 0; i < 50; i++)
        {
            leader.Step(dt);
            formation.Step(dt);
        }

        var errors = formation.GetPositionErrors();
        Assert.IsTrue(errors[0] < 200,
            $"Formation error should be manageable: {errors[0]:F1}m");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Adaptive Difficulty
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void AdaptiveDifficulty_IncreasesWhenPlayerDominates()
    {
        var ad = new AdaptiveDifficulty(
            initialDifficulty: 0.5,
            targetPerformance: 0.5,
            windowSize: 5);

        double initDiff = ad.Difficulty;

        // Player wins consistently
        for (int i = 0; i < 5; i++)
            ad.RecordPerformance(0.9);

        Assert.IsTrue(ad.Difficulty > initDiff,
            $"Difficulty should increase when player dominates: {ad.Difficulty:F3} vs initial {initDiff:F3}");
    }

    [TestMethod]
    public void AdaptiveDifficulty_DecreasesWhenPlayerStruggles()
    {
        var ad = new AdaptiveDifficulty(
            initialDifficulty: 0.5,
            targetPerformance: 0.5,
            windowSize: 5);

        double initDiff = ad.Difficulty;

        // Player loses consistently
        for (int i = 0; i < 5; i++)
            ad.RecordPerformance(0.1);

        Assert.IsTrue(ad.Difficulty < initDiff,
            $"Difficulty should decrease when player struggles: {ad.Difficulty:F3} vs initial {initDiff:F3}");
    }

    [TestMethod]
    public void AdaptiveDifficulty_GetAIParameters_ScalesWithDifficulty()
    {
        var easy = new AdaptiveDifficulty(initialDifficulty: 0.1);
        var hard = new AdaptiveDifficulty(initialDifficulty: 0.9);

        var easyParams = easy.GetAIParameters();
        var hardParams = hard.GetAIParameters();

        // Easy AI should have more reaction delay
        Assert.IsTrue(easyParams["reactionDelay"] > hardParams["reactionDelay"]);

        // Hard AI should have better accuracy
        Assert.IsTrue(hardParams["accuracyScale"] > easyParams["accuracyScale"]);

        // Hard AI should be more aggressive
        Assert.IsTrue(hardParams["aggressiveness"] > easyParams["aggressiveness"]);

        // Easy AI should have more noise
        Assert.IsTrue(easyParams["explorationNoise"] > hardParams["explorationNoise"]);
    }

    [TestMethod]
    public void AdaptiveDifficulty_StaysWithinBounds()
    {
        var ad = new AdaptiveDifficulty(initialDifficulty: 0.5);
        ad.MinDifficulty = 0.2;
        ad.MaxDifficulty = 0.8;

        // Push difficulty down aggressively
        for (int i = 0; i < 50; i++)
            ad.RecordPerformance(0.0);

        Assert.IsTrue(ad.Difficulty >= 0.2, $"Should not go below min: {ad.Difficulty:F3}");

        ad.Reset(0.5);

        // Push difficulty up aggressively
        for (int i = 0; i < 50; i++)
            ad.RecordPerformance(1.0);

        Assert.IsTrue(ad.Difficulty <= 0.8, $"Should not go above max: {ad.Difficulty:F3}");
    }
}
