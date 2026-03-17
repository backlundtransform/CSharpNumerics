using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Simulation;
using CSharpNumerics.ML.ReinforcementLearning.Experiment;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.ML.ReinforcementLearning.Policies;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Enums;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Engines.GIS.RL
{
    /// <summary>
    /// Fluent entry point for training an RL agent on a GIS plume scenario.
    /// Follows the same pattern as <see cref="Analysis.ScenarioClusterAnalyzer"/>
    /// but wraps <see cref="RLExperiment"/> instead of
    /// <see cref="ML.Clustering.ClusteringExperiment"/>.
    /// <para>
    /// Usage:
    /// <code>
    /// var result = ScenarioRLAnalyzer
    ///     .For(5.0, new Vector(0, 0, 50))
    ///     .WithWind(10, new Vector(1, 0, 0))
    ///     .WithStability(StabilityClass.D)
    ///     .OverGrid(new GeoGrid(-100, 100, -100, 100, 0, 0, 50))
    ///     .OverTime(0, 300, 10)
    ///     .WithAgent(new DQN { HiddenLayers = new[] { 64, 64 } })
    ///     .WithPolicy(new EpsilonGreedy(seed: 42))
    ///     .WithEpisodes(500)
    ///     .Run();
    /// </code>
    /// </para>
    /// </summary>
    public static class ScenarioRLAnalyzer
    {
        /// <summary>
        /// Start building an RL scenario analysis for a Gaussian plume.
        /// This convenience overload creates a <see cref="PlumeEnvironment"/> internally.
        /// </summary>
        /// <param name="emissionRate">Baseline emission rate Q (kg/s).</param>
        /// <param name="sourcePosition">Source position in world coordinates (z = stack height).</param>
        public static ScenarioRLBuilder For(double emissionRate, Vector sourcePosition)
            => new ScenarioRLBuilder(emissionRate, sourcePosition);

        /// <summary>
        /// Start building an RL scenario analysis for any GIS environment.
        /// Use this overload to plug in custom environment implementations
        /// (e.g. flood, fire, pollution models) that implement <see cref="IGISEnvironment"/>.
        /// </summary>
        /// <param name="environment">A pre-configured GIS environment.</param>
        public static ScenarioRLBuilder For(IGISEnvironment environment)
            => new ScenarioRLBuilder(environment);

        /// <summary>
        /// Start building an RL scenario analysis for any <see cref="IEnvironment"/>.
        /// Use this overload for environments that don't implement <see cref="IGISEnvironment"/>
        /// but should still benefit from the fluent RL configuration.
        /// </summary>
        /// <param name="environment">Any RL environment.</param>
        public static ScenarioRLBuilder For(IEnvironment environment)
            => new ScenarioRLBuilder(environment);
    }

    /// <summary>
    /// Builder that configures a GIS RL environment (plume or custom),
    /// wires it into <see cref="RLExperiment"/>, and returns an
    /// <see cref="RLScenarioResult"/>.
    /// </summary>
    public class ScenarioRLBuilder
    {
        // ── Pre-built environment (custom path) ──────────────────
        private readonly IEnvironment _prebuiltEnv;

        // ── Plume-specific config (convenience path) ─────────────
        private readonly double _emissionRate;
        private readonly Vector _sourcePosition;
        private double _windSpeed = 5;
        private Vector _windDirection = new Vector(1, 0, 0);
        private StabilityClass _stability = StabilityClass.D;
        private GeoGrid _grid;
        private TimeFrame _timeFrame;

        // ── Environment tuning (plume path) ──────────────────────
        private double _threshold = 1e-6;
        private double _actionCost = 0.05;
        private double _barrierEfficiency = 0.4;
        private double _filterEfficiency = 0.5;

        // ── RL config ────────────────────────────────────────────
        private IAgent _agent;
        private IPolicy _policy;
        private Func<IPolicy> _policyFactory;
        private int _episodes = 500;
        private int _maxStepsPerEpisode = 200;
        private int? _replayCapacity;
        private int? _replaySeed;
        private int? _seed;
        private int _evalEpisodes;
        private int _evalInterval;
        private RLPipelineGrid _rlGrid;
        private int _gridEvalEpisodes = 10;

        /// <summary>Plume convenience constructor.</summary>
        internal ScenarioRLBuilder(double emissionRate, Vector sourcePosition)
        {
            _emissionRate = emissionRate;
            _sourcePosition = sourcePosition;
        }

        /// <summary>Pre-built environment constructor.</summary>
        internal ScenarioRLBuilder(IEnvironment environment)
        {
            _prebuiltEnv = environment ?? throw new ArgumentNullException(nameof(environment));
        }

        // ═══════════════════════════════════════════════════════════
        //  Physics configuration
        // ═══════════════════════════════════════════════════════════

        /// <summary>Set baseline wind speed and direction (plume path only).</summary>
        public ScenarioRLBuilder WithWind(double speed, Vector direction)
        {
            EnsurePlumePath(nameof(WithWind));
            _windSpeed = speed;
            _windDirection = direction;
            return this;
        }

        /// <summary>Set atmospheric stability class (plume path only).</summary>
        public ScenarioRLBuilder WithStability(StabilityClass stability)
        {
            EnsurePlumePath(nameof(WithStability));
            _stability = stability;
            return this;
        }

        /// <summary>Define the spatial grid (plume path only).</summary>
        public ScenarioRLBuilder OverGrid(GeoGrid grid)
        {
            EnsurePlumePath(nameof(OverGrid));
            _grid = grid ?? throw new ArgumentNullException(nameof(grid));
            return this;
        }

        /// <summary>Define the time range (plume path only).</summary>
        public ScenarioRLBuilder OverTime(double startSeconds, double endSeconds, double stepSeconds)
        {
            EnsurePlumePath(nameof(OverTime));
            _timeFrame = new TimeFrame(startSeconds, endSeconds, stepSeconds);
            return this;
        }

        // ═══════════════════════════════════════════════════════════
        //  Environment tuning
        // ═══════════════════════════════════════════════════════════

        /// <summary>Set the concentration threshold (plume path only).</summary>
        public ScenarioRLBuilder WithThreshold(double threshold)
        {
            EnsurePlumePath(nameof(WithThreshold));
            _threshold = threshold;
            return this;
        }

        /// <summary>Set the per-action cost (plume path only).</summary>
        public ScenarioRLBuilder WithActionCost(double cost)
        {
            EnsurePlumePath(nameof(WithActionCost));
            _actionCost = cost;
            return this;
        }

        /// <summary>Configure barrier and filter efficiencies (plume path only).</summary>
        public ScenarioRLBuilder WithMitigationEfficiency(
            double barrierEfficiency = 0.4,
            double filterEfficiency = 0.5)
        {
            EnsurePlumePath(nameof(WithMitigationEfficiency));
            _barrierEfficiency = barrierEfficiency;
            _filterEfficiency = filterEfficiency;
            return this;
        }

        // ═══════════════════════════════════════════════════════════
        //  RL configuration
        // ═══════════════════════════════════════════════════════════

        /// <summary>Set the RL agent (e.g. DQN, PPO).</summary>
        public ScenarioRLBuilder WithAgent(IAgent agent)
        {
            _agent = agent ?? throw new ArgumentNullException(nameof(agent));
            return this;
        }

        /// <summary>Set the exploration policy.</summary>
        public ScenarioRLBuilder WithPolicy(IPolicy policy)
        {
            _policy = policy;
            return this;
        }

        /// <summary>Set a policy factory (required for grid search).</summary>
        public ScenarioRLBuilder WithPolicyFactory(Func<IPolicy> factory)
        {
            _policyFactory = factory;
            return this;
        }

        /// <summary>Set number of training episodes and max steps per episode.</summary>
        public ScenarioRLBuilder WithEpisodes(int episodes, int maxStepsPerEpisode = 200)
        {
            _episodes = episodes;
            _maxStepsPerEpisode = maxStepsPerEpisode;
            return this;
        }

        /// <summary>Add a replay buffer (for DQN and similar).</summary>
        public ScenarioRLBuilder WithReplayBuffer(int capacity, int? seed = null)
        {
            _replayCapacity = capacity;
            _replaySeed = seed;
            return this;
        }

        /// <summary>Add periodic evaluation during training.</summary>
        public ScenarioRLBuilder WithEvaluation(int evalEpisodes = 10, int evalInterval = 50)
        {
            _evalEpisodes = evalEpisodes;
            _evalInterval = evalInterval;
            return this;
        }

        /// <summary>Set the random seed for reproducibility.</summary>
        public ScenarioRLBuilder WithSeed(int seed)
        {
            _seed = seed;
            return this;
        }

        /// <summary>
        /// Use a hyperparameter grid for multi-agent comparison.
        /// Call <see cref="RunGrid"/> instead of <see cref="Run"/>.
        /// </summary>
        public ScenarioRLBuilder WithGrid(RLPipelineGrid grid, int evalEpisodes = 10)
        {
            _rlGrid = grid;
            _gridEvalEpisodes = evalEpisodes;
            return this;
        }

        // ═══════════════════════════════════════════════════════════
        //  Execute
        // ═══════════════════════════════════════════════════════════

        /// <summary>
        /// Train a single agent on the plume environment and return the result.
        /// </summary>
        public RLScenarioResult Run()
        {
            var env = BuildEnvironment();
            var builder = RLExperiment.For(env);

            if (_agent != null) builder = builder.WithAgent(_agent);
            if (_policy != null) builder = builder.WithPolicy(_policy);
            if (_replayCapacity.HasValue) builder = builder.WithReplayBuffer(_replayCapacity.Value, _replaySeed);
            if (_seed.HasValue) builder = builder.WithSeed(_seed.Value);
            if (_evalEpisodes > 0) builder = builder.WithEvaluation(_evalEpisodes, _evalInterval);

            builder = builder.WithEpisodes(_episodes, _maxStepsPerEpisode);

            var rlResult = builder.Run();

            return new RLScenarioResult(rlResult, env);
        }

        /// <summary>
        /// Run grid search over multiple agents/hyperparameters.
        /// Requires <see cref="WithGrid"/> to be called first.
        /// </summary>
        public RLGridScenarioResult RunGrid()
        {
            if (_rlGrid == null)
                throw new InvalidOperationException("Call WithGrid() before RunGrid().");

            var env = BuildEnvironment();
            var builder = RLExperiment.For(env)
                .WithGrid(_rlGrid, _gridEvalEpisodes)
                .WithEpisodes(_episodes, _maxStepsPerEpisode);

            if (_policyFactory != null) builder = builder.WithPolicyFactory(_policyFactory);
            else if (_policy != null) builder = builder.WithPolicyFactory(() => _policy.Clone());
            if (_replayCapacity.HasValue) builder = builder.WithReplayBuffer(_replayCapacity.Value, _replaySeed);
            if (_seed.HasValue) builder = builder.WithSeed(_seed.Value);

            var gridResult = builder.RunGrid();

            return new RLGridScenarioResult(gridResult, env);
        }

        private IEnvironment BuildEnvironment()
        {
            if (_prebuiltEnv != null)
                return _prebuiltEnv;

            if (_grid == null) throw new InvalidOperationException("Grid not set. Call OverGrid().");
            if (_timeFrame == null) throw new InvalidOperationException("Time frame not set. Call OverTime().");

            return new PlumeEnvironment(
                _emissionRate, _windSpeed, _windDirection,
                _sourcePosition.z, _sourcePosition,
                _grid, _timeFrame, _stability)
            {
                Threshold = _threshold,
                ActionCost = _actionCost,
                BarrierEfficiency = _barrierEfficiency,
                FilterEfficiency = _filterEfficiency
            };
        }

        private void EnsurePlumePath(string methodName)
        {
            if (_prebuiltEnv != null)
                throw new InvalidOperationException(
                    $"{methodName}() is not supported when using a pre-built environment. " +
                    "Configure the environment before passing it to ScenarioRLAnalyzer.For().");
        }
    }

    /// <summary>
    /// Result of training an RL agent on a GIS/RL scenario.
    /// Wraps <see cref="RLExperimentResult"/> with the scenario environment.
    /// </summary>
    public class RLScenarioResult
    {
        /// <summary>The RL training result.</summary>
        public RLExperimentResult ExperimentResult { get; }

        /// <summary>The environment that was used.</summary>
        public IEnvironment Environment { get; }

        /// <summary>Name of the trained agent.</summary>
        public string AgentName => ExperimentResult.AgentName;

        /// <summary>Average return over all training episodes.</summary>
        public double AverageReturn => ExperimentResult.AverageReturn;

        /// <summary>Best single-episode return.</summary>
        public double BestReturn => ExperimentResult.BestReturn;

        /// <summary>Wall-clock training time.</summary>
        public TimeSpan Duration => ExperimentResult.Duration;

        internal RLScenarioResult(RLExperimentResult rlResult, IEnvironment env)
        {
            ExperimentResult = rlResult;
            Environment = env;
        }
    }

    /// <summary>
    /// Result of RL grid search on a GIS/RL scenario.
    /// </summary>
    public class RLGridScenarioResult
    {
        /// <summary>The full grid search result.</summary>
        public RLGridSearchResult GridResult { get; }

        /// <summary>The environment that was used.</summary>
        public IEnvironment Environment { get; }

        /// <summary>Best evaluation score.</summary>
        public double BestScore => GridResult.BestScore;

        /// <summary>Best agent name.</summary>
        public string BestAgentName => GridResult.BestAgentName;

        internal RLGridScenarioResult(RLGridSearchResult gridResult, IEnvironment env)
        {
            GridResult = gridResult;
            Environment = env;
        }
    }
}
