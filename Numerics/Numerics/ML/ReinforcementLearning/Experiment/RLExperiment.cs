using CSharpNumerics.ML.ReinforcementLearning.Algorithms.ContinuousControl;
using CSharpNumerics.ML.ReinforcementLearning.Algorithms.PolicyGradient;
using CSharpNumerics.ML.ReinforcementLearning.Algorithms.Tabular;
using CSharpNumerics.ML.ReinforcementLearning.Algorithms.ValueBased;
using CSharpNumerics.ML.ReinforcementLearning.Buffers;
using CSharpNumerics.ML.ReinforcementLearning.Core;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.ML.ReinforcementLearning.Policies;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace CSharpNumerics.ML.ReinforcementLearning.Experiment;

/// <summary>
/// Fluent entry point for RL experiments.
/// <code>
/// var result = RLExperiment
///     .For(new GridWorld(5, 5))
///     .WithAgent(new QLearning(25, 4, env.StateToIndex))
///     .WithPolicy(new EpsilonGreedy())
///     .WithEpisodes(1000, 200)
///     .Run();
/// </code>
/// </summary>
public static class RLExperiment
{
    /// <summary>Start building an RL experiment for the given environment.</summary>
    public static RLExperimentBuilder For(IEnvironment environment)
        => new RLExperimentBuilder(environment);
}

/// <summary>
/// Builder that configures and runs an RL training loop.
/// </summary>
public class RLExperimentBuilder
{
    private readonly IEnvironment _env;
    private IAgent _agent;
    private IPolicy _policy;
    private IReplayBuffer _replayBuffer;
    private int _maxEpisodes = 1000;
    private int _maxStepsPerEpisode = 500;
    private int _evalEpisodes = 0;
    private int _evalInterval = 50;
    private int? _seed;

    // Grid search
    private RLPipelineGrid _grid;
    private Func<IPolicy> _policyFactory;
    private Func<IReplayBuffer> _bufferFactory;
    private int _gridEvalEpisodes = 10;

    // Monte Carlo evaluation
    private int _mcRuns;
    private int _mcEvalEpisodes;

    internal RLExperimentBuilder(IEnvironment environment)
    {
        _env = environment;
    }

    public RLExperimentBuilder WithAgent(IAgent agent)
    {
        _agent = agent;
        return this;
    }

    public RLExperimentBuilder WithPolicy(IPolicy policy)
    {
        _policy = policy;
        return this;
    }

    public RLExperimentBuilder WithReplayBuffer(int capacity, int? seed = null)
    {
        _replayBuffer = new ReplayBuffer(capacity, seed);
        return this;
    }

    public RLExperimentBuilder WithReplayBuffer(IReplayBuffer buffer)
    {
        if (buffer != null) _replayBuffer = buffer;
        return this;
    }

    public RLExperimentBuilder WithEpisodes(int maxEpisodes, int maxStepsPerEpisode = 500)
    {
        _maxEpisodes = maxEpisodes;
        _maxStepsPerEpisode = maxStepsPerEpisode;
        return this;
    }

    /// <summary>
    /// Configure periodic greedy evaluation during training.
    /// Every evalInterval episodes, run evalEpisodes greedy rollouts and record the average return.
    /// </summary>
    public RLExperimentBuilder WithEvaluation(int evalEpisodes, int evalInterval = 50)
    {
        _evalEpisodes = evalEpisodes;
        _evalInterval = evalInterval;
        return this;
    }

    public RLExperimentBuilder WithSeed(int seed)
    {
        _seed = seed;
        return this;
    }

    /// <summary>
    /// Configure a hyperparameter grid search.
    /// Each expanded configuration will be trained and evaluated.
    /// </summary>
    public RLExperimentBuilder WithGrid(RLPipelineGrid grid, int evalEpisodes = 10)
    {
        _grid = grid;
        _gridEvalEpisodes = evalEpisodes;
        return this;
    }

    /// <summary>
    /// Configure Monte Carlo evaluation: multiple independent training runs.
    /// </summary>
    public RLExperimentBuilder WithMonteCarloEvaluation(int runs, int evalEpisodesPerRun = 10)
    {
        _mcRuns = runs;
        _mcEvalEpisodes = evalEpisodesPerRun;
        return this;
    }

    /// <summary>
    /// Set a policy factory for grid search (fresh policy per config).
    /// </summary>
    public RLExperimentBuilder WithPolicyFactory(Func<IPolicy> policyFactory)
    {
        _policyFactory = policyFactory;
        return this;
    }

    /// <summary>
    /// Set a replay buffer factory for grid search (fresh buffer per config).
    /// </summary>
    public RLExperimentBuilder WithReplayBufferFactory(Func<IReplayBuffer> bufferFactory)
    {
        _bufferFactory = bufferFactory;
        return this;
    }

    /// <summary>Execute the training loop and return results.</summary>
    public RLExperimentResult Run()
    {
        if (_agent == null)
            throw new InvalidOperationException("Agent is required. Call WithAgent() before Run().");

        _policy ??= new EpsilonGreedy();

        // Auto-initialize deep RL agents
        InitializeDeepAgent(_agent);

        var training = new TrainingResult();
        var sw = Stopwatch.StartNew();

        for (int ep = 0; ep < _maxEpisodes; ep++)
        {
            var episode = RunEpisode(training: true);

            _agent.EndEpisode(episode);
            _policy.Decay();

            training.EpisodeReturns.Add(episode.TotalReturn);
            training.EpisodeSteps.Add(episode.Length);

            if (_policy is EpsilonGreedy eg)
                training.ExplorationValues.Add(eg.Epsilon);

            // Periodic evaluation
            if (_evalEpisodes > 0 && (ep + 1) % _evalInterval == 0)
            {
                double evalReturn = EvaluateGreedy(_evalEpisodes);
                training.EvalReturns.Add(evalReturn);
            }
        }

        sw.Stop();

        return new RLExperimentResult
        {
            Training = training,
            AgentName = _agent.Name,
            Parameters = _agent.GetHyperParameters(),
            Duration = sw.Elapsed
        };
    }

    /// <summary>
    /// Run a grid search: train each agent config and rank by evaluation return.
    /// </summary>
    public RLGridSearchResult RunGrid()
    {
        if (_grid == null)
            throw new InvalidOperationException("Grid is required. Call WithGrid() before RunGrid().");

        var configs = _grid.Expand();
        var entries = new List<RLGridSearchEntry>();
        var totalSw = Stopwatch.StartNew();

        foreach (var config in configs)
        {
            var agent = config.CreateAgent();
            var policy = _policyFactory?.Invoke() ?? _policy?.Clone() ?? new EpsilonGreedy();
            var buffer = _bufferFactory?.Invoke() ?? _replayBuffer;

            var singleResult = RLExperiment
                .For(_env)
                .WithAgent(agent)
                .WithPolicy(policy)
                .WithReplayBuffer(buffer)
                .WithEpisodes(_maxEpisodes, _maxStepsPerEpisode)
                .WithEvaluation(_evalEpisodes, _evalInterval)
                .MaybeWithSeed(_seed)
                .Run();

            var evaluator = new EpisodeEvaluator(_env, _maxStepsPerEpisode);
            var evalResult = evaluator.Evaluate(agent, _gridEvalEpisodes, _seed);

            entries.Add(new RLGridSearchEntry
            {
                AgentName = config.AgentTypeName,
                Parameters = config.Parameters,
                Description = config.Description,
                AverageEvalReturn = evalResult.MeanReturn,
                EvalResult = evalResult,
                TrainingResult = singleResult,
                Duration = singleResult.Duration
            });
        }

        totalSw.Stop();

        entries.Sort((a, b) => b.AverageEvalReturn.CompareTo(a.AverageEvalReturn));

        return new RLGridSearchResult
        {
            Rankings = entries,
            TotalDuration = totalSw.Elapsed
        };
    }

    /// <summary>
    /// Run multiple independent training runs with post-training evaluation.
    /// Useful for computing confidence intervals on performance.
    /// </summary>
    public RLMonteCarloResult RunMonteCarlo()
    {
        if (_agent == null)
            throw new InvalidOperationException("Agent is required. Call WithAgent() before RunMonteCarlo().");
        if (_mcRuns <= 0)
            throw new InvalidOperationException("Configure Monte Carlo runs via WithMonteCarloEvaluation() first.");

        var runs = new List<RLExperimentResult>();
        var evals = new List<EvalResult>();

        for (int run = 0; run < _mcRuns; run++)
        {
            var agent = _agent.Clone();
            var policy = _policyFactory?.Invoke() ?? _policy?.Clone() ?? new EpsilonGreedy();
            var buffer = _bufferFactory?.Invoke() ?? _replayBuffer;
            int runSeed = _seed.HasValue ? _seed.Value + run * 1000 : run;

            var result = RLExperiment
                .For(_env)
                .WithAgent(agent)
                .WithPolicy(policy)
                .WithReplayBuffer(buffer)
                .WithEpisodes(_maxEpisodes, _maxStepsPerEpisode)
                .WithEvaluation(_evalEpisodes, _evalInterval)
                .WithSeed(runSeed)
                .Run();

            var evaluator = new EpisodeEvaluator(_env, _maxStepsPerEpisode);
            var evalResult = evaluator.Evaluate(agent, _mcEvalEpisodes, runSeed);

            runs.Add(result);
            evals.Add(evalResult);
        }

        return new RLMonteCarloResult
        {
            Runs = runs,
            EvalResults = evals
        };
    }

    /// <summary>Internal: conditionally apply seed.</summary>
    private RLExperimentBuilder MaybeWithSeed(int? seed)
    {
        if (seed.HasValue) _seed = seed.Value;
        return this;
    }

    private Episode RunEpisode(bool training)
    {
        var episode = new Episode();
        var (state, _) = _env.Reset(_seed.HasValue ? _seed.Value + episode.GetHashCode() : null);

        // Reset OU process at the start of each episode
        if (_policy is OrnsteinUhlenbeck ou)
            ou.Reset(_env.ActionSize);

        bool continuous = !_env.IsDiscrete;

        for (int step = 0; step < _maxStepsPerEpisode; step++)
        {
            Transition transition;

            if (continuous)
            {
                var contAction = SelectContinuousAction(state, training);
                var (nextState, reward, done, _) = _env.Step(contAction);
                transition = new Transition(state, contAction, reward, nextState, done);
                episode.Transitions.Add(transition);
            }
            else
            {
                int action = SelectAction(state, training);
                var (nextState, reward, done, _) = _env.Step(action);
                transition = new Transition(state, action, reward, nextState, done);
                episode.Transitions.Add(transition);
            }

            if (training)
                _agent.Train(transition);

            state = transition.NextState;
            if (transition.Done) break;
        }

        return episode;
    }

    private int SelectAction(VectorN state, bool training)
    {
        if (training)
        {
            if (_agent is TabularAgent tabular)
                return _policy.SelectAction(tabular.GetQValues(state));
            if (_agent is DQN dqn)
                return _policy.SelectAction(dqn.GetQValues(state));
            if (_agent is DuelingDQN dueling)
                return _policy.SelectAction(dueling.GetQValues(state));
            return _agent.SelectAction(state);
        }
        else
        {
            // Greedy evaluation — no exploration
            return _agent.SelectAction(state);
        }
    }

    private VectorN SelectContinuousAction(VectorN state, bool training)
    {
        var meanAction = _agent.SelectContinuousAction(state);
        if (training && _policy != null)
            return _policy.SelectAction(meanAction, new VectorN(meanAction.Length));
        return meanAction;
    }

    private double EvaluateGreedy(int numEpisodes)
    {
        double total = 0;
        for (int i = 0; i < numEpisodes; i++)
        {
            var episode = RunEpisode(training: false);
            total += episode.TotalReturn;
        }
        return total / numEpisodes;
    }

    private void InitializeDeepAgent(IAgent agent)
    {
        if (agent is DQN dqn)
        {
            if (_replayBuffer != null) dqn.SetReplayBuffer(_replayBuffer);
            dqn.Initialize(_env.ObservationSize, _env.ActionSize, _seed);
        }
        else if (agent is DuelingDQN dueling)
        {
            if (_replayBuffer != null) dueling.SetReplayBuffer(_replayBuffer);
            dueling.Initialize(_env.ObservationSize, _env.ActionSize, _seed);
        }
        else if (agent is REINFORCE reinforce)
        {
            reinforce.Initialize(_env.ObservationSize, _env.ActionSize, _seed);
        }
        else if (agent is ActorCritic ac)
        {
            ac.Initialize(_env.ObservationSize, _env.ActionSize, _seed);
        }
        else if (agent is PPO ppo)
        {
            ppo.Initialize(_env.ObservationSize, _env.ActionSize, _seed);
        }
        else if (agent is DDPG ddpg)
        {
            if (_replayBuffer != null) ddpg.SetReplayBuffer(_replayBuffer);
            ddpg.Initialize(_env.ObservationSize, _env.ActionSize, _seed);
        }
    }
}
