using CSharpNumerics.ML.ReinforcementLearning.Core;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.ML.ReinforcementLearning.Policies;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Game.AI;

/// <summary>
/// Offline training loop for game AI agents.
/// Runs an RL environment headless (no rendering) and trains an agent
/// using the standard RL training loop with configurable exploration policy.
///
/// After training, the agent can be wrapped in a <see cref="GameAIAgent"/>
/// for real-time game use.
/// </summary>
public class AITrainer
{
    private readonly IEnvironment _env;
    private readonly IAgent _agent;
    private IPolicy _policy;
    private int _maxEpisodes;
    private int _maxStepsPerEpisode;
    private int? _seed;

    /// <summary>Training results from the most recent Run().</summary>
    public TrainingResult LastResult { get; private set; }

    /// <summary>Episode returns logged during training.</summary>
    public List<double> EpisodeReturns => LastResult?.EpisodeReturns ?? new List<double>();

    /// <summary>
    /// Creates a trainer.
    /// </summary>
    /// <param name="env">The RL environment to train in.</param>
    /// <param name="agent">The RL agent to train (PPO, DDPG, etc.).</param>
    public AITrainer(IEnvironment env, IAgent agent)
    {
        _env = env ?? throw new ArgumentNullException(nameof(env));
        _agent = agent ?? throw new ArgumentNullException(nameof(agent));
        _maxEpisodes = 500;
        _maxStepsPerEpisode = 1000;
    }

    /// <summary>Set the exploration policy.</summary>
    public AITrainer WithPolicy(IPolicy policy)
    {
        _policy = policy;
        return this;
    }

    /// <summary>Set training episode counts.</summary>
    public AITrainer WithEpisodes(int maxEpisodes, int maxStepsPerEpisode = 1000)
    {
        _maxEpisodes = maxEpisodes;
        _maxStepsPerEpisode = maxStepsPerEpisode;
        return this;
    }

    /// <summary>Set random seed for reproducibility.</summary>
    public AITrainer WithSeed(int seed)
    {
        _seed = seed;
        return this;
    }

    /// <summary>
    /// Runs the training loop and returns the trained agent wrapped as a GameAIAgent.
    /// </summary>
    /// <param name="agentName">Name for the resulting GameAIAgent.</param>
    /// <returns>A GameAIAgent wrapping the trained policy.</returns>
    public GameAIAgent Train(string agentName = null)
    {
        _policy ??= _env.IsDiscrete
            ? (IPolicy)new EpsilonGreedy(_seed)
            : new GaussianNoise(_seed);

        var result = new TrainingResult();
        bool continuous = !_env.IsDiscrete;

        for (int ep = 0; ep < _maxEpisodes; ep++)
        {
            var episode = new Episode();
            int episodeSeed = _seed.HasValue ? _seed.Value + ep : ep;
            var (state, _) = _env.Reset(episodeSeed);

            if (_policy is OrnsteinUhlenbeck ou)
                ou.Reset(_env.ActionSize);

            for (int step = 0; step < _maxStepsPerEpisode; step++)
            {
                Transition transition;

                if (continuous)
                {
                    var rawAction = _agent.SelectContinuousAction(state);
                    var exploredAction = _policy.SelectAction(rawAction, rawAction);
                    var (nextState, reward, done, _) = _env.Step(exploredAction);
                    transition = new Transition(state, exploredAction, reward, nextState, done);
                }
                else
                {
                    var qOrProb = new VectorN(new double[_env.ActionSize]); // placeholder
                    int action = _policy.SelectAction(qOrProb);
                    var (nextState, reward, done, _) = _env.Step(action);
                    transition = new Transition(state, action, reward, nextState, done);
                }

                episode.Transitions.Add(transition);
                _agent.Train(transition);

                state = transition.NextState;
                if (transition.Done) break;
            }

            _agent.EndEpisode(episode);
            _policy.Decay();

            result.EpisodeReturns.Add(episode.TotalReturn);
            result.EpisodeSteps.Add(episode.Length);
        }

        LastResult = result;

        return new GameAIAgent(_agent, continuous, agentName ?? _agent.Name);
    }

    /// <summary>
    /// Evaluate the trained agent on N episodes without exploration.
    /// </summary>
    /// <param name="numEpisodes">Number of evaluation episodes.</param>
    /// <returns>Average return across evaluation episodes.</returns>
    public double Evaluate(int numEpisodes = 10)
    {
        double totalReturn = 0;
        bool continuous = !_env.IsDiscrete;

        for (int ep = 0; ep < numEpisodes; ep++)
        {
            var (state, _) = _env.Reset(_seed.HasValue ? _seed.Value + 10000 + ep : (int?)null);
            double episodeReturn = 0;

            for (int step = 0; step < _maxStepsPerEpisode; step++)
            {
                VectorN nextState;
                double reward;
                bool done;

                if (continuous)
                {
                    var action = _agent.SelectContinuousAction(state);
                    (nextState, reward, done, _) = _env.Step(action);
                }
                else
                {
                    int action = _agent.SelectAction(state);
                    (nextState, reward, done, _) = _env.Step(action);
                }

                episodeReturn += reward;
                state = nextState;
                if (done) break;
            }

            totalReturn += episodeReturn;
        }

        return totalReturn / numEpisodes;
    }
}
