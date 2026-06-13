using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Engines.Game.AI;

/// <summary>
/// Wraps a trained RL policy for use as a game AI agent.
/// Each tick, the agent receives a game-state observation vector and
/// returns an action via the underlying RL agent's policy network.
///
/// Supports both discrete and continuous action spaces.
/// </summary>
public class GameAIAgent
{
    private readonly IAgent _agent;
    private readonly bool _continuous;
    private readonly string _name;

    /// <summary>Agent name.</summary>
    public string Name => _name;

    /// <summary>Whether the agent uses continuous actions.</summary>
    public bool IsContinuous => _continuous;

    /// <summary>The underlying RL agent.</summary>
    public IAgent Agent => _agent;

    /// <summary>Number of actions taken so far.</summary>
    public int ActionCount { get; private set; }

    /// <summary>
    /// Creates a game AI agent wrapping a trained RL policy.
    /// </summary>
    /// <param name="agent">A trained RL agent (PPO, DDPG, DQN, etc.).</param>
    /// <param name="continuous">True for continuous action space.</param>
    /// <param name="name">Agent name for identification.</param>
    public GameAIAgent(IAgent agent, bool continuous = true, string name = null)
    {
        _agent = agent ?? throw new ArgumentNullException(nameof(agent));
        _continuous = continuous;
        _name = name ?? agent.Name;
    }

    /// <summary>
    /// Given a state observation, return the best action.
    /// </summary>
    /// <param name="observation">Current state observation vector.</param>
    /// <returns>Continuous action vector or single-element vector with discrete action index.</returns>
    public VectorN Act(VectorN observation)
    {
        ActionCount++;
        if (_continuous)
            return _agent.SelectContinuousAction(observation);
        else
            return new VectorN(new double[] { _agent.SelectAction(observation) });
    }

    /// <summary>
    /// Select a discrete action from the observation.
    /// </summary>
    public int ActDiscrete(VectorN observation)
    {
        ActionCount++;
        return _agent.SelectAction(observation);
    }

    /// <summary>
    /// Reset the action counter (e.g. at start of a new game round).
    /// </summary>
    public void Reset()
    {
        ActionCount = 0;
    }
}
