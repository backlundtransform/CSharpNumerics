using CSharpNumerics.Engines.Game.AI;
using CSharpNumerics.Numerics.Objects;
using System;
using static CSharpNumerics.Engines.Game.Unity.UnityAdapter;

namespace CSharpNumerics.Engines.Game.Unity;

/// <summary>
/// Bridges Unity game state to a <see cref="GameAIAgent"/>.
///
/// Each frame:
///   1. Unity code calls <see cref="SetObservation"/> with game state
///   2. <see cref="GetAction"/> queries the AI agent and returns the action
///   3. Unity code applies the action
///
/// The bridge handles coordinate conversion and observation construction.
/// No Unity dependency — uses surrogate types.
/// </summary>
public class AIBridge
{
    private readonly GameAIAgent _agent;
    private VectorN _lastObservation;
    private VectorN _lastAction;

    /// <summary>The wrapped AI agent.</summary>
    public GameAIAgent Agent => _agent;

    /// <summary>Last observation passed to the agent.</summary>
    public VectorN LastObservation => _lastObservation;

    /// <summary>Last action returned by the agent.</summary>
    public VectorN LastAction => _lastAction;

    /// <summary>Whether the agent is enabled (if false, GetAction returns zero).</summary>
    public bool Enabled { get; set; } = true;

    /// <summary>
    /// Creates an AI bridge for a game AI agent.
    /// </summary>
    public AIBridge(GameAIAgent agent)
    {
        _agent = agent ?? throw new ArgumentNullException(nameof(agent));
    }

    /// <summary>
    /// Set the observation vector directly (raw doubles).
    /// </summary>
    public void SetObservation(VectorN observation)
    {
        _lastObservation = observation;
    }

    /// <summary>
    /// Build an observation vector from Unity-style game state.
    /// Converts positions from Unity Y-up to internal Z-up.
    /// </summary>
    /// <param name="agentPosition">Agent position in Unity coordinates.</param>
    /// <param name="agentVelocity">Agent velocity in Unity coordinates.</param>
    /// <param name="targetPosition">Target position in Unity coordinates.</param>
    /// <param name="extraFeatures">Additional observation features (health, ammo, etc.).</param>
    public void SetObservationFromUnity(
        UnityVector3 agentPosition, UnityVector3 agentVelocity,
        UnityVector3 targetPosition, double[] extraFeatures = null)
    {
        var aPos = FromUnityVector3(agentPosition);
        var aVel = FromUnityVector3(agentVelocity);
        var tPos = FromUnityVector3(targetPosition);

        int extraLen = extraFeatures?.Length ?? 0;
        var obs = new double[9 + extraLen];

        // Agent state
        obs[0] = aPos.x; obs[1] = aPos.y; obs[2] = aPos.z;
        obs[3] = aVel.x; obs[4] = aVel.y; obs[5] = aVel.z;

        // Target relative position
        obs[6] = tPos.x - aPos.x;
        obs[7] = tPos.y - aPos.y;
        obs[8] = tPos.z - aPos.z;

        // Extra features
        if (extraFeatures != null)
            Array.Copy(extraFeatures, 0, obs, 9, extraLen);

        _lastObservation = new VectorN(obs);
    }

    /// <summary>
    /// Query the AI agent for an action based on the last observation.
    /// </summary>
    /// <returns>Action vector from the agent, or zero vector if disabled.</returns>
    public VectorN GetAction()
    {
        if (!Enabled || _lastObservation.Length == 0)
        {
            _lastAction = default;
            return _lastAction;
        }

        _lastAction = _agent.Act(_lastObservation);
        return _lastAction;
    }

    /// <summary>
    /// Get the action as a Unity Vector3 (for movement commands).
    /// Takes the first 3 action components.
    /// </summary>
    public UnityVector3 GetActionAsUnityVector3()
    {
        var action = GetAction();
        if (action.Length == 0) return new UnityVector3(0, 0, 0);

        // Convert from internal Z-up to Unity Y-up
        var internal3d = new Vector(
            action.Length > 0 ? action[0] : 0,
            action.Length > 1 ? action[1] : 0,
            action.Length > 2 ? action[2] : 0);
        return ToUnityVector3(internal3d);
    }

    /// <summary>
    /// Get a discrete action index from the agent.
    /// </summary>
    public int GetDiscreteAction()
    {
        if (!Enabled || _lastObservation.Length == 0) return 0;
        return _agent.ActDiscrete(_lastObservation);
    }

    /// <summary>
    /// Reset the agent (e.g. at start of a new game round).
    /// </summary>
    public void Reset()
    {
        _agent.Reset();
        _lastObservation = default;
        _lastAction = default;
    }
}
