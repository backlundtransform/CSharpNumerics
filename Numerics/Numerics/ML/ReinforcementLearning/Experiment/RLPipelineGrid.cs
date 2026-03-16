using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.ReinforcementLearning.Experiment;

/// <summary>
/// Hyperparameter grid for RL agents.
/// Mirrors <see cref="CSharpNumerics.ML.PipelineGrid"/> but for RL-specific agent configurations.
///
/// <code>
/// var grid = new RLPipelineGrid()
///     .AddAgent&lt;DQN&gt;(g => g
///         .Add("LearningRate", 0.001, 0.0005)
///         .Add("Gamma", 0.95, 0.99))
///     .AddAgent&lt;REINFORCE&gt;(g => g
///         .Add("LearningRate", 0.01, 0.005));
/// </code>
/// </summary>
public class RLPipelineGrid
{
    private readonly List<AgentGrid> _agents = new();

    /// <summary>
    /// Add an agent type with hyperparameter grid.
    /// </summary>
    public RLPipelineGrid AddAgent<T>(Action<AgentGridBuilder> setup) where T : IAgent, new()
    {
        var ag = new AgentGrid
        {
            Factory = () => new T(),
            AgentTypeName = typeof(T).Name
        };
        var builder = new AgentGridBuilder(ag);
        setup(builder);
        _agents.Add(ag);
        return this;
    }

    /// <summary>
    /// Add an agent type with a custom factory and hyperparameter grid.
    /// Useful for agents that require constructor parameters.
    /// </summary>
    public RLPipelineGrid AddAgent<T>(Func<T> factory, Action<AgentGridBuilder> setup) where T : IAgent
    {
        var ag = new AgentGrid
        {
            Factory = () => factory(),
            AgentTypeName = typeof(T).Name
        };
        var builder = new AgentGridBuilder(ag);
        setup(builder);
        _agents.Add(ag);
        return this;
    }

    /// <summary>
    /// Expand the grid into all agent configurations (cartesian product of hyperparameters).
    /// </summary>
    public List<RLAgentConfig> Expand()
    {
        var configs = new List<RLAgentConfig>();

        foreach (var ag in _agents)
        {
            var paramCombinations = CartesianProduct(ag.Parameters);
            foreach (var paramSet in paramCombinations)
            {
                configs.Add(new RLAgentConfig
                {
                    Factory = ag.Factory,
                    AgentTypeName = ag.AgentTypeName,
                    Parameters = paramSet
                });
            }
        }

        return configs;
    }

    private static List<Dictionary<string, object>> CartesianProduct(
        Dictionary<string, object[]> grid)
    {
        var result = new List<Dictionary<string, object>>
        {
            new Dictionary<string, object>()
        };

        foreach (var kv in grid)
        {
            var temp = new List<Dictionary<string, object>>();
            foreach (var value in kv.Value)
            {
                foreach (var dict in result)
                {
                    var newDict = new Dictionary<string, object>(dict)
                    {
                        [kv.Key] = value
                    };
                    temp.Add(newDict);
                }
            }
            result = temp;
        }

        return result;
    }

    // ── Inner types ─────────────────────────────────────────────

    internal class AgentGrid
    {
        public Func<IAgent> Factory { get; set; }
        public string AgentTypeName { get; set; }
        public Dictionary<string, object[]> Parameters { get; } = new();
    }

    /// <summary>Builder for specifying hyperparameter values for an agent type.</summary>
    public class AgentGridBuilder
    {
        private readonly AgentGrid _grid;

        internal AgentGridBuilder(AgentGrid grid)
        {
            _grid = grid;
        }

        /// <summary>Add a hyperparameter with one or more values to search.</summary>
        public AgentGridBuilder Add(string name, params object[] values)
        {
            _grid.Parameters[name] = values;
            return this;
        }
    }
}

/// <summary>
/// A single agent configuration from the grid expansion.
/// </summary>
public class RLAgentConfig
{
    public Func<IAgent> Factory { get; set; }
    public string AgentTypeName { get; set; }
    public Dictionary<string, object> Parameters { get; set; } = new();

    /// <summary>Create and configure the agent with the specified parameters.</summary>
    public IAgent CreateAgent()
    {
        var agent = Factory();
        agent.SetHyperParameters(Parameters);
        return agent;
    }

    /// <summary>Human-readable description of this configuration.</summary>
    public string Description
    {
        get
        {
            if (Parameters.Count == 0) return AgentTypeName;
            var parts = Parameters.Select(kv => $"{kv.Key}={kv.Value}");
            return $"{AgentTypeName} ({string.Join(", ", parts)})";
        }
    }
}
