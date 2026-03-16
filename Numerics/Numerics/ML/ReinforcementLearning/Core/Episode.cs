using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.ReinforcementLearning.Core;

/// <summary>
/// A complete episode: list of transitions and computed statistics.
/// </summary>
public class Episode
{
    public List<Transition> Transitions { get; } = new();

    /// <summary>Undiscounted sum of rewards.</summary>
    public double TotalReturn => Transitions.Sum(t => t.Reward);

    /// <summary>Number of steps in the episode.</summary>
    public int Length => Transitions.Count;

    /// <summary>
    /// Compute discounted returns G_t for each step (from end to start).
    /// G_t = r_t + γ·G_{t+1}
    /// </summary>
    public double[] DiscountedReturns(double gamma)
    {
        var returns = new double[Transitions.Count];
        double g = 0;
        for (int t = Transitions.Count - 1; t >= 0; t--)
        {
            g = Transitions[t].Reward + gamma * g;
            returns[t] = g;
        }
        return returns;
    }
}
