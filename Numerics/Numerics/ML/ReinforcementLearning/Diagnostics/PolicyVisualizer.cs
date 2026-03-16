using CSharpNumerics.ML.ReinforcementLearning.Algorithms.PolicyGradient;
using CSharpNumerics.ML.ReinforcementLearning.Algorithms.Tabular;
using CSharpNumerics.ML.ReinforcementLearning.Environments;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Data;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.ReinforcementLearning.Diagnostics;

/// <summary>
/// Extracts policy data (action probabilities or greedy actions) for visualisation.
/// </summary>
public static class PolicyVisualizer
{
    /// <summary>
    /// Get action probabilities for each state in a GridWorld, from a policy gradient agent.
    /// Returns one List&lt;Serie&gt; per action — Serie.Index = stateIndex, Serie.Value = probability.
    /// </summary>
    public static List<List<Serie>> GetActionProbabilities(
        GridWorld env, Func<VectorN, VectorN> getProbs)
    {
        int numActions = env.ActionSize;
        var perAction = new List<List<Serie>>(numActions);
        for (int a = 0; a < numActions; a++)
            perAction.Add(new List<Serie>());

        for (int r = 0; r < env.Rows; r++)
        {
            for (int c = 0; c < env.Cols; c++)
            {
                var state = new VectorN(new double[] { r, c });
                int stateIdx = env.StateToIndex(state);
                var probs = getProbs(state);

                for (int a = 0; a < numActions; a++)
                {
                    perAction[a].Add(new Serie
                    {
                        Index = stateIdx,
                        Value = probs[a]
                    });
                }
            }
        }

        return perAction;
    }

    /// <summary>
    /// Get softmax action probabilities from a tabular agent's Q-values.
    /// Temperature controls exploration: low → greedy, high → uniform.
    /// </summary>
    public static List<List<Serie>> GetSoftmaxProbabilities(
        TabularAgent agent, GridWorld env, double temperature = 1.0)
    {
        return GetActionProbabilities(env, state =>
        {
            var qValues = agent.GetQValues(state);
            return Softmax(qValues, temperature);
        });
    }

    /// <summary>
    /// Get the entropy of the policy at each state (measures uncertainty).
    /// High entropy → uniform (exploring); low entropy → deterministic (exploiting).
    /// </summary>
    public static List<Serie> GetPolicyEntropy(
        GridWorld env, Func<VectorN, VectorN> getProbs)
    {
        var result = new List<Serie>();
        for (int r = 0; r < env.Rows; r++)
        {
            for (int c = 0; c < env.Cols; c++)
            {
                var state = new VectorN(new double[] { r, c });
                var probs = getProbs(state);

                double entropy = 0;
                for (int a = 0; a < probs.Length; a++)
                {
                    if (probs[a] > 1e-10)
                        entropy -= probs[a] * Math.Log(probs[a]);
                }

                result.Add(new Serie
                {
                    Index = env.StateToIndex(state),
                    Value = entropy
                });
            }
        }
        return result;
    }

    /// <summary>
    /// Get the dominant (most probable) action at each state.
    /// </summary>
    public static List<Serie> GetDominantAction(
        GridWorld env, Func<VectorN, VectorN> getProbs)
    {
        var result = new List<Serie>();
        for (int r = 0; r < env.Rows; r++)
        {
            for (int c = 0; c < env.Cols; c++)
            {
                var state = new VectorN(new double[] { r, c });
                var probs = getProbs(state);

                int best = 0;
                for (int a = 1; a < probs.Length; a++)
                    if (probs[a] > probs[best]) best = a;

                result.Add(new Serie
                {
                    Index = env.StateToIndex(state),
                    Value = best
                });
            }
        }
        return result;
    }

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
