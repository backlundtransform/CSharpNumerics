using CSharpNumerics.ML.ReinforcementLearning.Algorithms.Tabular;
using CSharpNumerics.ML.ReinforcementLearning.Environments;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Data;
using System.Collections.Generic;

namespace CSharpNumerics.ML.ReinforcementLearning.Diagnostics;

/// <summary>
/// Generates Q-value heatmap data for tabular agents on grid environments.
/// Returns <see cref="Serie"/>-based data ready for the existing export pipeline.
/// </summary>
public static class QValueHeatmap
{
    /// <summary>
    /// Get the max Q-value per state cell as a flat list of Serie (index=stateIndex, value=maxQ).
    /// </summary>
    public static List<Serie> GetMaxQValues(TabularAgent agent, GridWorld env)
    {
        var result = new List<Serie>();
        for (int r = 0; r < env.Rows; r++)
        {
            for (int c = 0; c < env.Cols; c++)
            {
                var state = new VectorN(new double[] { r, c });
                var qValues = agent.GetQValues(state);
                double maxQ = double.NegativeInfinity;
                for (int a = 0; a < qValues.Length; a++)
                    if (qValues[a] > maxQ) maxQ = qValues[a];

                result.Add(new Serie
                {
                    Index = env.StateToIndex(state),
                    Value = maxQ
                });
            }
        }
        return result;
    }

    /// <summary>
    /// Get Q-values for a specific action across all states.
    /// </summary>
    public static List<Serie> GetQValuesForAction(TabularAgent agent, GridWorld env, int action)
    {
        var result = new List<Serie>();
        for (int r = 0; r < env.Rows; r++)
        {
            for (int c = 0; c < env.Cols; c++)
            {
                var state = new VectorN(new double[] { r, c });
                var qValues = agent.GetQValues(state);
                result.Add(new Serie
                {
                    Index = env.StateToIndex(state),
                    Value = qValues[action]
                });
            }
        }
        return result;
    }

    /// <summary>
    /// Get the full Q-table as a Matrix (rows=states, columns=actions).
    /// </summary>
    public static Matrix GetQTableMatrix(TabularAgent agent) => agent.GetQTable();

    /// <summary>
    /// Get the greedy action for each state cell (index=stateIndex, value=actionIndex).
    /// </summary>
    public static List<Serie> GetGreedyPolicy(TabularAgent agent, GridWorld env)
    {
        var result = new List<Serie>();
        for (int r = 0; r < env.Rows; r++)
        {
            for (int c = 0; c < env.Cols; c++)
            {
                var state = new VectorN(new double[] { r, c });
                int action = agent.SelectAction(state);
                result.Add(new Serie
                {
                    Index = env.StateToIndex(state),
                    Value = action
                });
            }
        }
        return result;
    }
}
