using CSharpNumerics.ML.ReinforcementLearning.Algorithms.PolicyGradient;
using CSharpNumerics.ML.ReinforcementLearning.Algorithms.ValueBased;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Data;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.ReinforcementLearning.Diagnostics;

/// <summary>
/// Generates value-function surface data for continuous-state environments.
/// Samples V(s) or max-Q(s) over a grid of state-space points.
/// Returns <see cref="Serie"/>-based data ready for the existing export pipeline.
/// </summary>
public static class ValueFunctionSurface
{
    /// <summary>
    /// Sample the value function over a 1D state range.
    /// </summary>
    public static List<Serie> Sample1D(
        Func<VectorN, double> valueFunction,
        double min, double max, int numPoints = 100)
    {
        var result = new List<Serie>(numPoints);
        double step = (max - min) / (numPoints - 1);

        for (int i = 0; i < numPoints; i++)
        {
            double x = min + i * step;
            var state = new VectorN(new[] { x });
            double value = valueFunction(state);
            result.Add(new Serie { Index = x, Value = value });
        }

        return result;
    }

    /// <summary>
    /// Sample the value function over a 2D state grid.
    /// Returns a flat list of Serie where Index encodes (row * numCols + col)
    /// and Value is V(s).
    /// </summary>
    public static ValueSurface2D Sample2D(
        Func<VectorN, double> valueFunction,
        double minX, double maxX, int numX,
        double minY, double maxY, int numY)
    {
        var points = new List<Serie>(numX * numY);
        double stepX = (maxX - minX) / (numX - 1);
        double stepY = (maxY - minY) / (numY - 1);

        for (int ix = 0; ix < numX; ix++)
        {
            for (int iy = 0; iy < numY; iy++)
            {
                double x = minX + ix * stepX;
                double y = minY + iy * stepY;
                var state = new VectorN(new[] { x, y });
                double value = valueFunction(state);
                points.Add(new Serie { Index = ix * numY + iy, Value = value });
            }
        }

        return new ValueSurface2D
        {
            Points = points,
            NumX = numX,
            NumY = numY,
            MinX = minX,
            MaxX = maxX,
            MinY = minY,
            MaxY = maxY
        };
    }

    /// <summary>
    /// Extract max-Q(s) as the value function from a DQN agent.
    /// </summary>
    public static Func<VectorN, double> MaxQFunction(DQN agent) =>
        state =>
        {
            var q = agent.GetQValues(state);
            double max = q[0];
            for (int i = 1; i < q.Length; i++)
                if (q[i] > max) max = q[i];
            return max;
        };

    /// <summary>
    /// Extract max-Q(s) from a DuelingDQN agent.
    /// </summary>
    public static Func<VectorN, double> MaxQFunction(DuelingDQN agent) =>
        state =>
        {
            var q = agent.GetQValues(state);
            double max = q[0];
            for (int i = 1; i < q.Length; i++)
                if (q[i] > max) max = q[i];
            return max;
        };

    /// <summary>
    /// Extract V(s) from an ActorCritic agent.
    /// </summary>
    public static Func<VectorN, double> ValueFunction(ActorCritic agent) =>
        state => agent.GetValue(state);

    /// <summary>
    /// Extract V(s) from a PPO agent.
    /// </summary>
    public static Func<VectorN, double> ValueFunction(PPO agent) =>
        state => agent.GetValue(state);
}

/// <summary>
/// 2D value-function surface data with grid metadata.
/// </summary>
public class ValueSurface2D
{
    /// <summary>Flat list of (gridIndex, value) points.</summary>
    public List<Serie> Points { get; set; }

    public int NumX { get; set; }
    public int NumY { get; set; }
    public double MinX { get; set; }
    public double MaxX { get; set; }
    public double MinY { get; set; }
    public double MaxY { get; set; }

    /// <summary>Get the value at grid position (ix, iy).</summary>
    public double ValueAt(int ix, int iy) => Points[ix * NumY + iy].Value;

    /// <summary>Convert to a Matrix (rows=X, cols=Y) for heatmap rendering.</summary>
    public Matrix ToMatrix()
    {
        var data = new double[NumX, NumY];
        for (int ix = 0; ix < NumX; ix++)
            for (int iy = 0; iy < NumY; iy++)
                data[ix, iy] = Points[ix * NumY + iy].Value;
        return new Matrix(data);
    }
}
