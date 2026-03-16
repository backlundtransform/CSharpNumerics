using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.ReinforcementLearning.Environments;

/// <summary>
/// A simple deterministic grid world.
/// The agent starts at (0,0) and must reach the goal at (Rows-1, Cols-1).
/// Actions: 0=Up, 1=Right, 2=Down, 3=Left.
/// Reward: +1 on reaching the goal, -0.01 per step (encourages shortest path).
/// Walls are impassable; attempting to move into a wall or off the grid stays in place.
/// State is a 2-element VectorN: [row, col].
/// </summary>
public class GridWorld : IEnvironment
{
    public int Rows { get; }
    public int Cols { get; }

    private readonly HashSet<(int row, int col)> _walls;
    private readonly (int row, int col) _goal;
    private int _row;
    private int _col;
    private Random _rng;

    // Directions: Up, Right, Down, Left
    private static readonly int[] _dRow = { -1, 0, 1, 0 };
    private static readonly int[] _dCol = { 0, 1, 0, -1 };

    public int ObservationSize => 2;
    public int ActionSize => 4;
    public bool IsDiscrete => true;

    /// <summary>Total number of distinct states (Rows × Cols minus walls).</summary>
    public int StateCount => Rows * Cols;

    public GridWorld(int rows = 5, int cols = 5, HashSet<(int row, int col)> walls = null,
        (int row, int col)? goal = null)
    {
        Rows = rows;
        Cols = cols;
        _walls = walls ?? new HashSet<(int, int)>();
        _goal = goal ?? (rows - 1, cols - 1);
        _rng = new Random();
    }

    public (VectorN state, Dictionary<string, object> info) Reset(int? seed = null)
    {
        _rng = seed.HasValue ? new Random(seed.Value) : new Random();
        _row = 0;
        _col = 0;
        return (CurrentState(), new Dictionary<string, object>());
    }

    public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(int action)
    {
        int newRow = _row + _dRow[action];
        int newCol = _col + _dCol[action];

        if (newRow >= 0 && newRow < Rows && newCol >= 0 && newCol < Cols
            && !_walls.Contains((newRow, newCol)))
        {
            _row = newRow;
            _col = newCol;
        }

        bool done = (_row == _goal.row && _col == _goal.col);
        double reward = done ? 1.0 : -0.01;

        return (CurrentState(), reward, done, new Dictionary<string, object>());
    }

    public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(VectorN action)
    {
        return Step((int)action[0]);
    }

    /// <summary>Flat state index = row * Cols + col. Useful for tabular methods.</summary>
    public int StateToIndex(VectorN state) => (int)state[0] * Cols + (int)state[1];

    private VectorN CurrentState() => new VectorN(new double[] { _row, _col });
}
