using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Numerics.FiniteDifference;
using System;

namespace CSharpNumerics.Engines.Multiphysics.Snapshots;

/// <summary>
/// Immutable snapshot of a 2D scalar field at a single point in time.
/// Used for HeatPlate (temperature) and ElectricField (potential) simulations.
/// </summary>
public class FieldSnapshot
{
    private readonly double[] _flat;

    /// <summary>Simulation type that produced this snapshot.</summary>
    public MultiphysicsType Type { get; }

    /// <summary>Number of cells in x.</summary>
    public int Nx { get; }

    /// <summary>Number of cells in y.</summary>
    public int Ny { get; }

    /// <summary>Cell spacing in x (m).</summary>
    public double Dx { get; }

    /// <summary>Cell spacing in y (m).</summary>
    public double Dy { get; }

    /// <summary>Simulation time of this snapshot (s).</summary>
    public double Time { get; }

    /// <summary>Zero-based step index.</summary>
    public int StepIndex { get; }

    /// <summary>Total number of cells.</summary>
    public int Count => _flat.Length;

    /// <summary>
    /// Creates a snapshot from a 2D field array [ix, iy].
    /// </summary>
    public FieldSnapshot(double[,] field, MultiphysicsType type, double time, int stepIndex,
        double dx, double dy)
    {
        Nx = field.GetLength(0);
        Ny = field.GetLength(1);
        Dx = dx;
        Dy = dy;
        Type = type;
        Time = time;
        StepIndex = stepIndex;

        _flat = new double[Nx * Ny];
        for (int iy = 0; iy < Ny; iy++)
            for (int ix = 0; ix < Nx; ix++)
                _flat[iy * Nx + ix] = field[ix, iy];
    }

    /// <summary>Access value by flat row-major index.</summary>
    public double this[int flatIndex] => _flat[flatIndex];

    /// <summary>Access value by 2D index (ix, iy).</summary>
    public double this[int ix, int iy] => _flat[iy * Nx + ix];

    /// <summary>Returns a copy of the flat value array.</summary>
    public double[] GetValues() => (double[])_flat.Clone();

    /// <summary>Reconstructs the 2D array [ix, iy].</summary>
    public double[,] ToArray()
    {
        var arr = new double[Nx, Ny];
        for (int iy = 0; iy < Ny; iy++)
            for (int ix = 0; ix < Nx; ix++)
                arr[ix, iy] = _flat[iy * Nx + ix];
        return arr;
    }

    /// <summary>Maximum value in the field.</summary>
    public double Max()
    {
        double max = double.MinValue;
        for (int i = 0; i < _flat.Length; i++)
            if (_flat[i] > max) max = _flat[i];
        return max;
    }

    /// <summary>Minimum value in the field.</summary>
    public double Min()
    {
        double min = double.MaxValue;
        for (int i = 0; i < _flat.Length; i++)
            if (_flat[i] < min) min = _flat[i];
        return min;
    }
}
