using CSharpNumerics.Engines.Multiphysics.Enums;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Multiphysics.Snapshots;

/// <summary>
/// Ordered collection of <see cref="FieldSnapshot"/> objects representing a
/// time-evolving simulation. Supports indexing and linear interpolation between frames.
/// </summary>
public class SimulationTimeline
{
    private readonly List<FieldSnapshot> _snapshots;

    /// <summary>Simulation type.</summary>
    public MultiphysicsType Type { get; }

    /// <summary>Number of snapshots in the timeline.</summary>
    public int Count => _snapshots.Count;

    /// <summary>Time step between snapshots (s). Zero for steady-state.</summary>
    public double Dt { get; }

    /// <summary>Time of the first snapshot.</summary>
    public double StartTime => _snapshots.Count > 0 ? _snapshots[0].Time : 0;

    /// <summary>Time of the last snapshot.</summary>
    public double EndTime => _snapshots.Count > 0 ? _snapshots[_snapshots.Count - 1].Time : 0;

    /// <summary>Access a snapshot by index.</summary>
    public FieldSnapshot this[int index] => _snapshots[index];

    public SimulationTimeline(MultiphysicsType type, double dt)
    {
        Type = type;
        Dt = dt;
        _snapshots = new List<FieldSnapshot>();
    }

    /// <summary>Add a snapshot to the end of the timeline.</summary>
    internal void Add(FieldSnapshot snapshot) => _snapshots.Add(snapshot);

    /// <summary>
    /// Build a timeline from a <see cref="SimulationResult"/> that contains
    /// a <see cref="SimulationResult.Timeline"/>.
    /// </summary>
    public static SimulationTimeline FromResult(SimulationResult result, double dt, double dx, double dy)
    {
        if (result.Timeline == null || result.Timeline.Count == 0)
            throw new InvalidOperationException("Result contains no timeline data.");

        var timeline = new SimulationTimeline(result.Type, dt);
        for (int i = 0; i < result.Timeline.Count; i++)
        {
            var snap = new FieldSnapshot(
                result.Timeline[i], result.Type,
                time: i * dt, stepIndex: i,
                dx: dx, dy: dy);
            timeline.Add(snap);
        }
        return timeline;
    }

    /// <summary>
    /// Returns the snapshot list as an array.
    /// </summary>
    public FieldSnapshot[] ToArray() => _snapshots.ToArray();

    /// <summary>
    /// Linearly interpolate the field at an arbitrary time between snapshots.
    /// Returns the nearest snapshot field if time is outside range.
    /// </summary>
    /// <param name="time">Query time in seconds.</param>
    /// <returns>Interpolated 2D field [ix, iy].</returns>
    public double[,] InterpolateAt(double time)
    {
        if (_snapshots.Count == 0)
            throw new InvalidOperationException("Timeline is empty.");

        if (_snapshots.Count == 1 || time <= _snapshots[0].Time)
            return _snapshots[0].ToArray();

        if (time >= _snapshots[_snapshots.Count - 1].Time)
            return _snapshots[_snapshots.Count - 1].ToArray();

        // Find bounding snapshots
        int lo = 0;
        for (int i = 1; i < _snapshots.Count; i++)
        {
            if (_snapshots[i].Time >= time) { lo = i - 1; break; }
        }

        var a = _snapshots[lo];
        var b = _snapshots[lo + 1];
        double alpha = (time - a.Time) / (b.Time - a.Time);

        int nx = a.Nx, ny = a.Ny;
        var result = new double[nx, ny];
        for (int iy = 0; iy < ny; iy++)
            for (int ix = 0; ix < nx; ix++)
                result[ix, iy] = a[ix, iy] * (1.0 - alpha) + b[ix, iy] * alpha;

        return result;
    }
}
