using System;
using System.Collections.Generic;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Mechanics;
using CSharpNumerics.Physics.OrbitalMechanics;

namespace CSharpNumerics.Engines.Game.Rocket;

/// <summary>
/// Records telemetry data (state vectors, orbital elements, events) at a configurable sample rate.
/// </summary>
public class TelemetryRecorder
{
    private readonly List<TelemetryFrame> _frames = new List<TelemetryFrame>();
    private double _lastRecordTime = double.NegativeInfinity;

    /// <summary>Minimum interval between recorded frames (seconds). Default 1.0.</summary>
    public double SampleInterval { get; set; } = 1.0;

    /// <summary>All recorded telemetry frames.</summary>
    public IReadOnlyList<TelemetryFrame> Frames => _frames;

    /// <summary>Number of recorded frames.</summary>
    public int FrameCount => _frames.Count;

    /// <summary>
    /// Records a telemetry frame if enough time has elapsed since the last recording.
    /// </summary>
    /// <param name="state">Current rocket state.</param>
    /// <param name="time">Current simulation time (s).</param>
    /// <param name="useECIElements">If true, compute orbital elements from ECI state.</param>
    public void Record(RocketState state, double time, bool useECIElements = false)
    {
        if (time - _lastRecordTime < SampleInterval) return;

        var frame = new TelemetryFrame
        {
            Time = time,
            Position = state.Position,
            Velocity = state.Velocity,
            Mass = state.Mass,
            Altitude = state.Altitude,
            Speed = state.Speed
        };

        if (useECIElements)
        {
            frame.Elements = StateToElements.FromStateVector(state.Position, state.Velocity, EarthModel.GM);
            frame.HasElements = true;
        }

        _frames.Add(frame);
        _lastRecordTime = time;
    }

    /// <summary>Clears all recorded telemetry.</summary>
    public void Clear()
    {
        _frames.Clear();
        _lastRecordTime = double.NegativeInfinity;
    }

    /// <summary>
    /// Gets the maximum altitude reached across all recorded frames.
    /// </summary>
    public double MaxAltitude
    {
        get
        {
            double max = 0;
            for (int i = 0; i < _frames.Count; i++)
                if (_frames[i].Altitude > max)
                    max = _frames[i].Altitude;
            return max;
        }
    }

    /// <summary>
    /// Gets the maximum speed reached across all recorded frames.
    /// </summary>
    public double MaxSpeed
    {
        get
        {
            double max = 0;
            for (int i = 0; i < _frames.Count; i++)
                if (_frames[i].Speed > max)
                    max = _frames[i].Speed;
            return max;
        }
    }
}

/// <summary>
/// A single telemetry snapshot at a point in time.
/// </summary>
public struct TelemetryFrame
{
    /// <summary>Simulation time (seconds).</summary>
    public double Time;

    /// <summary>Position vector.</summary>
    public Vector Position;

    /// <summary>Velocity vector.</summary>
    public Vector Velocity;

    /// <summary>Vehicle mass (kg).</summary>
    public double Mass;

    /// <summary>Altitude (meters).</summary>
    public double Altitude;

    /// <summary>Speed (m/s).</summary>
    public double Speed;

    /// <summary>Whether orbital elements are computed for this frame.</summary>
    public bool HasElements;

    /// <summary>Orbital elements (only valid if HasElements is true).</summary>
    public OrbitalElements Elements;
}
