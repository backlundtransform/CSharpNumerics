using System;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Mechanics;
using CSharpNumerics.Physics.OrbitalMechanics;

namespace CSharpNumerics.Engines.Game.Rocket;

/// <summary>
/// Provides per-frame telemetry data for HUD display at a configurable delivery rate.
/// Buffers and interpolates between simulation steps to maintain smooth frame delivery
/// independent of simulation timestep.
///
/// Typical usage: simulation runs at variable dt, TelemetryStream delivers data at 60fps.
/// </summary>
public class TelemetryStream
{
    /// <summary>Target frame delivery rate (frames/sec). Default 60.</summary>
    public double TargetFrameRate { get; set; } = 60.0;

    /// <summary>Last delivered telemetry snapshot.</summary>
    public TelemetrySnapshot Current { get; private set; } = new TelemetrySnapshot();

    /// <summary>Number of frames delivered since last reset.</summary>
    public long FramesDelivered { get; private set; }

    /// <summary>Actual achieved delivery rate (frames/sim-second).</summary>
    public double ActualFrameRate { get; private set; }

    private TelemetrySnapshot _prev = new TelemetrySnapshot();
    private TelemetrySnapshot _next = new TelemetrySnapshot();
    private double _simTime;
    private double _lastDeliveryTime = double.NegativeInfinity;
    private double _deliveryInterval;
    private double _rateWindowStart;
    private long _rateWindowFrames;

    /// <summary>
    /// Pushes a new simulation state snapshot. Call this each simulation step.
    /// </summary>
    /// <param name="state">Current rocket state.</param>
    /// <param name="simTime">Current simulation time (s).</param>
    /// <param name="thrustMagnitude">Current thrust (N).</param>
    /// <param name="fuelFraction">Remaining fuel as fraction (0–1).</param>
    /// <param name="mu">Gravitational parameter for orbital elements.</param>
    public void Push(RocketState state, double simTime, double thrustMagnitude = 0,
        double fuelFraction = 1.0, double mu = 3.986004418e14)
    {
        _prev = _next;

        double altitude = state.Altitude;
        double speed = state.Speed;
        double r = state.Position.GetMagnitude();
        double accel = thrustMagnitude / Math.Max(state.Mass, 1.0);

        // Compute orbital elements if in space (altitude > 100 km or high speed)
        double apoAlt = 0, periAlt = 0;
        if (r > EarthModel.SemiMajorAxis + 50000 && speed > 1000)
        {
            double energy = 0.5 * speed * speed - mu / r;
            if (energy < 0)
            {
                double sma = -mu / (2.0 * energy);
                Vector h = new Vector(
                    state.Position.y * state.Velocity.z - state.Position.z * state.Velocity.y,
                    state.Position.z * state.Velocity.x - state.Position.x * state.Velocity.z,
                    state.Position.x * state.Velocity.y - state.Position.y * state.Velocity.x);
                double hMag = h.GetMagnitude();
                double ecc = Math.Sqrt(Math.Max(0, 1.0 + 2.0 * energy * hMag * hMag / (mu * mu)));
                apoAlt = sma * (1.0 + ecc) - EarthModel.SemiMajorAxis;
                periAlt = sma * (1.0 - ecc) - EarthModel.SemiMajorAxis;
            }
            else
            {
                apoAlt = double.PositiveInfinity;
                periAlt = r - EarthModel.SemiMajorAxis;
            }
        }

        _next = new TelemetrySnapshot
        {
            SimTime = simTime,
            Altitude = altitude,
            Speed = speed,
            Velocity = state.Velocity,
            Acceleration = accel,
            FuelPercent = fuelFraction * 100.0,
            ApoapsisAltitude = apoAlt,
            PeriapsisAltitude = periAlt,
            Mass = state.Mass
        };

        _simTime = simTime;
    }

    /// <summary>
    /// Attempts to deliver a telemetry frame. Returns true if a new frame is available
    /// (based on target frame rate), false if it's too soon.
    /// </summary>
    /// <param name="simTime">Current simulation time (s).</param>
    /// <returns>True if a new frame was delivered (accessible via <see cref="Current"/>).</returns>
    public bool TryDeliver(double simTime)
    {
        _deliveryInterval = 1.0 / TargetFrameRate;

        if (simTime - _lastDeliveryTime < _deliveryInterval)
            return false;

        // Interpolate between prev and next
        double t = 0;
        double span = _next.SimTime - _prev.SimTime;
        if (span > 0)
            t = Math.Min(1.0, (simTime - _prev.SimTime) / span);

        var snapshot = Interpolate(_prev, _next, t);
        snapshot.SimTime = simTime;
        Current = snapshot;
        _lastDeliveryTime = simTime;
        FramesDelivered++;

        // Update rate measurement
        _rateWindowFrames++;
        double windowDuration = simTime - _rateWindowStart;
        if (windowDuration >= 1.0)
        {
            ActualFrameRate = _rateWindowFrames / windowDuration;
            _rateWindowStart = simTime;
            _rateWindowFrames = 0;
        }

        return true;
    }

    /// <summary>Resets the stream state.</summary>
    public void Reset()
    {
        Current = new TelemetrySnapshot();
        _prev = new TelemetrySnapshot();
        _next = new TelemetrySnapshot();
        FramesDelivered = 0;
        _lastDeliveryTime = double.NegativeInfinity;
        ActualFrameRate = 0;
        _rateWindowStart = 0;
        _rateWindowFrames = 0;
    }

    private static TelemetrySnapshot Interpolate(TelemetrySnapshot a, TelemetrySnapshot b, double t)
    {
        return new TelemetrySnapshot
        {
            SimTime = a.SimTime + (b.SimTime - a.SimTime) * t,
            Altitude = a.Altitude + (b.Altitude - a.Altitude) * t,
            Speed = a.Speed + (b.Speed - a.Speed) * t,
            Velocity = new Vector(
                a.Velocity.x + (b.Velocity.x - a.Velocity.x) * t,
                a.Velocity.y + (b.Velocity.y - a.Velocity.y) * t,
                a.Velocity.z + (b.Velocity.z - a.Velocity.z) * t),
            Acceleration = a.Acceleration + (b.Acceleration - a.Acceleration) * t,
            FuelPercent = a.FuelPercent + (b.FuelPercent - a.FuelPercent) * t,
            ApoapsisAltitude = a.ApoapsisAltitude + (b.ApoapsisAltitude - a.ApoapsisAltitude) * t,
            PeriapsisAltitude = a.PeriapsisAltitude + (b.PeriapsisAltitude - a.PeriapsisAltitude) * t,
            Mass = a.Mass + (b.Mass - a.Mass) * t
        };
    }
}

/// <summary>
/// A single telemetry snapshot for HUD display.
/// </summary>
public struct TelemetrySnapshot
{
    /// <summary>Simulation time (s).</summary>
    public double SimTime;
    /// <summary>Altitude above surface (m).</summary>
    public double Altitude;
    /// <summary>Speed magnitude (m/s).</summary>
    public double Speed;
    /// <summary>Velocity vector (m/s).</summary>
    public Vector Velocity;
    /// <summary>Acceleration magnitude (m/s²).</summary>
    public double Acceleration;
    /// <summary>Fuel remaining (%).</summary>
    public double FuelPercent;
    /// <summary>Apoapsis altitude (m). Infinity if hyperbolic.</summary>
    public double ApoapsisAltitude;
    /// <summary>Periapsis altitude (m).</summary>
    public double PeriapsisAltitude;
    /// <summary>Vehicle mass (kg).</summary>
    public double Mass;
}
