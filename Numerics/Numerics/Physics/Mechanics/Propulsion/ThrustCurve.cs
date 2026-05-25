using System;
using System.Collections.Generic;

namespace CSharpNumerics.Physics.Mechanics.Propulsion;

/// <summary>
/// Thrust vs. time profile for rocket engines.
/// Supports solid boosters (grain geometry curves) and liquid engines (throttle schedules).
/// Uses linear interpolation between data points.
/// </summary>
public class ThrustCurve
{
    private readonly double[] _times;
    private readonly double[] _thrustFractions;

    /// <summary>Total burn duration in seconds.</summary>
    public double BurnDuration => _times[_times.Length - 1];

    /// <summary>
    /// Creates a thrust curve from time-thrust fraction pairs.
    /// Thrust fractions are in range [0, 1] representing fraction of maximum thrust.
    /// Times must be monotonically increasing starting from 0.
    /// </summary>
    /// <param name="times">Time points in seconds.</param>
    /// <param name="thrustFractions">Thrust fraction at each time point (0–1).</param>
    public ThrustCurve(double[] times, double[] thrustFractions)
    {
        if (times == null || times.Length < 2)
            throw new ArgumentException("At least two time points required.", nameof(times));
        if (thrustFractions == null || thrustFractions.Length != times.Length)
            throw new ArgumentException("Thrust fractions must match time points.", nameof(thrustFractions));

        _times = (double[])times.Clone();
        _thrustFractions = (double[])thrustFractions.Clone();
    }

    /// <summary>
    /// Evaluates the thrust fraction at a given time via linear interpolation.
    /// Returns 0 if time is past burn duration.
    /// </summary>
    /// <param name="time">Time since ignition in seconds.</param>
    public double Evaluate(double time)
    {
        if (time <= 0) return _thrustFractions[0];
        if (time >= BurnDuration) return 0;

        // Binary search for the interval
        int lo = 0, hi = _times.Length - 1;
        while (hi - lo > 1)
        {
            int mid = (lo + hi) / 2;
            if (_times[mid] <= time) lo = mid;
            else hi = mid;
        }

        double t0 = _times[lo], t1 = _times[hi];
        double f0 = _thrustFractions[lo], f1 = _thrustFractions[hi];
        double alpha = (time - t0) / (t1 - t0);
        return f0 + alpha * (f1 - f0);
    }

    /// <summary>
    /// Creates a constant thrust curve (liquid engine at fixed throttle).
    /// </summary>
    /// <param name="burnDuration">Duration in seconds.</param>
    /// <param name="thrustFraction">Constant thrust level (default 1.0).</param>
    public static ThrustCurve Constant(double burnDuration, double thrustFraction = 1.0)
    {
        return new ThrustCurve(
            new[] { 0.0, burnDuration },
            new[] { thrustFraction, thrustFraction });
    }

    /// <summary>
    /// Creates a typical solid booster profile: ramp-up, sustain, tail-off.
    /// </summary>
    /// <param name="burnDuration">Total burn time in seconds.</param>
    /// <param name="rampFraction">Fraction of burn time for ramp-up (default 0.05).</param>
    /// <param name="tailFraction">Fraction of burn time for tail-off (default 0.1).</param>
    public static ThrustCurve SolidBooster(double burnDuration, double rampFraction = 0.05, double tailFraction = 0.1)
    {
        double tRamp = burnDuration * rampFraction;
        double tTail = burnDuration * tailFraction;
        double tSustain = burnDuration - tRamp - tTail;

        return new ThrustCurve(
            new[] { 0.0, tRamp, tRamp + tSustain, burnDuration },
            new[] { 0.0, 1.0, 1.0, 0.0 });
    }

    /// <summary>
    /// Creates a throttle bucket profile for Max-Q reduction.
    /// Full thrust → reduced thrust → full thrust.
    /// </summary>
    /// <param name="burnDuration">Total burn time in seconds.</param>
    /// <param name="bucketStart">Time when throttle reduction begins (seconds).</param>
    /// <param name="bucketEnd">Time when throttle returns to full (seconds).</param>
    /// <param name="bucketLevel">Throttle fraction during bucket (e.g. 0.7).</param>
    public static ThrustCurve ThrottleBucket(double burnDuration, double bucketStart, double bucketEnd, double bucketLevel)
    {
        return new ThrustCurve(
            new[] { 0.0, bucketStart, bucketStart + 0.1, bucketEnd - 0.1, bucketEnd, burnDuration },
            new[] { 1.0, 1.0, bucketLevel, bucketLevel, 1.0, 1.0 });
    }
}
