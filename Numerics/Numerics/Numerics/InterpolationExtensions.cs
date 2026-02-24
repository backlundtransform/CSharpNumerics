using CSharpNumerics.Numerics.Enums;
using CSharpNumerics.Statistics.Data;
using System;
using System.Collections.Generic;
using System.Linq;

using PolyInterp = CSharpNumerics.Numerics.Interpolation.PolynomialInterpolation;
using SplineInterp = CSharpNumerics.Numerics.Interpolation.CubicSplineInterpolation;
using TrigInterp = CSharpNumerics.Numerics.Interpolation.TrigonometricInterpolation;

namespace CSharpNumerics.Numerics;

public static class InterpolationExtensions
{
    /// <summary>
    /// Interpolates a value at the specified x-coordinate using the requested interpolation type.
    /// </summary>
    /// <typeparam name="T">The source item type.</typeparam>
    /// <param name="ts">The sequence of source items.</param>
    /// <param name="func">Mapping from source item to an (x,y) point.</param>
    /// <param name="index">The x-coordinate at which to interpolate.</param>
    /// <param name="type">The interpolation type.</param>
    /// <returns>The interpolated y-value at <paramref name="index"/>.</returns>
    public static double Interpolate<T>(
        this IEnumerable<T> ts,
        Func<T, (double x, double y)> func,
        double index,
        InterpolationType type)
    {
        return type switch
        {
            InterpolationType.Linear => ts.LinearInterpolation(func, index),
            InterpolationType.Logarithmic => ts.LogarithmicInterpolation(func, index),
            InterpolationType.LinLog => ts.LinLogInterpolation(func, index),
            InterpolationType.LogLin => ts.LogLinInterpolation(func, index),
            InterpolationType.Polynomial => ts.PolynomialInterpolation(func, index),
            InterpolationType.CubicSpline => ts.CubicSplineInterpolation(func, index),
            InterpolationType.Trigonometric => ts.TrigonometricInterpolation(func, index),
            _ => throw new ArgumentOutOfRangeException(nameof(type))
        };
    }

    /// <summary>
    /// Performs linear interpolation on a time series at the specified timestamp.
    /// The interpolation is performed on x = TimeStamp.Ticks and y = Value.
    /// </summary>
    /// <param name="ts">The time series points.</param>
    /// <param name="timeStamp">The timestamp to interpolate at.</param>
    /// <returns>The interpolated value at <paramref name="timeStamp"/>.</returns>
    public static double LinearInterpolationTimeSerie(this IEnumerable<TimeSerie> ts, DateTime timeStamp)
    {
        return ts.LinearInterpolation(p => (p.TimeStamp.Ticks, p.Value), timeStamp.Ticks);
    }

    /// <summary>
    /// Performs logarithmic interpolation (log-log) on a time series at the specified timestamp.
    /// The interpolation is performed on x = TimeStamp.Ticks and y = Value.
    /// </summary>
    /// <param name="ts">The time series points.</param>
    /// <param name="timeStamp">The timestamp to interpolate at.</param>
    /// <returns>The interpolated value at <paramref name="timeStamp"/>.</returns>
    public static double LogarithmicInterpolationTimeSerie(this IEnumerable<TimeSerie> ts, DateTime timeStamp)
    {
        return ts.LogarithmicInterpolation(p => (p.TimeStamp.Ticks, p.Value), timeStamp.Ticks);
    }

    /// <summary>
    /// Performs Lin-Log interpolation on a time series at the specified timestamp.
    /// Linear in x and logarithmic in y.
    /// </summary>
    /// <param name="ts">The time series points.</param>
    /// <param name="timeStamp">The timestamp to interpolate at.</param>
    /// <returns>The interpolated value at <paramref name="timeStamp"/>.</returns>
    public static double LinLogInterpolationTimeSerie(this IEnumerable<TimeSerie> ts, DateTime timeStamp)
    {
        return ts.LinLogInterpolation(p => (p.TimeStamp.Ticks, p.Value), timeStamp.Ticks);
    }

    /// <summary>
    /// Performs Log-Lin interpolation on a time series at the specified timestamp.
    /// Logarithmic in x and linear in y.
    /// </summary>
    /// <param name="ts">The time series points.</param>
    /// <param name="timeStamp">The timestamp to interpolate at.</param>
    /// <returns>The interpolated value at <paramref name="timeStamp"/>.</returns>
    public static double LogLinInterpolationTimeSerie(this IEnumerable<TimeSerie> ts, DateTime timeStamp)
    {
        return ts.LogLinInterpolation(p => (p.TimeStamp.Ticks, p.Value), timeStamp.Ticks);
    }

    /// <summary>
    /// Performs linear interpolation.
    /// Given points (x1,y1) and (x2,y2) around x, returns:
    /// y = y1 + (y2 - y1) * (x - x1) / (x2 - x1).
    /// </summary>
    /// <typeparam name="T">The source item type.</typeparam>
    /// <param name="ts">The sequence of source items.</param>
    /// <param name="func">Mapping from source item to an (x,y) point.</param>
    /// <param name="index">The x-coordinate at which to interpolate.</param>
    /// <returns>The linearly interpolated y-value at <paramref name="index"/>.</returns>
    public static double LinearInterpolation<T>(this IEnumerable<T> ts, Func<T, (double x, double y)> func, double index)
    {
        var points = ts.Select(func).OrderBy(p => p.x).ToList();
        var (x, y) = points.LastOrDefault(p => p.x < index);
        var next = points.FirstOrDefault(p => p.x > index);
        var previousValue = y;
        var nextValue = next.y;

        var previousIndex = x;
        var nextIndex = next.x;

        var currentIndex = index;

        return previousValue + (nextValue - previousValue) * (currentIndex - previousIndex) / (nextIndex - previousIndex);
    }

    /// <summary>
    /// Performs logarithmic interpolation (log-log).
    /// Interpolates linearly in log-space for both x and y.
    /// Requires x &gt; 0 and y &gt; 0 for the bracketing points.
    /// </summary>
    /// <typeparam name="T">The source item type.</typeparam>
    /// <param name="ts">The sequence of source items.</param>
    /// <param name="func">Mapping from source item to an (x,y) point.</param>
    /// <param name="index">The x-coordinate at which to interpolate.</param>
    /// <returns>The interpolated y-value at <paramref name="index"/>.</returns>
    /// <exception cref="ArgumentException">Thrown when x or y values are not positive.</exception>
    public static double LogarithmicInterpolation<T>(
        this IEnumerable<T> ts,
        Func<T, (double x, double y)> func,
        double index)
    {
        var points = ts.Select(func).OrderBy(p => p.x).ToList();
        var (x, y) = points.LastOrDefault(p => p.x < index);
        var next = points.FirstOrDefault(p => p.x > index);

        if (x == 0 || next.x == 0 || y == 0 || next.y == 0)
            throw new ArgumentException("Error is not valid x and y should be > 0");

        double lx1 = Math.Log(x);
        double ly1 = Math.Log(y);
        double lx2 = Math.Log(next.x);
        double ly2 = Math.Log(next.y);
        double lxi = Math.Log(index);

        double t = (lxi - lx1) / (lx2 - lx1);
        double lyi = ly1 + t * (ly2 - ly1);

        return Math.Exp(lyi);
    }

    /// <summary>
    /// Performs Lin-Log interpolation.
    /// Linear in x and logarithmic in y.
    /// Requires y &gt; 0 for the bracketing points.
    /// </summary>
    /// <typeparam name="T">The source item type.</typeparam>
    /// <param name="ts">The sequence of source items.</param>
    /// <param name="func">Mapping from source item to an (x,y) point.</param>
    /// <param name="index">The x-coordinate at which to interpolate.</param>
    /// <returns>The interpolated y-value at <paramref name="index"/>.</returns>
    /// <exception cref="ArgumentException">Thrown when y values are not positive.</exception>
    public static double LinLogInterpolation<T>(
        this IEnumerable<T> ts,
        Func<T, (double x, double y)> func,
        double index)
    {
        var points = ts.Select(func).OrderBy(p => p.x).ToList();
        var (x, y) = points.LastOrDefault(p => p.x < index);
        var next = points.FirstOrDefault(p => p.x > index);

        if (y <= 0 || next.y <= 0)
            throw new ArgumentException("Error is not valid  y should be > 0");

        double x1 = x;
        double x2 = next.x;
        double t = (index - x1) / (x2 - x1);

        double ly1 = Math.Log(y);
        double ly2 = Math.Log(next.y);

        double lyi = ly1 + t * (ly2 - ly1);

        return Math.Exp(lyi);
    }

    /// <summary>
    /// Performs Log-Lin interpolation.
    /// Logarithmic in x and linear in y.
    /// Requires x &gt; 0 for the bracketing points.
    /// </summary>
    /// <typeparam name="T">The source item type.</typeparam>
    /// <param name="ts">The sequence of source items.</param>
    /// <param name="func">Mapping from source item to an (x,y) point.</param>
    /// <param name="index">The x-coordinate at which to interpolate.</param>
    /// <returns>The interpolated y-value at <paramref name="index"/>.</returns>
    /// <exception cref="ArgumentException">Thrown when x values are not positive.</exception>
    public static double LogLinInterpolation<T>(
        this IEnumerable<T> ts,
        Func<T, (double x, double y)> func,
        double index)
    {
        var points = ts.Select(func).OrderBy(p => p.x).ToList();
        var (x, y) = points.LastOrDefault(p => p.x < index);
        var next = points.FirstOrDefault(p => p.x > index);

        if (x <= 0 || next.x <= 0)
            throw new ArgumentException("Error is not valid  x should be > 0");

        double lx1 = Math.Log(x);
        double lx2 = Math.Log(next.x);
        double lxi = Math.Log(index);

        double t = (lxi - lx1) / (lx2 - lx1);

        return y + t * (next.y - y);
    }

    /// <summary>
    /// Linearly interpolates a value based on tick positions.
    /// Returns firstValue when currentTicks == firstTicks, lastValue when currentTicks == lastTicks,
    /// and a proportional blend in between.
    /// </summary>
    /// <param name="currentTicks">The tick position to interpolate at.</param>
    /// <param name="firstTicks">The tick position of the first known value.</param>
    /// <param name="lastTicks">The tick position of the last known value.</param>
    /// <param name="firstValue">The value at <paramref name="firstTicks"/>.</param>
    /// <param name="lastValue">The value at <paramref name="lastTicks"/>.</param>
    /// <returns>The interpolated value at <paramref name="currentTicks"/>.</returns>
    public static double Interpolation(this long currentTicks,
        long firstTicks,
        long lastTicks,
        double firstValue,
        double lastValue)
    {
        if (lastTicks == firstTicks)
        {
            return firstValue;
        }

        return firstValue + (lastValue - firstValue) * (currentTicks - firstTicks) / (lastTicks - firstTicks);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Polynomial interpolation (Lagrange)
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Performs polynomial (Lagrange) interpolation using all data points.
    /// Suitable for small to moderate N. For large N, prefer splines.
    /// </summary>
    /// <typeparam name="T">The source item type.</typeparam>
    /// <param name="ts">The sequence of source items.</param>
    /// <param name="func">Mapping from source item to an (x,y) point.</param>
    /// <param name="index">The x-coordinate at which to interpolate.</param>
    /// <returns>The interpolated y-value.</returns>
    public static double PolynomialInterpolation<T>(
        this IEnumerable<T> ts,
        Func<T, (double x, double y)> func,
        double index)
    {
        var points = ts.Select(func).OrderBy(p => p.x).ToList();
        var x = points.Select(p => p.x).ToArray();
        var y = points.Select(p => p.y).ToArray();
        return PolyInterp.Lagrange(x, y, index);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Cubic spline interpolation (natural)
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Performs natural cubic spline interpolation through all data points.
    /// Produces a smooth (C²-continuous) interpolant.
    /// </summary>
    /// <typeparam name="T">The source item type.</typeparam>
    /// <param name="ts">The sequence of source items.</param>
    /// <param name="func">Mapping from source item to an (x,y) point.</param>
    /// <param name="index">The x-coordinate at which to interpolate.</param>
    /// <returns>The interpolated y-value.</returns>
    public static double CubicSplineInterpolation<T>(
        this IEnumerable<T> ts,
        Func<T, (double x, double y)> func,
        double index)
    {
        var points = ts.Select(func).OrderBy(p => p.x).ToList();
        var x = points.Select(p => p.x).ToArray();
        var y = points.Select(p => p.y).ToArray();
        var spline = new SplineInterp(x, y);
        return spline.Evaluate(index);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Trigonometric interpolation
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Performs trigonometric interpolation. Best for periodic data.
    /// The period is auto-detected from the data range.
    /// </summary>
    /// <typeparam name="T">The source item type.</typeparam>
    /// <param name="ts">The sequence of source items.</param>
    /// <param name="func">Mapping from source item to an (x,y) point.</param>
    /// <param name="index">The x-coordinate at which to interpolate.</param>
    /// <returns>The interpolated y-value.</returns>
    public static double TrigonometricInterpolation<T>(
        this IEnumerable<T> ts,
        Func<T, (double x, double y)> func,
        double index)
    {
        var points = ts.Select(func).OrderBy(p => p.x).ToList();
        var x = points.Select(p => p.x).ToArray();
        var y = points.Select(p => p.y).ToArray();
        var trig = new TrigInterp(x, y);
        return trig.Evaluate(index);
    }
}
