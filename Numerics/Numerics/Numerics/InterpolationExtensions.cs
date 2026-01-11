using CSharpNumerics.Enums;
using Numerics.Models;
using System.Collections.Generic;

namespace System.Linq;

public static class InterpolationExtensions
{

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
            _ => throw new ArgumentOutOfRangeException(nameof(type))
        };
    }

    public static double LinearInterpolationTimeSerie(this IEnumerable<TimeSerie> ts, DateTime timeStamp)
    {

        return ts.LinearInterpolation(p => (p.TimeStamp.Ticks, p.Value), timeStamp.Ticks);
    }


    public static double LogarithmicInterpolationTimeSerie(this IEnumerable<TimeSerie> ts, DateTime timeStamp)
    {

        return ts.LogarithmicInterpolation(p => (p.TimeStamp.Ticks, p.Value), timeStamp.Ticks);
    }


    public static double LinLogInterpolationTimeSerie(this IEnumerable<TimeSerie> ts, DateTime timeStamp)
    {

        return ts.LinLogInterpolation(p => (p.TimeStamp.Ticks, p.Value), timeStamp.Ticks);
    }

    public static double LogLinInterpolationTimeSerie(this IEnumerable<TimeSerie> ts, DateTime timeStamp)
    {

        return ts.LogLinInterpolation(p => (p.TimeStamp.Ticks, p.Value), timeStamp.Ticks);
    }




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

}
