using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Data;
using System.Collections.Generic;
using System.Linq;

namespace System;

/// <summary>
/// Provides extension methods for numerical integration:
/// trapezoidal rule, Simpson's rule, Gauss-Legendre quadrature,
/// Romberg integration, adaptive Simpson, Monte Carlo, and discrete data integration.
/// </summary>
public static class IntegrationExtensions
{
    /// <summary>
    /// Default step size used by the trapezoidal integration method.
    /// </summary>
    public const double h = 0.000001;

    #region Trapezoidal Rule

    /// <summary>
    /// Numerically integrates a scalar function using the trapezoidal rule with step size <see cref="h"/>.
    /// </summary>
    /// <param name="func">The function to integrate.</param>
    /// <param name="lowerLimit">Lower integration limit.</param>
    /// <param name="upperLimit">Upper integration limit.</param>
    /// <returns>An approximation of the definite integral.</returns>
    public static double Integrate(this Func<double, double> func, double lowerLimit, double upperLimit)
    {
        var result = 0.0;
        for (var i = lowerLimit; i <= upperLimit; i += h)
        {
            result += h * (func(i) + func(i + h)) / 2;
        }
        return result;
    }

    #endregion

    #region Simpson's Rule

    /// <summary>
    /// Numerically integrates a scalar function using the composite Simpson's 1/3 rule.
    /// Accuracy: O(h⁴) — significantly more accurate than the trapezoidal rule.
    /// </summary>
    /// <param name="func">The function to integrate.</param>
    /// <param name="lowerLimit">Lower integration limit.</param>
    /// <param name="upperLimit">Upper integration limit.</param>
    /// <param name="subintervals">Number of subintervals (must be even; will be rounded up if odd).</param>
    /// <returns>An approximation of the definite integral.</returns>
    public static double IntegrateSimpson(this Func<double, double> func, double lowerLimit, double upperLimit, int subintervals = 1000)
    {
        if (subintervals % 2 != 0) subintervals++;

        double step = (upperLimit - lowerLimit) / subintervals;
        double sum = func(lowerLimit) + func(upperLimit);

        for (int i = 1; i < subintervals; i++)
        {
            double x = lowerLimit + i * step;
            sum += (i % 2 == 0 ? 2 : 4) * func(x);
        }

        return sum * step / 3.0;
    }

    /// <summary>
    /// Numerically integrates a scalar function using the composite Simpson's 3/8 rule.
    /// Uses cubic interpolation on groups of three subintervals. Accuracy: O(h⁴).
    /// </summary>
    /// <param name="func">The function to integrate.</param>
    /// <param name="lowerLimit">Lower integration limit.</param>
    /// <param name="upperLimit">Upper integration limit.</param>
    /// <param name="subintervals">Number of subintervals (must be divisible by 3; will be adjusted if not).</param>
    /// <returns>An approximation of the definite integral.</returns>
    public static double IntegrateSimpson38(this Func<double, double> func, double lowerLimit, double upperLimit, int subintervals = 999)
    {
        while (subintervals % 3 != 0) subintervals++;

        double step = (upperLimit - lowerLimit) / subintervals;
        double sum = func(lowerLimit) + func(upperLimit);

        for (int i = 1; i < subintervals; i++)
        {
            double x = lowerLimit + i * step;
            sum += (i % 3 == 0 ? 2 : 3) * func(x);
        }

        return sum * 3.0 * step / 8.0;
    }

    #endregion

    #region Gauss-Legendre Quadrature

    /// <summary>
    /// Numerically integrates a scalar function using composite Gauss-Legendre quadrature.
    /// Uses 5-point quadrature on each subinterval for very high accuracy on smooth functions.
    /// </summary>
    /// <param name="func">The function to integrate.</param>
    /// <param name="lowerLimit">Lower integration limit.</param>
    /// <param name="upperLimit">Upper integration limit.</param>
    /// <param name="subintervals">Number of subintervals (default 10).</param>
    /// <returns>An approximation of the definite integral.</returns>
    public static double IntegrateGaussLegendre(this Func<double, double> func, double lowerLimit, double upperLimit, int subintervals = 10)
    {
        // 5-point Gauss-Legendre nodes and weights on [-1, 1]
        double[] nodes =
        {
            -0.9061798459386640,
            -0.5384693101056831,
             0.0,
             0.5384693101056831,
             0.9061798459386640
        };
        double[] weights =
        {
            0.2369268850561891,
            0.4786286704993665,
            0.5688888888888889,
            0.4786286704993665,
            0.2369268850561891
        };

        double totalSum = 0.0;
        double subWidth = (upperLimit - lowerLimit) / subintervals;

        for (int s = 0; s < subintervals; s++)
        {
            double a = lowerLimit + s * subWidth;
            double b = a + subWidth;

            double halfWidth = (b - a) / 2.0;
            double midpoint = (a + b) / 2.0;

            double subSum = 0.0;
            for (int i = 0; i < nodes.Length; i++)
            {
                subSum += weights[i] * func(halfWidth * nodes[i] + midpoint);
            }

            totalSum += halfWidth * subSum;
        }

        return totalSum;
    }

    #endregion

    #region Romberg Integration

    /// <summary>
    /// Numerically integrates a scalar function using Romberg integration.
    /// Applies Richardson extrapolation to the trapezoidal rule for rapid convergence.
    /// Accuracy improves exponentially with each level.
    /// </summary>
    /// <param name="func">The function to integrate.</param>
    /// <param name="lowerLimit">Lower integration limit.</param>
    /// <param name="upperLimit">Upper integration limit.</param>
    /// <param name="maxLevel">Maximum number of refinement levels (default 10).</param>
    /// <returns>An approximation of the definite integral.</returns>
    public static double IntegrateRomberg(this Func<double, double> func, double lowerLimit, double upperLimit, int maxLevel = 10)
    {
        double[,] R = new double[maxLevel, maxLevel];

        // R[0,0] = trapezoidal with 1 interval
        R[0, 0] = (upperLimit - lowerLimit) * (func(lowerLimit) + func(upperLimit)) / 2.0;

        for (int n = 1; n < maxLevel; n++)
        {
            // R[n,0] = trapezoidal with 2^n intervals
            int intervals = 1 << n; // 2^n
            double step = (upperLimit - lowerLimit) / intervals;

            double sum = 0.0;
            for (int k = 1; k <= intervals / 2; k++)
            {
                sum += func(lowerLimit + (2 * k - 1) * step);
            }

            R[n, 0] = R[n - 1, 0] / 2.0 + step * sum;

            // Richardson extrapolation
            for (int m = 1; m <= n; m++)
            {
                double factor = Math.Pow(4, m);
                R[n, m] = (factor * R[n, m - 1] - R[n - 1, m - 1]) / (factor - 1.0);
            }
        }

        return R[maxLevel - 1, maxLevel - 1];
    }

    #endregion

    #region Adaptive Simpson

    /// <summary>
    /// Numerically integrates a scalar function using adaptive Simpson's rule.
    /// Recursively subdivides intervals where the error exceeds the specified tolerance.
    /// </summary>
    /// <param name="func">The function to integrate.</param>
    /// <param name="lowerLimit">Lower integration limit.</param>
    /// <param name="upperLimit">Upper integration limit.</param>
    /// <param name="tolerance">Desired absolute accuracy (default 1e-10).</param>
    /// <param name="maxDepth">Maximum recursion depth (default 50).</param>
    /// <returns>An approximation of the definite integral within the specified tolerance.</returns>
    public static double IntegrateAdaptive(this Func<double, double> func, double lowerLimit, double upperLimit, double tolerance = 1e-10, int maxDepth = 50)
    {
        double whole = SimpsonStep(func, lowerLimit, upperLimit);
        return AdaptiveRecursive(func, lowerLimit, upperLimit, tolerance, whole, maxDepth);
    }

    private static double AdaptiveRecursive(Func<double, double> func, double a, double b, double tolerance, double whole, int depth)
    {
        double mid = (a + b) / 2.0;
        double left = SimpsonStep(func, a, mid);
        double right = SimpsonStep(func, mid, b);
        double combined = left + right;

        if (depth <= 0 || Math.Abs(combined - whole) <= 15.0 * tolerance)
        {
            return combined + (combined - whole) / 15.0;
        }

        return AdaptiveRecursive(func, a, mid, tolerance / 2.0, left, depth - 1)
             + AdaptiveRecursive(func, mid, b, tolerance / 2.0, right, depth - 1);
    }

    private static double SimpsonStep(Func<double, double> func, double a, double b)
    {
        double mid = (a + b) / 2.0;
        return (b - a) / 6.0 * (func(a) + 4.0 * func(mid) + func(b));
    }

    #endregion

    #region Multi-Dimensional Integration

    /// <summary>
    /// Numerically integrates a bivariate function over a rectangular region.
    /// Implemented by lifting to a 3D <see cref="Vector"/> function and using the vector integration method.
    /// </summary>
    /// <param name="func">The bivariate function.</param>
    /// <param name="xlimit">X-axis limits.</param>
    /// <param name="ylimit">Y-axis limits.</param>
    /// <returns>An approximation of the double integral over the rectangle.</returns>
    public static double Integrate(
        this Func<(double x, double y), double> func,
        (double lowerLimit, double upperLimit) xlimit,
        (double lowerLimit, double upperLimit) ylimit)
    {
        Func<Vector, double> funcv = (Vector v) => func((v.x, v.y));

        return funcv.Integrate(new Vector(xlimit.lowerLimit, ylimit.lowerLimit, 1), new Vector(xlimit.upperLimit, ylimit.upperLimit, 2));
    }

    /// <summary>
    /// Numerically integrates a function over a 3D rectangular box using Monte Carlo sampling.
    /// </summary>
    /// <param name="func">The function to integrate.</param>
    /// <param name="lowerLimit">Lower bounds (x,y,z).</param>
    /// <param name="upperLimit">Upper bounds (x,y,z).</param>
    /// <returns>An approximation of the triple integral over the box.</returns>
    public static double Integrate(this Func<Vector, double> func, Vector lowerLimit, Vector upperLimit)
    {
        var rnd = new Random();

        var result = 0.0;
        var throws = 999999;

        for (var i = 0; i < throws; i++)
        {
            var x = rnd.NextDouble() * (upperLimit.x - lowerLimit.x) + lowerLimit.x;
            var y = rnd.NextDouble() * (upperLimit.y - lowerLimit.y) + lowerLimit.y;
            var z = rnd.NextDouble() * (upperLimit.z - lowerLimit.z) + lowerLimit.z;

            result += func(new Vector(x, y, z));
        }
        return (upperLimit.x - lowerLimit.x) * (upperLimit.y - lowerLimit.y) * (upperLimit.z - lowerLimit.z) * result / 999999;
    }

    /// <summary>
    /// Numerically integrates a complex function over a rectangular region.
    /// Integrates real and imaginary components separately and returns a complex result.
    /// </summary>
    /// <param name="func">The complex function.</param>
    /// <param name="xlimit">X-axis limits.</param>
    /// <param name="ylimit">Y-axis limits.</param>
    /// <returns>An approximation of the complex integral over the rectangle.</returns>
    public static ComplexNumber Integrate(
        this ComplexFunction func,
        (double lowerLimit, double upperLimit) xlimit,
        (double lowerLimit, double upperLimit) ylimit)
    {
        Func<Vector, double> funcu = (Vector u) => func.u((u.x, u.y));
        Func<Vector, double> funcv = (Vector v) => func.v((v.x, v.y));

        var re = funcu.Integrate(new Vector(xlimit.lowerLimit, ylimit.lowerLimit, 1), new Vector(xlimit.upperLimit, ylimit.upperLimit, 2));

        var img = funcv.Integrate(new Vector(xlimit.lowerLimit, ylimit.lowerLimit, 1), new Vector(xlimit.upperLimit, ylimit.upperLimit, 2));

        return new ComplexNumber(re, img);
    }

    #endregion

    #region Discrete Data Integration

    /// <summary>
    /// Integrates a list of timestamped values by first mapping them to <see cref="TimeSerie"/> and then integrating.
    /// </summary>
    /// <typeparam name="T">The source item type.</typeparam>
    /// <param name="data">The source data.</param>
    /// <param name="func">Mapping from source item to (timestamp, value).</param>
    /// <returns>An approximation of the integral.</returns>
    public static double Integrate<T>(this List<T> data, Func<T, (DateTime timeStamp, double value)> func)
    {
        return data.Select(func).Select(p => new TimeSerie() { TimeStamp = p.timeStamp, Value = p.value }).Integrate();
    }

    /// <summary>
    /// Integrates a list of indexed values by first mapping them to <see cref="Serie"/> and then integrating.
    /// </summary>
    /// <typeparam name="T">The source item type.</typeparam>
    /// <param name="data">The source data.</param>
    /// <param name="func">Mapping from source item to (index, value).</param>
    /// <returns>An approximation of the integral.</returns>
    public static double Integrate<T>(this List<T> data, Func<T, (double index, double value)> func)
    {
        return data.Select(func).Select(p => new Serie() { Index = p.index, Value = p.value }).Integrate();
    }

    /// <summary>
    /// Integrates a time series over the interval [start, end] by converting to seconds-relative indices.
    /// </summary>
    /// <param name="data">The time series points.</param>
    /// <param name="start">Start time.</param>
    /// <param name="end">End time.</param>
    /// <returns>An approximation of the integral over the interval.</returns>
    public static double Integrate(this List<TimeSerie> data, DateTime start, DateTime end)
    {
        return data.Select(p => new Serie() { Index = (p.TimeStamp - start).TotalSeconds, Value = p.Value }).ToList().Integrate(0.0, (end - start).TotalSeconds);
    }

    /// <summary>
    /// Integrates an indexed series over [start, end] using a trapezoidal-like weighting scheme.
    /// </summary>
    /// <param name="data">The series points.</param>
    /// <param name="start">Lower integration limit.</param>
    /// <param name="end">Upper integration limit.</param>
    /// <returns>An approximation of the integral over the interval.</returns>
    public static double Integrate(this List<Serie> data, double start, double end)
    {
        var sum = 0.0;

        if (data.Count == 1)
        {
            return (end - start) * data[0].Value;
        }

        sum += ((data[0].Index - start) + (data[1].Index - data[0].Index) / 2) * data[0].Value;

        for (var i = 1; i < data.Count - 1; i++)
        {
            sum += ((data[i].Index - data[i - 1].Index + (data[i + 1].Index - data[i].Index)) / 2) * data[i].Value;
        }

        sum += (data[data.Count - 1].Index - data[data.Count - 2].Index) / 2 + (end - data[data.Count - 1].Index) * data[data.Count - 1].Value;

        return sum;
    }

    /// <summary>
    /// Integrates a sequence of <see cref="TimeSerie"/> by converting to seconds-relative indices and applying trapezoidal integration.
    /// </summary>
    /// <param name="data">The time series points.</param>
    /// <returns>An approximation of the integral.</returns>
    public static double Integrate(this IEnumerable<TimeSerie> data)
    {
        return data.Select(p => new Serie() { Index = (p.TimeStamp - data.First().TimeStamp).TotalSeconds, Value = p.Value }).ToList().Integrate();
    }

    /// <summary>
    /// Integrates a sequence of <see cref="Serie"/> using the trapezoidal rule.
    /// </summary>
    /// <param name="data">The series points.</param>
    /// <returns>An approximation of the integral.</returns>
    public static double Integrate(this IEnumerable<Serie> data)
    {
        var sum = 0.0;
        var dataList = data.ToList();

        for (var i = 0; i < dataList.Count - 1; i++)
        {
            sum += (dataList[i + 1].Index - dataList[i].Index) * (dataList[i].Value + dataList[i + 1].Value) / 2;
        }

        return sum;
    }

    #endregion
}
