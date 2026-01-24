using CSharpNumerics.Numerics.Enums;
using Numerics.Models;
using Numerics.Objects;
using System.Collections.Generic;
using System.Linq;

namespace System;

public static class CalculusExtensions
{
    /// <summary>
    /// Default step size used by numerical methods in this class.
    /// </summary>
    public const double h = 0.000001;

    /// <summary>
    /// Approximates the left-hand limit of a function at x by evaluating f(x - ε).
    /// </summary>
    /// <param name="function">The function to evaluate.</param>
    /// <param name="x">Point at which to evaluate the left-hand limit.</param>
    /// <returns>An approximation of the left-hand limit at <paramref name="x"/>.</returns>
    public static double LeftLimit(this Func<double, double> function, double x)
    {
        return function(x - double.Epsilon);
    }

    /// <summary>
    /// Approximates the right-hand limit of a function at x by evaluating f(x + ε).
    /// </summary>
    /// <param name="function">The function to evaluate.</param>
    /// <param name="x">Point at which to evaluate the right-hand limit.</param>
    /// <returns>An approximation of the right-hand limit at <paramref name="x"/>.</returns>
    public static double RightLimit(this Func<double, double> function, double x)
    {
        return function(x + double.Epsilon);
    }

    /// <summary>
    /// Approximates the (two-sided) limit of a function at x by comparing left and right evaluations.
    /// If the left and right evaluations are equal, returns that value; otherwise returns NaN.
    /// </summary>
    /// <param name="function">The function to evaluate.</param>
    /// <param name="x">Point at which to evaluate the limit.</param>
    /// <returns>The limit if the left and right evaluations match; otherwise <see cref="double.NaN"/>.</returns>
    public static double Limit(this Func<double, double> function, double x)
    {
        var right = function.LeftLimit(x);

        var left = function.RightLimit(x);

        return (right == left) ? right : double.NaN;
    }

    /// <summary>
    /// Numerically differentiates a scalar function at x.
    /// Uses a finite-difference stencil based on Pascal coefficients to approximate the derivative of the given order.
    /// </summary>
    /// <param name="func">The function to differentiate.</param>
    /// <param name="x">Point at which to differentiate.</param>
    /// <param name="order">Derivative order (default 1).</param>
    /// <returns>An approximation of the derivative of <paramref name="func"/> at <paramref name="x"/>.</returns>
    public static double Derivate(this Func<double, double> func, double x, int order = 1)
    {
        var h0 = Math.Pow(10, order) * h;

        var d = 0.0;

        var pascalMatrix = new Matrix(new double[order + 1, order + 1]).Pascal();

        for (var i = 1; i <= order + 1; i++)
        {
            d += Sign(i) * pascalMatrix.values[order, i - 1] * func(x + (order - i) * h0);
        }

        return d / Math.Pow(h0, order);
    }

    /// <summary>
    /// Numerically computes a partial derivative of a multivariate function.
    /// Differentiates with respect to the variable at <paramref name="index"/>.
    /// </summary>
    /// <param name="func">The multivariate function.</param>
    /// <param name="variables">The point (vector) at which to evaluate the derivative.</param>
    /// <param name="index">Index of the variable to differentiate with respect to.</param>
    /// <param name="order">Derivative order (default 1).</param>
    /// <returns>An approximation of the partial derivative.</returns>
    public static double Derivate(this Func<double[], double> func, double[] variables, int index, int order = 1)
    {
        Func<double, double> funcIndex = (double v) => func(variables.Select((c, i) => { if (i == index) { return v; } return c; }).ToArray());
        return funcIndex.Derivate(variables[index], order);
    }

    /// <summary>
    /// Computes the derivative of a composed expression using the specified operator.
    /// Supports Product rule, Quotient rule, and Chain rule.
    /// </summary>
    /// <param name="funcF">First function.</param>
    /// <param name="funcG">Second function.</param>
    /// <param name="x">Point at which to evaluate the derivative.</param>
    /// <param name="derivateOperator">The derivative operator (product/quotient/chain).</param>
    /// <param name="order">Derivative order (default 1).</param>
    /// <returns>An approximation of the resulting derivative at <paramref name="x"/>.</returns>
    public static double Derivate(
        this Func<double, double> funcF,
        Func<double, double> funcG,
        double x,
        DerivateOperator derivateOperator = DerivateOperator.Chain,
        int order = 1)
    {
        switch (derivateOperator)
        {
            case DerivateOperator.Product:
                return funcF.Derivate(x, order) * funcG(x) + funcF(x) * funcG.Derivate(x, order);
            case DerivateOperator.Quotient:
                return (funcF.Derivate(x, order) * funcG(x) - funcF(x) * funcG.Derivate(x, order)) / Math.Pow(funcG(x), 2);
            case DerivateOperator.Chain:
            default:
                return funcF.Derivate(funcG(x), order) * funcG.Derivate(x, order);
        }
    }

    /// <summary>
    /// Numerically computes a partial derivative of a scalar field f(x,y,z) represented as a function of <see cref="Vector"/>.
    /// </summary>
    /// <param name="func">The scalar field.</param>
    /// <param name="variables">The point at which to evaluate the derivative.</param>
    /// <param name="cartesian">Which coordinate to differentiate with respect to (x, y, or z).</param>
    /// <param name="order">Derivative order (default 1).</param>
    /// <returns>An approximation of the partial derivative.</returns>
    public static double Derivate(this Func<Vector, double> func, Vector variables, Cartesian cartesian, int order = 1)
    {
        switch (cartesian)
        {
            case Cartesian.x:
                Func<double, double> funcX = (double v) => func(new Vector(v, variables.y, variables.z));
                return funcX.Derivate(variables.x, order);

            case Cartesian.y:
                Func<double, double> funcY = (double v) => func(new Vector(variables.x, v, variables.z));
                return funcY.Derivate(variables.y, order);
            case Cartesian.z:
                Func<double, double> funcZ = (double v) => func(new Vector(variables.x, variables.y, v));
                return funcZ.Derivate(variables.z, order);
        }

        return func(variables);
    }

    /// <summary>
    /// Numerically computes a partial derivative of a bivariate function f(x,y).
    /// </summary>
    /// <param name="func">The function to differentiate.</param>
    /// <param name="variables">The point at which to evaluate the derivative.</param>
    /// <param name="cartesian">Which coordinate to differentiate with respect to (x or y).</param>
    /// <param name="order">Derivative order (default 1).</param>
    /// <returns>An approximation of the partial derivative.</returns>
    public static double Derivate(this Func<(double x, double y), double> func, (double x, double y) variables, Cartesian cartesian, int order = 1)
    {
        switch (cartesian)
        {
            case Cartesian.x:
                Func<double, double> funcX = (double v) => func((v, variables.y));
                return funcX.Derivate(variables.x, order);

            case Cartesian.y:
                Func<double, double> funcY = (double v) => func((variables.x, v));
                return funcY.Derivate(variables.y, order);
        }

        return func(variables);
    }

    /// <summary>
    /// Numerically differentiates a complex function.
    /// The derivative is computed by differentiating the real and imaginary components.
    /// </summary>
    /// <param name="func">The complex function.</param>
    /// <param name="variables">The complex point at which to evaluate the derivative.</param>
    /// <param name="order">Derivative order (default 1).</param>
    /// <returns>An approximation of the complex derivative.</returns>
    public static ComplexNumber Derivate(this ComplexFunction func, ComplexNumber variables, int order = 1)
    {
        var du = func.u.Derivate((variables.realPart, variables.imaginaryPart), Cartesian.x, order);

        var dv = func.v.Derivate((variables.realPart, variables.imaginaryPart), Cartesian.x, order);

        return (new ComplexNumber(du, dv));
    }

    /// <summary>
    /// Differentiates a list of timestamped values by first mapping them to <see cref="TimeSerie"/> and then applying numerical differentiation.
    /// </summary>
    /// <typeparam name="T">The source item type.</typeparam>
    /// <param name="data">The source data.</param>
    /// <param name="func">Mapping from source item to (timestamp, value).</param>
    /// <returns>A sequence of derived <see cref="TimeSerie"/> points.</returns>
    public static IEnumerable<TimeSerie> Derivate<T>(this List<T> data, Func<T, (DateTime timeStamp, double value)> func)
    {
        return data.Select(func).Select(p => new TimeSerie() { TimeStamp = p.timeStamp, Value = p.value }).Derivate();
    }

    /// <summary>
    /// Differentiates a list of indexed values by first mapping them to <see cref="Serie"/> and then applying numerical differentiation.
    /// </summary>
    /// <typeparam name="T">The source item type.</typeparam>
    /// <param name="data">The source data.</param>
    /// <param name="func">Mapping from source item to (index, value).</param>
    /// <returns>A sequence of derived <see cref="Serie"/> points.</returns>
    public static IEnumerable<Serie> Derivate<T>(this List<T> data, Func<T, (double index, double value)> func)
    {
        return data.Select(func).Select(p => new Serie() { Index = p.index, Value = p.value }).Derivate();
    }

    /// <summary>
    /// Numerically differentiates a sequence of <see cref="Serie"/> using a forward-difference approximation.
    /// </summary>
    /// <param name="series">The input series.</param>
    /// <returns>A sequence of <see cref="Serie"/> representing the derivative.</returns>
    public static IEnumerable<Serie> Derivate(this IEnumerable<Serie> series)
    {
        var derivativesSeries = new List<Serie>();

        foreach (var serie in series)
        {
            var deltaY = series.LinearInterpolation(p => (p.Index, p.Value), serie.Index + h) - serie.Value;

            derivativesSeries.Add(new Serie() { Value = deltaY / h, Index = serie.Index });
        }
        return derivativesSeries;
    }

    /// <summary>
    /// Numerically differentiates a sequence of <see cref="TimeSerie"/> using a forward-difference approximation.
    /// </summary>
    /// <param name="series">The input time series.</param>
    /// <returns>A sequence of <see cref="TimeSerie"/> representing the derivative.</returns>
    public static IEnumerable<TimeSerie> Derivate(this IEnumerable<TimeSerie> series)
    {
        var derivativesSeries = new List<TimeSerie>();
        foreach (var serie in series)
        {
            var deltaY = series.LinearInterpolation(p => (p.TimeStamp.Ticks, p.Value), serie.TimeStamp.Ticks + h) - serie.Value;

            derivativesSeries.Add(new TimeSerie() { Value = deltaY / h, TimeStamp = serie.TimeStamp });
        }
        return derivativesSeries;
    }

    /// <summary>
    /// Creates a gradient vector field over a diagonal line of points starting at (xmin, ymin, zmin).
    /// Each step evaluates ∇f at (xmin+i, ymin+i, zmin+i).
    /// </summary>
    /// <param name="func">The scalar field.</param>
    /// <param name="xmin">Starting x coordinate.</param>
    /// <param name="ymin">Starting y coordinate.</param>
    /// <param name="zmin">Starting z coordinate.</param>
    /// <param name="stepSize">Increment per step.</param>
    /// <param name="maxSteps">Maximum value of the step parameter.</param>
    /// <returns>A dictionary mapping points to their gradient vectors.</returns>
    public static IDictionary<Vector, Vector> Gradient(this Func<Vector, double> func, double xmin, double ymin, double zmin, double stepSize, double maxSteps)
    {
        var vectorField = new Dictionary<Vector, Vector>();
        for (var i = 0.0; i < maxSteps; i += stepSize)
        {
            vectorField.Add(new Vector((xmin + i, ymin + i, zmin + i)), func.Gradient((xmin + i, ymin + i, zmin + i)));
        }
        return vectorField;
    }

    /// <summary>
    /// Computes the gradient ∇f at a point: (∂f/∂x, ∂f/∂y, ∂f/∂z).
    /// </summary>
    /// <param name="func">The scalar field.</param>
    /// <param name="points">Point at which to compute the gradient.</param>
    /// <returns>The gradient vector at <paramref name="points"/>.</returns>
    public static Vector Gradient(this Func<Vector, double> func, (double, double, double) points)
    {
        var dx = func.Derivate(new Vector(points), Cartesian.x);
        var dy = func.Derivate(new Vector(points), Cartesian.y);
        var dz = func.Derivate(new Vector(points), Cartesian.z);
        return new Vector(dx, dy, dz);
    }

    /// <summary>
    /// Computes the Laplacian ∇²f at a point: ∂²f/∂x² + ∂²f/∂y² + ∂²f/∂z².
    /// </summary>
    /// <param name="func">The scalar field.</param>
    /// <param name="points">Point at which to compute the Laplacian.</param>
    /// <returns>The Laplacian at <paramref name="points"/>.</returns>
    public static double Laplacian(this Func<Vector, double> func, (double, double, double) points)
    {
        var dx2 = func.Derivate(new Vector(points), Cartesian.x, 2);
        var dy2 = func.Derivate(new Vector(points), Cartesian.y, 2);
        var dz2 = func.Derivate(new Vector(points), Cartesian.z, 2);
        return dx2 + dy2 + dz2;
    }

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

    private static int Sign(int index) => index % 2 == 0 ? -1 : 1;
}

