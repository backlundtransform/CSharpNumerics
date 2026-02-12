using CSharpNumerics.Numerics.Enums;
using Numerics.Models;
using Numerics.Objects;
using System.Collections.Generic;
using System.Linq;

namespace System;

/// <summary>
/// Provides extension methods for numerical differentiation of scalar, multivariate,
/// vector, complex, and discrete data functions.
/// </summary>
public static class DerivativeExtensions
{
    /// <summary>
    /// Default step size used by finite-difference methods.
    /// </summary>
    public const double h = 0.000001;

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

    internal static int Sign(int index) => index % 2 == 0 ? -1 : 1;
}
