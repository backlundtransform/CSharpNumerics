using System.Collections.Generic;

namespace System;

public static class AlgebraExtensions
{
    /// <summary>
    /// Computes the factorial of an integer.
    /// n! = n * (n-1) * ... * 2 * 1.
    /// For n &lt;= 1, returns 1.
    /// </summary>
    /// <param name="number">The integer input.</param>
    /// <returns>The factorial of <paramref name="number"/> as a <see cref="double"/>.</returns>
    public static double Factorial(this int number)
    {
        if (number <= 1)
        {
            return 1;
        }

        return number * (number - 1).Factorial();
    }

    /// <summary>
    /// Computes the factorial of a double by converting it to an integer.
    /// This method casts <paramref name="number"/> to <see cref="int"/> and calls the integer overload.
    /// </summary>
    /// <param name="number">The numeric input. The fractional part is discarded by the cast to <see cref="int"/>.</param>
    /// <returns>The factorial of <paramref name="number"/> (after casting to <see cref="int"/>) as a <see cref="double"/>.</returns>
    public static double Factorial(this double number)
    {
        return ((int)number).Factorial();
    }

    /// <summary>
    /// Computes the dot product between two lists of doubles.
    /// dot(a,b) = sum(a[i] * b[i]).
    /// </summary>
    /// <param name="firstList">The first vector.</param>
    /// <param name="secondList">The second vector.</param>
    /// <returns>The dot product of <paramref name="firstList"/> and <paramref name="secondList"/>.</returns>
    public static double Dot(this List<double> firstList, List<double> secondList)
    {
        var value = 0.0;

        for (var j = 0; j < firstList.Count; j++)
        {
            value += firstList[j] * secondList[j];
        }

        return value;
    }

    /// <summary>
    /// Finds a root of a function using the Newton-Raphson method.
    /// The update rule is: x_{n+1} = x_n - f(x_n) / f'(x_n).
    /// </summary>
    /// <param name="func">Function for which to find a root.</param>
    /// <param name="xZero">Initial guess.</param>
    /// <returns>An approximation of x such that func(x) is near zero.</returns>
    public static double NewtonRaphson(this Func<double, double> func, double xZero = 1.0)
    {
        var value = xZero;

        for (var j = 0; j < 100; j++)
        {
            var y = func(value);
            var yPrime = func.Derivate(value);

            value -= y / yPrime;
        }

        return value;
    }
}
