using System.Collections.Generic;


namespace System;

public static class NumericExtensions
{
    public static bool IsPrime(this int n)
    {
        if (n <= 1)
        {
            return false;
        }

        for (int i = 2; i < n; i++)
        {
            if (n % i == 0)
            {
                return false;
            }
        }

        return true;
    }

    public static List<int> GetPrimeFactors(this int n)
    {
        var factors = new List<int>();
        for (int i = 2; i <= n / i; i++)
        {
            while (n % i == 0)
            {
                factors.Add(i);
                n /= i;
            }
        }
        if (n > 1)
        {
            factors.Add(n);
        }
        return factors;
    }

    public static bool IsHappy(this int n)
    {
        var visited = new HashSet<int>();
        while (n != 1 && !visited.Contains(n))
        {
            visited.Add(n);
            int sum = 0;
            while (n != 0)
            {
                int digit = n % 10;
                sum += digit * digit;
                n /= 10;
            }
            n = sum;
        }
        return n == 1;
    }

    public static bool IsPerfectNumber(this int number)
    {
        var sum = 1;

        for (var i = 2; i <= Math.Sqrt(number); i++)
        {
            if (number % i == 0)
            {
                sum += i;
                sum += number / i;
            }
        }

        return (sum == number && number != 1);

    }

    public static int GetDecimalPlaces(this double n)
    {
        n = Math.Abs(n);
        n -= (int)n;
        var decimalPlaces = 0;
        while (n > 0)
        {
            decimalPlaces++;
            n *= 10;
            n -= (int)n;
        }
        return decimalPlaces;
    }


    public static double Abs(this double value)
    {
      
        return Math.Abs(value);
    }


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
