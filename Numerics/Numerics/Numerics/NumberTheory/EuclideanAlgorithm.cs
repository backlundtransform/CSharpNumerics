using System;
using System.Collections.Generic;
using System.Numerics;

namespace CSharpNumerics.Numerics.NumberTheory;

/// <summary>
/// Greatest common divisor, least common multiple and the extended Euclidean algorithm.
/// </summary>
/// <remarks>
/// <para>
/// All results are non-negative, and gcd(0, 0) = 0. For <see cref="long"/> inputs, the absolute
/// value of <see cref="long.MinValue"/> cannot be represented, so it throws
/// <see cref="OverflowException"/>.
/// </para>
/// <para>
/// Not constant-time: do not use for cryptographic purposes.
/// </para>
/// </remarks>
public static class EuclideanAlgorithm
{
    /// <summary>
    /// Greatest common divisor of <paramref name="a"/> and <paramref name="b"/> (Euclid's algorithm).
    /// </summary>
    /// <exception cref="OverflowException">An input is <see cref="long.MinValue"/>.</exception>
    public static long Gcd(long a, long b)
    {
        a = Math.Abs(a);
        b = Math.Abs(b);
        while (b != 0)
        {
            long t = a % b;
            a = b;
            b = t;
        }
        return a;
    }

    /// <summary>
    /// Greatest common divisor of <paramref name="a"/> and <paramref name="b"/>
    /// (non-negative).
    /// </summary>
    public static BigInteger Gcd(BigInteger a, BigInteger b) => BigInteger.GreatestCommonDivisor(a, b);

    /// <summary>
    /// Greatest common divisor of a sequence of integers. An empty sequence gives 0.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="values"/> is null.</exception>
    /// <exception cref="OverflowException">A value is <see cref="long.MinValue"/>.</exception>
    public static long Gcd(IEnumerable<long> values)
    {
        if (values == null) throw new ArgumentNullException(nameof(values));

        long g = 0;
        foreach (long v in values)
        {
            g = Gcd(g, v);
            if (g == 1) break;
        }
        return g;
    }

    /// <summary>
    /// Greatest common divisor using Stein's binary algorithm, which replaces division with
    /// shifts and subtraction. Gives the same result as <see cref="Gcd(long, long)"/>.
    /// </summary>
    /// <exception cref="OverflowException">An input is <see cref="long.MinValue"/>.</exception>
    public static long BinaryGcd(long a, long b)
    {
        ulong u = (ulong)Math.Abs(a);
        ulong v = (ulong)Math.Abs(b);
        if (u == 0) return (long)v;
        if (v == 0) return (long)u;

        int shift = TrailingZeros(u | v);
        u >>= TrailingZeros(u);
        while (v != 0)
        {
            v >>= TrailingZeros(v);
            if (u > v)
            {
                ulong t = u;
                u = v;
                v = t;
            }
            v -= u;
        }
        return (long)(u << shift);
    }

    /// <summary>
    /// Least common multiple of <paramref name="a"/> and <paramref name="b"/> (non-negative).
    /// lcm(a, 0) = 0.
    /// </summary>
    /// <exception cref="OverflowException">An input is <see cref="long.MinValue"/>, or the result does not fit in a <see cref="long"/>.</exception>
    public static long Lcm(long a, long b)
    {
        if (a == 0 || b == 0) return 0;
        a = Math.Abs(a);
        b = Math.Abs(b);
        return checked(a / Gcd(a, b) * b);
    }

    /// <summary>
    /// Least common multiple of <paramref name="a"/> and <paramref name="b"/> (non-negative).
    /// </summary>
    public static BigInteger Lcm(BigInteger a, BigInteger b)
    {
        if (a.IsZero || b.IsZero) return BigInteger.Zero;
        return BigInteger.Abs(a / BigInteger.GreatestCommonDivisor(a, b) * b);
    }

    /// <summary>
    /// Least common multiple of a sequence of integers. An empty sequence gives 1.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="values"/> is null.</exception>
    /// <exception cref="OverflowException">A value is <see cref="long.MinValue"/>, or the result does not fit in a <see cref="long"/>.</exception>
    public static long Lcm(IEnumerable<long> values)
    {
        if (values == null) throw new ArgumentNullException(nameof(values));

        long l = 1;
        foreach (long v in values)
            l = Lcm(l, v);
        return l;
    }

    /// <summary>
    /// Extended Euclidean algorithm: finds g = gcd(a, b) and Bézout coefficients x, y with
    /// a·x + b·y = g.
    /// </summary>
    /// <remarks>
    /// g is non-negative. The coefficients are the minimal ones produced by the algorithm:
    /// |x| ≤ |b / g| and |y| ≤ |a / g|, so they always fit in a <see cref="long"/>.
    /// </remarks>
    /// <exception cref="OverflowException">An input is <see cref="long.MinValue"/>.</exception>
    public static (long Gcd, long X, long Y) ExtendedGcd(long a, long b)
    {
        if (a == long.MinValue || b == long.MinValue)
            throw new OverflowException("The absolute value of long.MinValue cannot be represented.");

        long oldR = a, r = b;
        long oldS = 1, s = 0;
        long oldT = 0, t = 1;
        while (r != 0)
        {
            long q = oldR / r;
            (oldR, r) = (r, oldR - q * r);
            (oldS, s) = (s, oldS - q * s);
            (oldT, t) = (t, oldT - q * t);
        }

        return oldR < 0 ? (-oldR, -oldS, -oldT) : (oldR, oldS, oldT);
    }

    /// <summary>
    /// Extended Euclidean algorithm for <see cref="BigInteger"/>: finds g = gcd(a, b) ≥ 0 and
    /// x, y with a·x + b·y = g.
    /// </summary>
    public static (BigInteger Gcd, BigInteger X, BigInteger Y) ExtendedGcd(BigInteger a, BigInteger b)
    {
        BigInteger oldR = a, r = b;
        BigInteger oldS = BigInteger.One, s = BigInteger.Zero;
        BigInteger oldT = BigInteger.Zero, t = BigInteger.One;
        while (!r.IsZero)
        {
            BigInteger q = BigInteger.Divide(oldR, r);
            (oldR, r) = (r, oldR - q * r);
            (oldS, s) = (s, oldS - q * s);
            (oldT, t) = (t, oldT - q * t);
        }

        return oldR.Sign < 0 ? (-oldR, -oldS, -oldT) : (oldR, oldS, oldT);
    }

    private static int TrailingZeros(ulong x)
    {
        // netstandard2.1 has no BitOperations.TrailingZeroCount. x is never 0 here.
        int n = 0;
        while ((x & 1) == 0)
        {
            x >>= 1;
            n++;
        }
        return n;
    }
}
