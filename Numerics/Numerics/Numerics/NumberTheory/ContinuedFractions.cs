using System;
using System.Collections.Generic;
using System.Numerics;

namespace CSharpNumerics.Numerics.NumberTheory;

/// <summary>
/// Simple continued fractions x = a₀ + 1/(a₁ + 1/(a₂ + …)), written [a₀; a₁, a₂, …]:
/// expansion, convergents and best rational approximation.
/// </summary>
/// <remarks>
/// <para>
/// a₀ = ⌊x⌋ may be negative or zero. Every later term is positive. A rational number has a finite
/// expansion, and the last term is never 1 unless the expansion is [1].
/// </para>
/// <para>
/// A <see cref="double"/> is expanded exactly: it is first converted to the rational number it
/// represents, mantissa / 2ᵏ. Its terms agree with those of the real number it approximates only
/// as far as the double's precision allows. For Math.PI that is the first 13 terms.
/// </para>
/// </remarks>
public static class ContinuedFractions
{
    /// <summary>
    /// The continued fraction of <paramref name="numerator"/> / <paramref name="denominator"/>.
    /// </summary>
    /// <example>Expand(415, 93) → [4; 2, 6, 7]</example>
    /// <exception cref="DivideByZeroException"><paramref name="denominator"/> is 0.</exception>
    /// <exception cref="OverflowException">The quotient is 2⁶³ (long.MinValue / −1), so a₀ does not fit in a <see cref="long"/>.</exception>
    public static long[] Expand(long numerator, long denominator)
    {
        if (denominator == 0) throw new DivideByZeroException("The denominator must be nonzero.");

        var terms = new List<long>();
        foreach (BigInteger a in Terms(numerator, denominator))
            terms.Add((long)a); // later terms are at most |denominator|, so only a₀ can overflow
        return terms.ToArray();
    }

    /// <summary>
    /// The continued fraction of the exact value of <paramref name="x"/>, up to
    /// <paramref name="maxTerms"/> terms. The expansion also stops before a term that does not fit
    /// in a <see cref="long"/>, which only happens when the remaining tail is negligible.
    /// </summary>
    /// <remarks>
    /// <paramref name="maxTerms"/> is a <see cref="long"/> so that <c>Expand(415, 93)</c> with integer
    /// literals binds to <see cref="Expand(long, long)"/> instead of being ambiguous.
    /// </remarks>
    /// <example>Expand(Math.PI, 5) → [3; 7, 15, 1, 292]</example>
    /// <exception cref="ArgumentException"><paramref name="x"/> is NaN or infinite.</exception>
    /// <exception cref="ArgumentOutOfRangeException">⌊x⌋ does not fit in a <see cref="long"/>, or <paramref name="maxTerms"/> &lt; 1.</exception>
    public static long[] Expand(double x, long maxTerms = 64)
    {
        if (maxTerms < 1) throw new ArgumentOutOfRangeException(nameof(maxTerms), maxTerms, "maxTerms must be at least 1.");
        var (numerator, denominator) = ToExactRational(x);

        var terms = new List<long>();
        foreach (BigInteger a in Terms(numerator, denominator))
        {
            if (terms.Count == 0 && (a < long.MinValue || a > long.MaxValue))
                throw new ArgumentOutOfRangeException(nameof(x), x, "⌊x⌋ must fit in a long.");
            if (terms.Count == maxTerms || a > long.MaxValue) break;
            terms.Add((long)a);
        }
        return terms.ToArray();
    }

    /// <summary>
    /// The convergents hₙ/kₙ = [a₀; a₁, …, aₙ] of a continued fraction, in order. Each fraction is
    /// in lowest terms with a positive denominator. The last one equals the whole expansion.
    /// </summary>
    /// <example>Convergents([3; 7, 15, 1]) → 3/1, 22/7, 333/106, 355/113</example>
    /// <exception cref="ArgumentNullException"><paramref name="terms"/> is null.</exception>
    /// <exception cref="ArgumentException">A term after the first is not positive.</exception>
    /// <exception cref="OverflowException">A numerator or denominator does not fit in a <see cref="long"/>.</exception>
    public static (long Numerator, long Denominator)[] Convergents(IReadOnlyList<long> terms)
    {
        if (terms == null) throw new ArgumentNullException(nameof(terms));

        var result = new (long, long)[terms.Count];
        long h1 = 1, h2 = 0; // h₋₁, h₋₂
        long k1 = 0, k2 = 1; // k₋₁, k₋₂
        for (int i = 0; i < terms.Count; i++)
        {
            long a = terms[i];
            if (i > 0 && a <= 0)
                throw new ArgumentException($"Term {i} is {a}; every term after the first must be positive.", nameof(terms));

            long h = checked(a * h1 + h2);
            long k = checked(a * k1 + k2);
            result[i] = (h, k);
            (h2, h1) = (h1, h);
            (k2, k1) = (k1, k);
        }
        return result;
    }

    /// <summary>
    /// The fraction p/q closest to <paramref name="numerator"/> / <paramref name="denominator"/>
    /// among all fractions with 1 ≤ q ≤ <paramref name="maxDenominator"/>. On a tie, the one with
    /// the smaller denominator. This is a convergent or a semiconvergent of the continued fraction.
    /// </summary>
    /// <example>BestRationalApproximation(314159, 100000, 1000) → 355/113</example>
    /// <exception cref="DivideByZeroException"><paramref name="denominator"/> is 0.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxDenominator"/> &lt; 1.</exception>
    public static (long Numerator, long Denominator) BestRationalApproximation(long numerator, long denominator, long maxDenominator)
    {
        if (denominator == 0) throw new DivideByZeroException("The denominator must be nonzero.");
        return BestApproximation(numerator, denominator, maxDenominator);
    }

    /// <summary>
    /// The fraction p/q closest to <paramref name="x"/> among all fractions with
    /// 1 ≤ q ≤ <paramref name="maxDenominator"/>. On a tie, the one with the smaller denominator.
    /// </summary>
    /// <example>BestRationalApproximation(Math.PI, 1000) → 355/113</example>
    /// <exception cref="ArgumentException"><paramref name="x"/> is NaN or infinite.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxDenominator"/> &lt; 1.</exception>
    /// <exception cref="OverflowException">The numerator of the result does not fit in a <see cref="long"/>.</exception>
    public static (long Numerator, long Denominator) BestRationalApproximation(double x, long maxDenominator)
    {
        var (numerator, denominator) = ToExactRational(x);
        return BestApproximation(numerator, denominator, maxDenominator);
    }

    private static (long, long) BestApproximation(BigInteger numerator, BigInteger denominator, long maxDenominator)
    {
        if (maxDenominator < 1)
            throw new ArgumentOutOfRangeException(nameof(maxDenominator), maxDenominator, "maxDenominator must be at least 1.");
        if (denominator.Sign < 0)
        {
            numerator = -numerator;
            denominator = -denominator;
        }

        BigInteger h1 = 1, h2 = 0, k1 = 0, k2 = 1;
        foreach (BigInteger a in Terms(numerator, denominator))
        {
            BigInteger k = a * k1 + k2;
            if (k > maxDenominator)
            {
                // Largest semiconvergent (m·h₁ + h₂)/(m·k₁ + k₂) with m < a still within the bound.
                BigInteger m = (maxDenominator - k2) / k1;
                BigInteger hs = m * h1 + h2, ks = m * k1 + k2;

                // Keep the semiconvergent only if it is strictly closer than the last convergent.
                // |x − p/q| compares as |numerator·q − p·denominator| / q.
                BigInteger errorSemi = BigInteger.Abs(numerator * ks - hs * denominator) * k1;
                BigInteger errorConv = BigInteger.Abs(numerator * k1 - h1 * denominator) * ks;
                return errorSemi < errorConv ? ((long)hs, (long)ks) : ((long)h1, (long)k1);
            }

            (h2, h1) = (h1, a * h1 + h2);
            (k2, k1) = (k1, k);
        }

        // The expansion ended within the bound: x itself has a small enough denominator.
        return ((long)h1, (long)k1);
    }

    /// <summary>The terms of numerator / denominator via the Euclidean algorithm, with a₀ = ⌊x⌋.</summary>
    private static IEnumerable<BigInteger> Terms(BigInteger numerator, BigInteger denominator)
    {
        if (denominator.Sign < 0)
        {
            numerator = -numerator;
            denominator = -denominator;
        }

        while (!denominator.IsZero)
        {
            BigInteger a = BigInteger.DivRem(numerator, denominator, out BigInteger remainder);
            if (remainder.Sign < 0)
            {
                // Floor division: only possible for a₀, when x is negative.
                a -= 1;
                remainder += denominator;
            }
            yield return a;
            numerator = denominator;
            denominator = remainder;
        }
    }

    /// <summary>The exact value of a finite double as numerator / denominator with denominator a positive power of two.</summary>
    private static (BigInteger Numerator, BigInteger Denominator) ToExactRational(double x)
    {
        if (double.IsNaN(x) || double.IsInfinity(x))
            throw new ArgumentException("x must be finite.", nameof(x));

        long bits = BitConverter.DoubleToInt64Bits(x);
        bool negative = bits < 0;
        int exponent = (int)((bits >> 52) & 0x7FF);
        long mantissa = bits & 0xFFFFFFFFFFFFFL;

        if (exponent == 0) exponent = 1;           // subnormal
        else mantissa |= 1L << 52;                 // implicit leading bit
        exponent -= 1075;                          // value = mantissa · 2^exponent

        BigInteger numerator = negative ? -mantissa : mantissa;
        if (exponent >= 0) return (numerator << exponent, BigInteger.One);
        return (numerator, BigInteger.One << -exponent);
    }
}
