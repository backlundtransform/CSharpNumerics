using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Numerics.NumberTheory;

/// <summary>
/// Prime factorisation and the arithmetic functions built on it: divisors, τ(n), σ(n) and μ(n).
/// </summary>
/// <remarks>
/// <para>
/// <see cref="Factor"/> removes small factors by trial division and splits what remains with
/// Pollard's rho (Brent's variant), testing each part with <see cref="PrimalityTesting.IsPrime(long)"/>.
/// Any <see cref="long"/> factors in well under a second. The slowest case is a product of two
/// primes of about 31 bits each.
/// </para>
/// <para>
/// All functions take n ≥ 1 and throw <see cref="ArgumentOutOfRangeException"/> otherwise.
/// </para>
/// </remarks>
public static class IntegerFactorization
{
    /// <summary>Trial division covers primes below this bound before Pollard's rho takes over.</summary>
    private const int TrialDivisionBound = 1000;

    private static readonly int[] SmallPrimes = PrimeSieve.PrimesUpTo(TrialDivisionBound);

    /// <summary>
    /// The prime factorisation of n as (prime, exponent) pairs with the primes in ascending order.
    /// Factor(1) is empty.
    /// </summary>
    /// <example>Factor(360) → [(2, 3), (3, 2), (5, 1)]</example>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> &lt; 1.</exception>
    public static IReadOnlyList<(long Prime, int Exponent)> Factor(long n)
    {
        ValidatePositive(n);

        var factors = new SortedDictionary<long, int>();
        foreach (int p in SmallPrimes)
        {
            if ((long)p * p > n) break;
            int k = 0;
            while (n % p == 0)
            {
                n /= p;
                k++;
            }
            if (k > 0) factors[p] = k;
        }

        // What is left has no prime factor below the bound, so below bound² it is 1 or prime.
        if (n > 1)
        {
            if (n < (long)TrialDivisionBound * TrialDivisionBound)
                Add(factors, n, 1);
            else
                SplitLarge(n, factors);
        }

        return factors.Select(kv => (kv.Key, kv.Value)).ToList();
    }

    /// <summary>
    /// The prime factors of n listed with repetition in ascending order: 12 → [2, 2, 3].
    /// PrimeFactors(1) is empty.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> &lt; 1.</exception>
    public static long[] PrimeFactors(long n)
    {
        var result = new List<long>();
        foreach (var (p, k) in Factor(n))
            for (int i = 0; i < k; i++) result.Add(p);
        return result.ToArray();
    }

    /// <summary>All positive divisors of n in ascending order.</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> &lt; 1.</exception>
    public static long[] Divisors(long n)
    {
        var divisors = new List<long> { 1 };
        foreach (var (p, k) in Factor(n))
        {
            int count = divisors.Count;
            long power = 1;
            for (int e = 1; e <= k; e++)
            {
                power *= p;
                for (int i = 0; i < count; i++) divisors.Add(divisors[i] * power);
            }
        }
        divisors.Sort();
        return divisors.ToArray();
    }

    /// <summary>The number of positive divisors τ(n) (also written d(n) or σ₀(n)).</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> &lt; 1.</exception>
    public static long DivisorCount(long n)
    {
        long count = 1;
        foreach (var (_, k) in Factor(n)) count *= k + 1;
        return count;
    }

    /// <summary>The sum of the positive divisors σ(n) (also written σ₁(n)).</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> &lt; 1.</exception>
    /// <exception cref="OverflowException">σ(n) does not fit in a <see cref="long"/>. σ(n) can exceed n by a factor of about 5 near long.MaxValue.</exception>
    public static long DivisorSum(long n)
    {
        long sum = 1;
        foreach (var (p, k) in Factor(n))
        {
            // 1 + p + … + pᵏ
            long term = 1, power = 1;
            for (int e = 1; e <= k; e++)
            {
                power *= p; // pᵉ divides n, so this cannot overflow
                term = checked(term + power);
            }
            sum = checked(sum * term);
        }
        return sum;
    }

    /// <summary>
    /// The Möbius function μ(n): 0 if n has a squared prime factor, otherwise (−1)ᵏ for k distinct
    /// prime factors. μ(1) = 1.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> &lt; 1.</exception>
    public static int Mobius(long n)
    {
        int sign = 1;
        foreach (var (_, k) in Factor(n))
        {
            if (k > 1) return 0;
            sign = -sign;
        }
        return sign;
    }

    /// <summary>Factors n, which has no prime factor below the trial-division bound, into <paramref name="factors"/>.</summary>
    private static void SplitLarge(long n, SortedDictionary<long, int> factors)
    {
        var pending = new Stack<long>();
        pending.Push(n);
        while (pending.Count > 0)
        {
            long m = pending.Pop();
            if (PrimalityTesting.IsPrime(m))
            {
                Add(factors, m, 1);
                continue;
            }

            long d = FindDivisor(m);
            pending.Push(d);
            pending.Push(m / d);
        }
    }

    /// <summary>A nontrivial divisor of the odd composite m, by Pollard's rho with Brent's cycle detection.</summary>
    private static long FindDivisor(long m)
    {
        for (long c = 1; c < 100; c++)
        {
            long d = Brent(m, c);
            if (d != m) return d;
        }

        // Unreachable in practice: each c fails only with small probability.
        throw new InvalidOperationException($"Pollard's rho failed to split {m}.");
    }

    /// <summary>
    /// One run of Brent's variant of Pollard's rho with f(x) = x² + c mod m, starting at x = 2.
    /// Products of |x − y| are batched and checked with one gcd per batch. Returns a divisor of m,
    /// which is m itself when this c fails.
    /// </summary>
    private static long Brent(long m, long c)
    {
        const int batch = 128;

        long y = 2, x = 2, ys = 2, q = 1, g = 1;
        for (long r = 1; g == 1; r *= 2)
        {
            x = y;
            for (long i = 0; i < r; i++) y = Step(y, c, m);

            for (long k = 0; k < r && g == 1; k += batch)
            {
                ys = y;
                long steps = Math.Min(batch, r - k);
                for (long i = 0; i < steps; i++)
                {
                    y = Step(y, c, m);
                    q = ModularArithmetic.MulReduced(q, Math.Abs(x - y), m);
                }
                g = EuclideanAlgorithm.Gcd(q, m);
            }
        }

        if (g == m)
        {
            // The batch overshot (or q hit 0). Step through it one value at a time.
            do
            {
                ys = Step(ys, c, m);
                g = EuclideanAlgorithm.Gcd(Math.Abs(x - ys), m);
            } while (g == 1);
        }
        return g;
    }

    private static long Step(long y, long c, long m) =>
        ModularArithmetic.AddReduced(ModularArithmetic.MulReduced(y, y, m), c, m);

    private static void Add(SortedDictionary<long, int> factors, long p, int k)
    {
        factors.TryGetValue(p, out int existing);
        factors[p] = existing + k;
    }

    private static void ValidatePositive(long n)
    {
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), n, "n must be positive.");
    }
}
