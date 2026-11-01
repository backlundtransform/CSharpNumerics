using System;
using System.Collections.Generic;

namespace CSharpNumerics.Numerics.NumberTheory;

/// <summary>
/// Lists primes with the sieve of Eratosthenes, either all primes up to a limit or the primes in
/// a range [low, high] anywhere below <see cref="long.MaxValue"/> (segmented sieve).
/// </summary>
/// <remarks>
/// Both sieves store only odd numbers, one byte each, so they use about half a byte of memory
/// per number in the range.
/// </remarks>
public static class PrimeSieve
{
    /// <summary>Largest allowed high − low for <see cref="PrimesInRange"/> (2³¹, about 1 GB of sieve).</summary>
    public const long MaxRangeLength = 1L << 31;

    /// <summary>⌊√long.MaxValue⌋.</summary>
    private const long SqrtLongMax = 3037000499;

    /// <summary>Block length (in odd numbers) used when streaming the base primes of a segmented sieve.</summary>
    private const int BlockSize = 1 << 16;

    /// <summary>
    /// All primes p ≤ <paramref name="limit"/> in ascending order. Empty when limit &lt; 2.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="limit"/> is int.MaxValue (the sieve would not fit in an array).</exception>
    public static int[] PrimesUpTo(int limit)
    {
        if (limit < 2) return Array.Empty<int>();
        if (limit == int.MaxValue)
            throw new ArgumentOutOfRangeException(nameof(limit), limit, "limit must be less than int.MaxValue.");

        // composite[i] ↔ 2i + 1
        int size = (limit - 1) / 2 + 1;
        var composite = new bool[size];
        composite[0] = true; // 1 is not prime
        for (long i = 1; (2 * i + 1) * (2 * i + 1) <= limit; i++)
        {
            if (composite[i]) continue;
            long p = 2 * i + 1;
            for (long j = p * p / 2; j < size; j += p)
                composite[j] = true;
        }

        var primes = new List<int>(EstimateCount(limit)) { 2 };
        for (int i = 1; i < size; i++)
            if (!composite[i]) primes.Add(2 * i + 1);
        return primes.ToArray();
    }

    /// <summary>
    /// All primes p with <paramref name="low"/> ≤ p ≤ <paramref name="high"/> in ascending order
    /// (segmented sieve). Works anywhere below <see cref="long.MaxValue"/>.
    /// </summary>
    /// <remarks>
    /// The sieve needs every prime up to √high, which it generates block by block instead of
    /// storing them. When the range is very short compared with √high, each candidate is tested
    /// with <see cref="PrimalityTesting.IsPrime(long)"/> instead, which is faster.
    /// </remarks>
    /// <exception cref="ArgumentException"><paramref name="high"/> &lt; <paramref name="low"/>, or the range is longer than <see cref="MaxRangeLength"/>.</exception>
    public static long[] PrimesInRange(long low, long high)
    {
        if (high < low) throw new ArgumentException("high must be ≥ low.", nameof(high));
        if (low < 0) low = 0; // also keeps high − low from overflowing
        if (high - low > MaxRangeLength)
            throw new ArgumentException($"The range may span at most {MaxRangeLength} numbers.", nameof(high));

        var primes = new List<long>();
        if (high < 2) return primes.ToArray();
        if (low <= 2) primes.Add(2);

        long lo = Math.Max(low, 3);
        if ((lo & 1) == 0) lo++;
        if (lo > high) return primes.ToArray();

        // composite[i] ↔ lo + 2i
        long size = (high - lo) / 2 + 1;
        long root = IntegerSqrt(high);

        if (size < root / 1024)
        {
            // Streaming every prime up to √high would dominate; test the few candidates directly.
            for (long i = 0; i < size; i++)
            {
                long n = lo + 2 * i;
                if (PrimalityTesting.IsPrime(n)) primes.Add(n);
            }
            return primes.ToArray();
        }

        var composite = new bool[size];
        foreach (long p in OddPrimesUpTo(root))
            CrossOff(composite, lo, p);

        for (long i = 0; i < size; i++)
            if (!composite[i]) primes.Add(lo + 2 * i);
        return primes.ToArray();
    }

    /// <summary>⌊√n⌋ for n ≥ 0.</summary>
    internal static long IntegerSqrt(long n)
    {
        long r = (long)Math.Sqrt(n);
        while (r * r > n) r--;
        while (r < SqrtLongMax && (r + 1) * (r + 1) <= n) r++;
        return r;
    }

    /// <summary>
    /// The odd primes up to <paramref name="limit"/> ≤ ⌊√long.MaxValue⌋, generated block by block
    /// so that they are never all held in memory.
    /// </summary>
    private static IEnumerable<long> OddPrimesUpTo(long limit)
    {
        int[] basePrimes = PrimesUpTo((int)IntegerSqrt(limit));
        var block = new bool[BlockSize];

        for (long lo = 3; lo <= limit; lo += 2L * BlockSize)
        {
            long hi = Math.Min(limit, lo + 2L * (BlockSize - 1));
            int length = (int)((hi - lo) / 2 + 1);
            Array.Clear(block, 0, length);

            foreach (int p in basePrimes)
            {
                if ((long)p * p > hi) break;
                if (p != 2) CrossOff(block, lo, p, length);
            }

            for (int i = 0; i < length; i++)
                if (!block[i]) yield return lo + 2L * i;
        }
    }

    private static void CrossOff(bool[] composite, long lo, long p) => CrossOff(composite, lo, p, composite.LongLength);

    /// <summary>
    /// Marks the odd multiples m ≥ p² of the odd prime p in the window composite[i] ↔ lo + 2i
    /// (lo odd, i &lt; length). Written so that no intermediate value can exceed the window.
    /// </summary>
    private static void CrossOff(bool[] composite, long lo, long p, long length)
    {
        long hi = lo + 2 * (length - 1);
        long first = p * p; // p ≤ ⌊√long.MaxValue⌋, so this does not overflow
        if (first > hi) return;
        if (first < lo)
        {
            long offset = (p - lo % p) % p;
            if (offset > hi - lo) return;
            first = lo + offset;
            if (((first - lo) & 1) != 0)
            {
                // first is even: the next multiple is odd.
                if (hi - first < p) return;
                first += p;
            }
        }

        for (long i = (first - lo) / 2; i < length; i += p)
            composite[i] = true;
    }

    /// <summary>Upper estimate of π(n) (Rosser–Schoenfeld) used to size the result list.</summary>
    private static int EstimateCount(int n) =>
        n < 17 ? 6 : (int)Math.Min(int.MaxValue, 1.25506 * n / Math.Log(n));
}
