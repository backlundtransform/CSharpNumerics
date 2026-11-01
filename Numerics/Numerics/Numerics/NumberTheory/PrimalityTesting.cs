using System;
using System.Numerics;

namespace CSharpNumerics.Numerics.NumberTheory;

/// <summary>
/// Primality tests (Miller–Rabin) and the next and previous prime.
/// </summary>
/// <remarks>
/// <para>
/// <see cref="IsPrime(long)"/> is deterministic: it uses the first twelve primes as witnesses,
/// which is proven correct for every n &lt; 3.3·10²⁴ and so for every <see cref="long"/>.
/// <see cref="IsProbablePrime"/> handles <see cref="BigInteger"/> with random witnesses. A composite
/// passes each round with probability at most 1/4.
/// </para>
/// <para>
/// Not constant-time: do not use for cryptographic purposes.
/// </para>
/// </remarks>
public static class PrimalityTesting
{
    /// <summary>Largest prime below 2⁶³: 2⁶³ − 25.</summary>
    public const long LargestLongPrime = 9223372036854775783;

    private static readonly long[] Witnesses = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };

    /// <summary>
    /// True when <paramref name="n"/> is prime. Deterministic for every <see cref="long"/>.
    /// Values below 2 are not prime.
    /// </summary>
    public static bool IsPrime(long n)
    {
        if (n < 2) return false;
        foreach (long p in Witnesses)
            if (n % p == 0) return n == p;
        if (n < 37 * 37) return true;

        long d = n - 1;
        int s = 0;
        while ((d & 1) == 0)
        {
            d >>= 1;
            s++;
        }

        foreach (long a in Witnesses)
            if (!IsStrongProbablePrime(n, a, d, s)) return false;
        return true;
    }

    /// <summary>
    /// The Miller–Rabin test of <paramref name="n"/> to the single base <paramref name="witness"/>.
    /// True when n is a strong probable prime to that base. Every odd prime passes. A composite
    /// that passes is a strong pseudoprime to the base, such as 2047 = 23·89 for base 2.
    /// </summary>
    /// <exception cref="ArgumentException"><paramref name="n"/> is not an odd number ≥ 3.</exception>
    public static bool IsStrongProbablePrime(long n, long witness)
    {
        if (n < 3 || (n & 1) == 0) throw new ArgumentException("n must be odd and ≥ 3.", nameof(n));

        long d = n - 1;
        int s = 0;
        while ((d & 1) == 0)
        {
            d >>= 1;
            s++;
        }
        return IsStrongProbablePrime(n, witness, d, s);
    }

    /// <summary>
    /// Probabilistic Miller–Rabin test for <see cref="BigInteger"/>. A prime always returns true.
    /// A composite returns true with probability at most 4^−rounds. Values that fit in a
    /// <see cref="long"/> use the deterministic <see cref="IsPrime(long)"/>.
    /// </summary>
    /// <param name="n">The number to test.</param>
    /// <param name="rounds">Number of random witnesses.</param>
    /// <param name="seed">Seed for the witnesses, for reproducible results. Null uses an unseeded generator.</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="rounds"/> &lt; 1.</exception>
    public static bool IsProbablePrime(BigInteger n, int rounds = 40, int? seed = null)
    {
        if (rounds < 1) throw new ArgumentOutOfRangeException(nameof(rounds), rounds, "rounds must be at least 1.");
        if (n <= long.MaxValue) return IsPrime((long)BigInteger.Max(n, BigInteger.Zero));

        foreach (long p in Witnesses)
            if ((n % p).IsZero) return false;

        BigInteger d = n - 1;
        int s = 0;
        while (d.IsEven)
        {
            d >>= 1;
            s++;
        }

        var random = seed.HasValue ? new Random(seed.Value) : new Random();
        BigInteger nMinus1 = n - 1;
        for (int round = 0; round < rounds; round++)
        {
            BigInteger a = round == 0 ? new BigInteger(2) : RandomWitness(n, random);
            BigInteger x = BigInteger.ModPow(a, d, n);
            if (x.IsOne || x == nMinus1) continue;

            bool composite = true;
            for (int r = 1; r < s; r++)
            {
                x = BigInteger.ModPow(x, 2, n);
                if (x == nMinus1)
                {
                    composite = false;
                    break;
                }
            }
            if (composite) return false;
        }
        return true;
    }

    /// <summary>The smallest prime greater than <paramref name="n"/>.</summary>
    /// <exception cref="OverflowException">n ≥ <see cref="LargestLongPrime"/>, so the next prime does not fit in a <see cref="long"/>.</exception>
    public static long NextPrime(long n)
    {
        if (n < 2) return 2;
        if (n >= LargestLongPrime)
            throw new OverflowException($"There is no prime above {n} that fits in a long.");

        long candidate = (n & 1) == 0 ? n + 1 : n + 2;
        while (!IsPrime(candidate)) candidate += 2;
        return candidate;
    }

    /// <summary>The largest prime less than <paramref name="n"/>.</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> ≤ 2, so there is no smaller prime.</exception>
    public static long PreviousPrime(long n)
    {
        if (n <= 2) throw new ArgumentOutOfRangeException(nameof(n), n, "There is no prime below 2.");
        if (n == 3) return 2;

        long candidate = (n & 1) == 0 ? n - 1 : n - 2;
        while (!IsPrime(candidate)) candidate -= 2;
        return candidate;
    }

    /// <summary>The Miller–Rabin round for odd n ≥ 3 with n − 1 = d·2ˢ, d odd.</summary>
    private static bool IsStrongProbablePrime(long n, long witness, long d, int s)
    {
        long a = ModularArithmetic.Mod(witness, n);
        if (a == 0) return true; // the base is a multiple of n and says nothing

        long x = ModularArithmetic.PowReduced(a, (ulong)d, n);
        if (x == 1 || x == n - 1) return true;
        for (int r = 1; r < s; r++)
        {
            x = ModularArithmetic.MulReduced(x, x, n);
            if (x == n - 1) return true;
        }
        return false;
    }

    /// <summary>A uniformly random witness in [2, n − 2].</summary>
    private static BigInteger RandomWitness(BigInteger n, Random random)
    {
        BigInteger range = n - 3; // number of values in [2, n − 2] minus one
        byte[] bytes = range.ToByteArray(); // little-endian, top byte < 0x80 since range > 0

        // Mask the top byte to the bit length of range, so at least half the draws are accepted.
        byte top = bytes[bytes.Length - 1];
        byte mask = 0;
        while (mask < top) mask = (byte)((mask << 1) | 1);

        BigInteger a;
        do
        {
            random.NextBytes(bytes);
            bytes[bytes.Length - 1] &= mask;
            a = new BigInteger(bytes);
        } while (a > range);
        return a + 2;
    }
}
