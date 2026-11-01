using System;
using System.Numerics;

namespace CSharpNumerics.Numerics.NumberTheory;

/// <summary>
/// Modular arithmetic: reduction, overflow-safe multiplication, exponentiation, inverses,
/// Euler's totient, Carmichael's function, multiplicative order, primitive roots,
/// Legendre and Jacobi symbols, and modular square roots.
/// </summary>
/// <remarks>
/// <para>
/// Results are always in [0, m), including for negative inputs: Mod(−3, 7) = 4.
/// A modulus m ≤ 0 throws <see cref="ArgumentOutOfRangeException"/>. m = 1 is allowed and every
/// residue is 0.
/// </para>
/// <para>
/// <see cref="long"/> arithmetic never overflows for any modulus up to <see cref="long.MaxValue"/>.
/// </para>
/// <para>
/// <see cref="Totient"/>, <see cref="CarmichaelLambda"/>, <see cref="MultiplicativeOrder"/> and
/// <see cref="PrimitiveRoot"/> factor their argument with <see cref="IntegerFactorization.Factor"/>.
/// </para>
/// <para>
/// Not constant-time: do not use for cryptographic purposes.
/// </para>
/// </remarks>
public static class ModularArithmetic
{
    /// <summary>
    /// Largest m for which a·b with a, b &lt; m cannot overflow a <see cref="long"/>: ⌊√long.MaxValue⌋.
    /// </summary>
    private const long DirectMultiplyLimit = 3037000499;

    #region Basic operations

    /// <summary>The residue of <paramref name="a"/> modulo <paramref name="m"/>, in [0, m).</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="m"/> ≤ 0.</exception>
    public static long Mod(long a, long m)
    {
        ValidateModulus(m);
        return Reduce(a, m);
    }

    /// <summary>The residue of <paramref name="a"/> modulo <paramref name="m"/>, in [0, m).</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="m"/> ≤ 0.</exception>
    public static BigInteger Mod(BigInteger a, BigInteger m)
    {
        ValidateModulus(m);
        return Reduce(a, m);
    }

    /// <summary>(a + b) mod m, in [0, m).</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="m"/> ≤ 0.</exception>
    public static long AddMod(long a, long b, long m)
    {
        ValidateModulus(m);
        return AddReduced(Reduce(a, m), Reduce(b, m), m);
    }

    /// <summary>(a − b) mod m, in [0, m).</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="m"/> ≤ 0.</exception>
    public static long SubtractMod(long a, long b, long m)
    {
        ValidateModulus(m);
        long d = Reduce(a, m) - Reduce(b, m);
        return d < 0 ? d + m : d;
    }

    /// <summary>(a · b) mod m, in [0, m), without overflow for any m.</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="m"/> ≤ 0.</exception>
    public static long MulMod(long a, long b, long m)
    {
        ValidateModulus(m);
        return MulReduced(Reduce(a, m), Reduce(b, m), m);
    }

    #endregion

    #region Exponentiation and inverses

    /// <summary>
    /// bᵉ mod m by square-and-multiply. A negative exponent means (b⁻¹)^|e|.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="m"/> ≤ 0.</exception>
    /// <exception cref="ArgumentException"><paramref name="exponent"/> is negative and b has no inverse modulo m.</exception>
    public static long Pow(long b, long exponent, long m)
    {
        ValidateModulus(m);
        if (m == 1) return 0;

        long baseValue = Reduce(b, m);
        ulong e;
        if (exponent < 0)
        {
            baseValue = Inverse(baseValue, m);
            // −long.MinValue does not fit in a long, but it does in a ulong.
            e = (ulong)(-(exponent + 1)) + 1;
        }
        else
        {
            e = (ulong)exponent;
        }

        return PowReduced(baseValue, e, m);
    }

    /// <summary>
    /// bᵉ mod m. A negative exponent means (b⁻¹)^|e|.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="m"/> ≤ 0.</exception>
    /// <exception cref="ArgumentException"><paramref name="exponent"/> is negative and b has no inverse modulo m.</exception>
    public static BigInteger Pow(BigInteger b, BigInteger exponent, BigInteger m)
    {
        ValidateModulus(m);
        if (m.IsOne) return BigInteger.Zero;

        BigInteger baseValue = Reduce(b, m);
        if (exponent.Sign < 0)
        {
            baseValue = Inverse(baseValue, m);
            exponent = -exponent;
        }
        return BigInteger.ModPow(baseValue, exponent, m);
    }

    /// <summary>
    /// The multiplicative inverse a⁻¹ mod m, in [0, m).
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="m"/> ≤ 0.</exception>
    /// <exception cref="ArgumentException">gcd(a, m) ≠ 1, so no inverse exists.</exception>
    public static long Inverse(long a, long m)
    {
        if (!TryInverse(a, m, out long inverse))
            throw new ArgumentException($"{a} has no inverse modulo {m} because gcd({a}, {m}) ≠ 1.", nameof(a));
        return inverse;
    }

    /// <summary>
    /// The multiplicative inverse a⁻¹ mod m, in [0, m).
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="m"/> ≤ 0.</exception>
    /// <exception cref="ArgumentException">gcd(a, m) ≠ 1, so no inverse exists.</exception>
    public static BigInteger Inverse(BigInteger a, BigInteger m)
    {
        if (!TryInverse(a, m, out BigInteger inverse))
            throw new ArgumentException($"{a} has no inverse modulo {m} because gcd({a}, {m}) ≠ 1.", nameof(a));
        return inverse;
    }

    /// <summary>
    /// Tries to compute a⁻¹ mod m. Returns false when gcd(a, m) ≠ 1.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="m"/> ≤ 0.</exception>
    public static bool TryInverse(long a, long m, out long inverse)
    {
        ValidateModulus(m);
        var (g, x, _) = EuclideanAlgorithm.ExtendedGcd(Reduce(a, m), m);
        if (g != 1)
        {
            inverse = 0;
            return false;
        }
        inverse = Reduce(x, m);
        return true;
    }

    /// <summary>
    /// Tries to compute a⁻¹ mod m. Returns false when gcd(a, m) ≠ 1.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="m"/> ≤ 0.</exception>
    public static bool TryInverse(BigInteger a, BigInteger m, out BigInteger inverse)
    {
        ValidateModulus(m);
        var (g, x, _) = EuclideanAlgorithm.ExtendedGcd(Reduce(a, m), m);
        if (!g.IsOne)
        {
            inverse = BigInteger.Zero;
            return false;
        }
        inverse = Reduce(x, m);
        return true;
    }

    #endregion

    #region Multiplicative group

    /// <summary>
    /// Euler's totient φ(n): the number of integers in [1, n] coprime to n. φ(1) = 1.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> ≤ 0.</exception>
    public static long Totient(long n)
    {
        if (n <= 0) throw new ArgumentOutOfRangeException(nameof(n), n, "n must be positive.");

        long result = n;
        foreach (var (p, _) in IntegerFactorization.Factor(n))
            result = result / p * (p - 1);
        return result;
    }

    /// <summary>
    /// Carmichael's function λ(n): the exponent of the multiplicative group modulo n, that is the
    /// smallest m &gt; 0 with aᵐ ≡ 1 (mod n) for every a coprime to n. λ(1) = 1.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> ≤ 0.</exception>
    public static long CarmichaelLambda(long n)
    {
        if (n <= 0) throw new ArgumentOutOfRangeException(nameof(n), n, "n must be positive.");

        long result = 1;
        foreach (var (p, k) in IntegerFactorization.Factor(n))
        {
            long pk1 = 1;
            for (int i = 1; i < k; i++) pk1 *= p;

            long lambda;
            if (p == 2 && k >= 3)
                lambda = pk1 / 2;            // λ(2ᵏ) = 2ᵏ⁻² for k ≥ 3
            else
                lambda = pk1 * (p - 1);      // λ(pᵏ) = pᵏ⁻¹(p − 1)

            result = EuclideanAlgorithm.Lcm(result, lambda);
        }
        return result;
    }

    /// <summary>
    /// The multiplicative order of a modulo m: the smallest k &gt; 0 with aᵏ ≡ 1 (mod m).
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="m"/> ≤ 0.</exception>
    /// <exception cref="ArgumentException">gcd(a, m) ≠ 1, so a has no order.</exception>
    public static long MultiplicativeOrder(long a, long m)
    {
        ValidateModulus(m);
        long r = Reduce(a, m);
        if (EuclideanAlgorithm.Gcd(r, m) != 1)
            throw new ArgumentException($"gcd({a}, {m}) ≠ 1, so {a} has no multiplicative order modulo {m}.", nameof(a));
        if (m == 1) return 1;

        // The order divides λ(m). Strip prime factors from λ(m) while the power stays 1.
        long order = CarmichaelLambda(m);
        foreach (var (q, _) in IntegerFactorization.Factor(order))
        {
            while (order % q == 0 && PowReduced(r, (ulong)(order / q), m) == 1)
                order /= q;
        }
        return order;
    }

    /// <summary>
    /// The smallest primitive root modulo m: a generator of the multiplicative group modulo m.
    /// One exists exactly when m is 1, 2, 4, pᵏ or 2pᵏ for an odd prime p.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="m"/> ≤ 0.</exception>
    /// <exception cref="ArgumentException">No primitive root exists modulo m.</exception>
    public static long PrimitiveRoot(long m)
    {
        if (!TryPrimitiveRoot(m, out long root))
            throw new ArgumentException($"No primitive root exists modulo {m}.", nameof(m));
        return root;
    }

    /// <summary>
    /// Tries to find the smallest primitive root modulo m. Returns false when none exists.
    /// The primitive root modulo 1 is 0, the only residue.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="m"/> ≤ 0.</exception>
    public static bool TryPrimitiveRoot(long m, out long root)
    {
        ValidateModulus(m);
        root = 0;
        if (m == 1) return true;

        // The group is cyclic exactly when its exponent equals its order.
        long phi = Totient(m);
        if (CarmichaelLambda(m) != phi) return false;

        var primeFactors = IntegerFactorization.Factor(phi);
        for (long g = 1; g < m; g++)
        {
            if (EuclideanAlgorithm.Gcd(g, m) != 1) continue;

            bool isGenerator = true;
            foreach (var (q, _) in primeFactors)
            {
                if (PowReduced(g, (ulong)(phi / q), m) == 1)
                {
                    isGenerator = false;
                    break;
                }
            }
            if (isGenerator)
            {
                root = g;
                return true;
            }
        }
        return false;
    }

    #endregion

    #region Quadratic residues

    /// <summary>
    /// The Legendre symbol (a/p) for an odd prime p: 0 if p | a, 1 if a is a nonzero quadratic
    /// residue modulo p, and −1 otherwise. Computed with Euler's criterion.
    /// </summary>
    /// <exception cref="ArgumentException"><paramref name="p"/> is not an odd number ≥ 3, or Euler's criterion shows it is not prime.</exception>
    public static int Legendre(long a, long p)
    {
        if (p < 3 || (p & 1) == 0)
            throw new ArgumentException("p must be an odd prime.", nameof(p));

        long r = PowReduced(Reduce(a, p), (ulong)((p - 1) / 2), p);
        if (r == 0) return 0;
        if (r == 1) return 1;
        if (r == p - 1) return -1;
        throw new ArgumentException($"{p} is not prime.", nameof(p));
    }

    /// <summary>
    /// The Jacobi symbol (a/n) for an odd positive n. Equals the Legendre symbol when n is prime.
    /// (a/n) = −1 proves a is not a square modulo n, but (a/n) = 1 does not prove that it is.
    /// </summary>
    /// <exception cref="ArgumentException"><paramref name="n"/> is not odd and positive.</exception>
    public static int Jacobi(long a, long n)
    {
        if (n <= 0 || (n & 1) == 0)
            throw new ArgumentException("n must be odd and positive.", nameof(n));

        a = Reduce(a, n);
        int result = 1;
        while (a != 0)
        {
            while ((a & 1) == 0)
            {
                a >>= 1;
                long r = n & 7;
                if (r == 3 || r == 5) result = -result;
            }
            (a, n) = (n, a);
            if ((a & 3) == 3 && (n & 3) == 3) result = -result;
            a %= n;
        }
        return n == 1 ? result : 0;
    }

    /// <summary>
    /// A square root of a modulo the prime p (Tonelli–Shanks): r with r² ≡ a (mod p).
    /// Returns the smaller of the two roots r and p − r.
    /// </summary>
    /// <exception cref="ArgumentException"><paramref name="p"/> is not prime, or a is not a quadratic residue modulo p.</exception>
    public static long SquareRoot(long a, long p)
    {
        if (!TrySquareRoot(a, p, out long root))
            throw new ArgumentException($"{a} is not a quadratic residue modulo {p}.", nameof(a));
        return root;
    }

    /// <summary>
    /// Tries to compute a square root of a modulo the prime p (Tonelli–Shanks).
    /// Returns false when a is not a quadratic residue. On success the root is the smaller
    /// of r and p − r.
    /// </summary>
    /// <exception cref="ArgumentException"><paramref name="p"/> is less than 2, or the computation shows it is not prime.</exception>
    public static bool TrySquareRoot(long a, long p, out long root)
    {
        if (p < 2) throw new ArgumentException("p must be a prime.", nameof(p));

        root = 0;
        a = Reduce(a, p);
        if (a == 0) return true;
        if (p == 2)
        {
            root = a;
            return true;
        }
        if (Legendre(a, p) != 1) return false;

        long r;
        if ((p & 3) == 3)
        {
            r = PowReduced(a, (ulong)((p + 1) / 4), p);
        }
        else
        {
            // p − 1 = q · 2ˢ with q odd.
            long q = p - 1;
            int s = 0;
            while ((q & 1) == 0)
            {
                q >>= 1;
                s++;
            }

            // Any quadratic non-residue z. Half the residues are non-residues when p is prime.
            long z = 2;
            while (Legendre(z, p) != -1) z++;

            int mExp = s;
            long c = PowReduced(z, (ulong)q, p);
            long t = PowReduced(a, (ulong)q, p);
            r = PowReduced(a, (ulong)((q + 1) / 2), p);
            while (t != 1)
            {
                // Least i in (0, mExp) with t^(2^i) = 1.
                int i = 0;
                long t2 = t;
                while (t2 != 1)
                {
                    // For a prime p, i < mExp always holds. Reaching mExp means p is composite.
                    if (++i >= mExp) throw new ArgumentException($"{p} is not prime.", nameof(p));
                    t2 = MulReduced(t2, t2, p);
                }

                long b = c;
                for (int j = 0; j < mExp - i - 1; j++) b = MulReduced(b, b, p);

                mExp = i;
                c = MulReduced(b, b, p);
                t = MulReduced(t, c, p);
                r = MulReduced(r, b, p);
            }
        }

        if (MulReduced(r, r, p) != a) throw new ArgumentException($"{p} is not prime.", nameof(p));
        root = Math.Min(r, p - r);
        return true;
    }

    #endregion

    #region Helpers

    private static void ValidateModulus(long m)
    {
        if (m <= 0) throw new ArgumentOutOfRangeException(nameof(m), m, "Modulus must be positive.");
    }

    private static void ValidateModulus(BigInteger m)
    {
        if (m.Sign <= 0) throw new ArgumentOutOfRangeException(nameof(m), m, "Modulus must be positive.");
    }

    private static long Reduce(long a, long m)
    {
        long r = a % m;
        return r < 0 ? r + m : r;
    }

    private static BigInteger Reduce(BigInteger a, BigInteger m)
    {
        BigInteger r = BigInteger.Remainder(a, m);
        return r.Sign < 0 ? r + m : r;
    }

    /// <summary>(a + b) mod m for a, b already in [0, m). The sum is formed in a ulong, so it cannot overflow.</summary>
    internal static long AddReduced(long a, long b, long m)
    {
        ulong sum = (ulong)a + (ulong)b;
        return (long)(sum >= (ulong)m ? sum - (ulong)m : sum);
    }

    /// <summary>(a · b) mod m for a, b already in [0, m), without overflow.</summary>
    internal static long MulReduced(long a, long b, long m)
    {
        if (m <= DirectMultiplyLimit) return a * b % m;

#if NET7_0_OR_GREATER
        return (long)((UInt128)(ulong)a * (ulong)b % (ulong)m);
#else
        // Double-and-add. m ≤ 2⁶³ − 1, so x + x and r + x stay below 2⁶⁴.
        ulong x = (ulong)a, y = (ulong)b, mod = (ulong)m, r = 0;
        while (y != 0)
        {
            if ((y & 1) != 0)
            {
                r += x;
                if (r >= mod) r -= mod;
            }
            x += x;
            if (x >= mod) x -= mod;
            y >>= 1;
        }
        return (long)r;
#endif
    }

    /// <summary>bᵉ mod m for b in [0, m) and m &gt; 1.</summary>
    internal static long PowReduced(long b, ulong e, long m)
    {
        long result = 1;
        while (e != 0)
        {
            if ((e & 1) != 0) result = MulReduced(result, b, m);
            e >>= 1;
            if (e != 0) b = MulReduced(b, b, m);
        }
        return result;
    }

    #endregion
}
