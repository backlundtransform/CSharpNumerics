using System;
using System.Linq;
using System.Numerics;
using CSharpNumerics.Numerics;
using CSharpNumerics.Numerics.NumberTheory;

namespace NumericsTests;

[TestClass]
public class NumberTheoryTests
{
    // Primes used across the tests (verified with deterministic Miller–Rabin).
    private const long P998 = 998244353;            // 119·2²³ + 1, ≡ 1 (mod 8)
    private const long PLarge = 4179340454199820289; // 29·2⁵⁷ + 1, above the direct-multiply limit
    private const long Mersenne61 = 2305843009213693951; // 2⁶¹ − 1

    #region EuclideanAlgorithm

    [TestMethod]
    public void Gcd_KnownValuesAndSigns()
    {
        Assert.AreEqual(12, EuclideanAlgorithm.Gcd(84, 36));
        Assert.AreEqual(12, EuclideanAlgorithm.Gcd(-84, 36));
        Assert.AreEqual(12, EuclideanAlgorithm.Gcd(84, -36));
        Assert.AreEqual(7, EuclideanAlgorithm.Gcd(0, -7));
        Assert.AreEqual(0, EuclideanAlgorithm.Gcd(0, 0));
        Assert.AreEqual(1, EuclideanAlgorithm.Gcd(long.MaxValue, long.MaxValue - 1));
    }

    [TestMethod]
    public void Gcd_LongMinValue_Throws()
    {
        Assert.ThrowsException<OverflowException>(() => EuclideanAlgorithm.Gcd(long.MinValue, 3));
        Assert.ThrowsException<OverflowException>(() => EuclideanAlgorithm.ExtendedGcd(3, long.MinValue));
    }

    [TestMethod]
    public void BinaryGcd_MatchesEuclid()
    {
        var rng = new Random(1);
        for (int i = 0; i < 5000; i++)
        {
            long a = rng.NextInt64(-1_000_000_000_000, 1_000_000_000_000) * rng.Next(1, 64);
            long b = rng.NextInt64(-1_000_000_000_000, 1_000_000_000_000) * rng.Next(1, 64);
            Assert.AreEqual(EuclideanAlgorithm.Gcd(a, b), EuclideanAlgorithm.BinaryGcd(a, b), $"{a}, {b}");
        }
        Assert.AreEqual(0, EuclideanAlgorithm.BinaryGcd(0, 0));
        Assert.AreEqual(48, EuclideanAlgorithm.BinaryGcd(0, -48));
    }

    [TestMethod]
    public void Gcd_BigInteger_IsNonNegative()
    {
        Assert.AreEqual(new BigInteger(12), EuclideanAlgorithm.Gcd(new BigInteger(-84), new BigInteger(-36)));
    }

    [TestMethod]
    public void GcdAndLcm_OfSequences()
    {
        Assert.AreEqual(6, EuclideanAlgorithm.Gcd(new long[] { 12, 18, -30 }));
        Assert.AreEqual(0, EuclideanAlgorithm.Gcd(Array.Empty<long>()));
        Assert.AreEqual(60, EuclideanAlgorithm.Lcm(new long[] { 4, 6, 10 }));
        Assert.AreEqual(1, EuclideanAlgorithm.Lcm(Array.Empty<long>()));
        Assert.AreEqual(2520, EuclideanAlgorithm.Lcm(Enumerable.Range(1, 10).Select(i => (long)i)));
    }

    [TestMethod]
    public void Lcm_KnownValuesAndOverflow()
    {
        Assert.AreEqual(36, EuclideanAlgorithm.Lcm(12, -18));
        Assert.AreEqual(0, EuclideanAlgorithm.Lcm(0, 5));
        Assert.ThrowsException<OverflowException>(() => EuclideanAlgorithm.Lcm(Mersenne61, Mersenne61 - 2));
        Assert.AreEqual(BigInteger.Pow(2, 122) - 1, // (2⁶¹ − 1)(2⁶¹ + 1)
            EuclideanAlgorithm.Lcm(new BigInteger(Mersenne61), new BigInteger(Mersenne61 + 2)));
    }

    [TestMethod]
    public void ExtendedGcd_KnownValue()
    {
        Assert.AreEqual((2L, -9L, 47L), EuclideanAlgorithm.ExtendedGcd(240, 46));
        Assert.AreEqual((5L, -1L, 0L), EuclideanAlgorithm.ExtendedGcd(-5, 0));
        Assert.AreEqual((0L, 1L, 0L), EuclideanAlgorithm.ExtendedGcd(0, 0));
    }

    [TestMethod]
    public void ExtendedGcd_SatisfiesBezoutIdentity()
    {
        var rng = new Random(2);
        for (int i = 0; i < 5000; i++)
        {
            long a = rng.NextInt64(long.MinValue + 1, long.MaxValue);
            long b = rng.NextInt64(long.MinValue + 1, long.MaxValue);
            var (g, x, y) = EuclideanAlgorithm.ExtendedGcd(a, b);

            Assert.AreEqual(EuclideanAlgorithm.Gcd(a, b), g);
            Assert.AreEqual(new BigInteger(g), (BigInteger)a * x + (BigInteger)b * y, $"{a}, {b}");
        }
    }

    [TestMethod]
    public void ExtendedGcd_BigInteger_SatisfiesBezoutIdentity()
    {
        var a = BigInteger.Pow(3, 200) - 1;
        var b = -BigInteger.Pow(2, 150) * 3;
        var (g, x, y) = EuclideanAlgorithm.ExtendedGcd(a, b);
        Assert.AreEqual(BigInteger.GreatestCommonDivisor(a, b), g);
        Assert.AreEqual(g, a * x + b * y);
    }

    #endregion

    #region ModularArithmetic — basic operations

    [TestMethod]
    public void Mod_IsAlwaysNonNegative()
    {
        Assert.AreEqual(4, ModularArithmetic.Mod(-3, 7));
        Assert.AreEqual(0, ModularArithmetic.Mod(-14, 7));
        Assert.AreEqual(0, ModularArithmetic.Mod(123, 1));
        Assert.AreEqual(long.MaxValue - 1, ModularArithmetic.Mod(long.MinValue, long.MaxValue));
        Assert.AreEqual(new BigInteger(4), ModularArithmetic.Mod(new BigInteger(-3), new BigInteger(7)));
    }

    [TestMethod]
    public void NonPositiveModulus_Throws()
    {
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => ModularArithmetic.Mod(3, 0));
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => ModularArithmetic.Pow(3, 2, -5));
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => ModularArithmetic.Inverse(3, 0));
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => ModularArithmetic.Mod(BigInteger.One, BigInteger.Zero));
    }

    [TestMethod]
    public void AddSubtractMulMod_MatchBigIntegerNearLongMaxValue()
    {
        var rng = new Random(3);
        long[] moduli = { 2, 97, 3037000499, 3037000500, Mersenne61, PLarge, long.MaxValue - 24, long.MaxValue };
        foreach (long m in moduli)
        {
            for (int i = 0; i < 2000; i++)
            {
                long a = i < 2 ? (i == 0 ? m - 1 : long.MinValue) : rng.NextInt64(long.MinValue, long.MaxValue);
                long b = i < 2 ? m - 1 : rng.NextInt64(long.MinValue, long.MaxValue);
                BigInteger M = m;

                Assert.AreEqual((long)Ref(Ref(a, M) + Ref(b, M), M), ModularArithmetic.AddMod(a, b, m), $"add {a} {b} {m}");
                Assert.AreEqual((long)Ref(Ref(a, M) - Ref(b, M), M), ModularArithmetic.SubtractMod(a, b, m), $"sub {a} {b} {m}");
                Assert.AreEqual((long)Ref(Ref(a, M) * Ref(b, M), M), ModularArithmetic.MulMod(a, b, m), $"mul {a} {b} {m}");
            }
        }
    }

    #endregion

    #region ModularArithmetic — exponentiation and inverses

    [TestMethod]
    public void Pow_KnownValues()
    {
        Assert.AreEqual(445, ModularArithmetic.Pow(4, 13, 497));
        Assert.AreEqual(1, ModularArithmetic.Pow(5, 0, 7));
        Assert.AreEqual(0, ModularArithmetic.Pow(5, 0, 1));
        Assert.AreEqual(0, ModularArithmetic.Pow(0, 5, 7));
        Assert.AreEqual(4, ModularArithmetic.Pow(-3, 3, 31)); // (−27) mod 31
    }

    [TestMethod]
    public void Pow_MatchesBigIntegerModPow()
    {
        var rng = new Random(4);
        for (int i = 0; i < 2000; i++)
        {
            long m = rng.NextInt64(2, long.MaxValue);
            long b = rng.NextInt64(long.MinValue, long.MaxValue);
            long e = rng.NextInt64(0, long.MaxValue);
            Assert.AreEqual((long)BigInteger.ModPow(Ref(b, m), e, m), ModularArithmetic.Pow(b, e, m), $"{b}^{e} mod {m}");
        }
    }

    [TestMethod]
    public void Pow_FermatsLittleTheorem()
    {
        foreach (long p in new[] { 1000000007L, P998, Mersenne61, PLarge })
        {
            foreach (long a in new[] { 2L, 3L, 12345L, p - 1 })
                Assert.AreEqual(1, ModularArithmetic.Pow(a, p - 1, p), $"{a}^(p−1) mod {p}");
        }
    }

    [TestMethod]
    public void Pow_NegativeExponentUsesInverse()
    {
        Assert.AreEqual(ModularArithmetic.Pow(ModularArithmetic.Inverse(3, 11), 4, 11), ModularArithmetic.Pow(3, -4, 11));

        // e = long.MinValue: |e| = 2⁶³ does not fit in a long.
        long inv = ModularArithmetic.Inverse(5, Mersenne61);
        long expected = (long)BigInteger.ModPow(inv, BigInteger.Pow(2, 63), Mersenne61);
        Assert.AreEqual(expected, ModularArithmetic.Pow(5, long.MinValue, Mersenne61));

        Assert.ThrowsException<ArgumentException>(() => ModularArithmetic.Pow(4, -1, 8));
    }

    [TestMethod]
    public void Inverse_KnownValuesAndFailure()
    {
        Assert.AreEqual(4, ModularArithmetic.Inverse(3, 11));
        Assert.AreEqual(7, ModularArithmetic.Inverse(-3, 11)); // −3 ≡ 8, 8·7 = 56 ≡ 1
        Assert.AreEqual(0, ModularArithmetic.Inverse(5, 1));
        Assert.IsFalse(ModularArithmetic.TryInverse(6, 9, out long none));
        Assert.AreEqual(0, none);
        Assert.ThrowsException<ArgumentException>(() => ModularArithmetic.Inverse(6, 9));
    }

    [TestMethod]
    public void Inverse_IsInverseForLargeModuli()
    {
        var rng = new Random(5);
        foreach (long m in new[] { Mersenne61, PLarge, long.MaxValue })
        {
            for (int i = 0; i < 500; i++)
            {
                long a = rng.NextInt64(1, m);
                if (!ModularArithmetic.TryInverse(a, m, out long inv)) continue;
                Assert.AreEqual(1, ModularArithmetic.MulMod(a, inv, m), $"{a}⁻¹ mod {m}");
            }
        }
    }

    [TestMethod]
    public void BigInteger_PowAndInverse()
    {
        var p = BigInteger.Pow(2, 127) - 1; // Mersenne prime
        Assert.AreEqual(BigInteger.One, ModularArithmetic.Pow(new BigInteger(3), p - 1, p));

        var a = BigInteger.Pow(10, 30) + 7;
        var inv = ModularArithmetic.Inverse(a, p);
        Assert.AreEqual(BigInteger.One, a * inv % p);
        Assert.AreEqual(inv, ModularArithmetic.Pow(a, BigInteger.MinusOne, p));
        Assert.AreEqual(ModularArithmetic.Pow(new BigInteger(-3), new BigInteger(3), new BigInteger(31)), new BigInteger(4));
        Assert.IsFalse(ModularArithmetic.TryInverse(new BigInteger(6), new BigInteger(9), out _));
    }

    #endregion

    #region ModularArithmetic — multiplicative group

    [TestMethod]
    public void Totient_MatchesBruteForce()
    {
        for (long n = 1; n <= 300; n++)
        {
            long expected = n == 1 ? 1 : Enumerable.Range(1, (int)n).Count(k => EuclideanAlgorithm.Gcd(k, n) == 1);
            Assert.AreEqual(expected, ModularArithmetic.Totient(n), $"φ({n})");
        }
        Assert.AreEqual(1000000006, ModularArithmetic.Totient(1000000007));
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => ModularArithmetic.Totient(0));
    }

    [TestMethod]
    public void CarmichaelLambda_MatchesBruteForce()
    {
        for (long n = 1; n <= 300; n++)
            Assert.AreEqual(BruteForceExponent(n), ModularArithmetic.CarmichaelLambda(n), $"λ({n})");

        Assert.AreEqual(80, ModularArithmetic.CarmichaelLambda(561));
        Assert.AreEqual(1L << 18, ModularArithmetic.CarmichaelLambda(1L << 20));
    }

    [TestMethod]
    public void MultiplicativeOrder_MatchesBruteForce()
    {
        for (long m = 1; m <= 120; m++)
        {
            for (long a = 0; a < m; a++)
            {
                if (EuclideanAlgorithm.Gcd(a, m) != 1)
                {
                    Assert.ThrowsException<ArgumentException>(() => ModularArithmetic.MultiplicativeOrder(a, m));
                    continue;
                }
                Assert.AreEqual(BruteForceOrder(a, m), ModularArithmetic.MultiplicativeOrder(a, m), $"ord {a} mod {m}");
            }
        }
        Assert.AreEqual(P998 - 1, ModularArithmetic.MultiplicativeOrder(3, P998)); // 3 is a primitive root
    }

    [TestMethod]
    public void PrimitiveRoot_MatchesBruteForce()
    {
        for (long m = 1; m <= 200; m++)
        {
            long phi = ModularArithmetic.Totient(m);
            long expected = -1;
            for (long g = 0; g < m && expected < 0; g++)
            {
                if (EuclideanAlgorithm.Gcd(g, m) == 1 && BruteForceOrder(g, m) == phi)
                    expected = g;
            }

            bool found = ModularArithmetic.TryPrimitiveRoot(m, out long root);
            Assert.AreEqual(expected >= 0, found, $"existence mod {m}");
            if (found) Assert.AreEqual(expected, root, $"root mod {m}");
        }

        Assert.AreEqual(3, ModularArithmetic.PrimitiveRoot(P998));
        Assert.AreEqual(3, ModularArithmetic.PrimitiveRoot(2 * 49)); // 2·7²
        Assert.ThrowsException<ArgumentException>(() => ModularArithmetic.PrimitiveRoot(15));
    }

    #endregion

    #region ModularArithmetic — quadratic residues

    [TestMethod]
    public void Legendre_MatchesBruteForceSquares()
    {
        foreach (long p in new long[] { 3, 5, 7, 11, 13, 101, 997 })
        {
            var squares = Enumerable.Range(1, (int)p - 1).Select(x => (long)x * x % p).ToHashSet();
            for (long a = -p; a < 2 * p; a++)
            {
                long r = ModularArithmetic.Mod(a, p);
                int expected = r == 0 ? 0 : squares.Contains(r) ? 1 : -1;
                Assert.AreEqual(expected, ModularArithmetic.Legendre(a, p), $"({a}/{p})");
            }
        }
    }

    [TestMethod]
    public void Legendre_RejectsEvenAndDetectsComposite()
    {
        Assert.ThrowsException<ArgumentException>(() => ModularArithmetic.Legendre(3, 2));
        Assert.ThrowsException<ArgumentException>(() => ModularArithmetic.Legendre(3, 10));
        Assert.ThrowsException<ArgumentException>(() => ModularArithmetic.Legendre(2, 15)); // 2⁷ ≡ 8 (mod 15)
    }

    [TestMethod]
    public void Jacobi_KnownValuesAndProperties()
    {
        Assert.AreEqual(-1, ModularArithmetic.Jacobi(1001, 9907));
        Assert.AreEqual(1, ModularArithmetic.Jacobi(5, 1));
        Assert.AreEqual(0, ModularArithmetic.Jacobi(6, 9));
        Assert.AreEqual(1, ModularArithmetic.Jacobi(2, 15)); // = 1, yet 2 is not a square mod 15

        // Equals the Legendre symbol for primes, and is multiplicative in the denominator.
        foreach (long p in new long[] { 3, 5, 7, 997 })
            for (long a = 0; a < 50; a++)
                Assert.AreEqual(ModularArithmetic.Legendre(a, p), ModularArithmetic.Jacobi(a, p));

        for (long a = -30; a < 30; a++)
            Assert.AreEqual(ModularArithmetic.Jacobi(a, 21) * ModularArithmetic.Jacobi(a, 55),
                ModularArithmetic.Jacobi(a, 21 * 55), $"a = {a}");

        Assert.ThrowsException<ArgumentException>(() => ModularArithmetic.Jacobi(3, 8));
        Assert.ThrowsException<ArgumentException>(() => ModularArithmetic.Jacobi(3, -7));
    }

    [TestMethod]
    public void SquareRoot_SmallPrimesAgreeWithBruteForce()
    {
        foreach (long p in new long[] { 2, 3, 5, 13, 17, 41, 97, 193, 257 })
        {
            for (long a = 0; a < p; a++)
            {
                long expected = -1;
                for (long x = 0; x <= p / 2 && expected < 0; x++)
                    if (x * x % p == a) expected = x;

                bool found = ModularArithmetic.TrySquareRoot(a, p, out long root);
                Assert.AreEqual(expected >= 0, found, $"√{a} mod {p}");
                if (found) Assert.AreEqual(expected, root, $"√{a} mod {p}");
            }
        }
    }

    [TestMethod]
    public void SquareRoot_LargePrimes()
    {
        var rng = new Random(6);
        foreach (long p in new[] { 1000000007L, P998, Mersenne61, PLarge })
        {
            for (int i = 0; i < 200; i++)
            {
                long x = rng.NextInt64(0, p);
                long a = ModularArithmetic.MulMod(x, x, p);
                long r = ModularArithmetic.SquareRoot(a, p);
                Assert.IsTrue(r == x || r == p - x || (x == 0 && r == 0), $"√{a} mod {p}");
                Assert.IsTrue(r <= p - r || r == 0);
            }
        }

        Assert.IsFalse(ModularArithmetic.TrySquareRoot(3, P998, out _)); // 3 is a primitive root, so a non-residue
        Assert.ThrowsException<ArgumentException>(() => ModularArithmetic.SquareRoot(3, P998));
    }

    #endregion

    #region ChineseRemainderTheorem

    [TestMethod]
    public void Crt_ClassicExample()
    {
        Assert.AreEqual((23L, 105L), ChineseRemainderTheorem.Solve(new long[] { 2, 3, 2 }, new long[] { 3, 5, 7 }));
        Assert.AreEqual((23L, 105L), ChineseRemainderTheorem.Solve(new long[] { -1, -12, 9 }, new long[] { 3, 5, 7 }));
        Assert.AreEqual((0L, 1L), ChineseRemainderTheorem.Solve(Array.Empty<long>(), Array.Empty<long>()));
    }

    [TestMethod]
    public void Crt_Solve_RejectsInvalidInput()
    {
        Assert.ThrowsException<ArgumentException>(() => ChineseRemainderTheorem.Solve(new long[] { 1, 3 }, new long[] { 4, 6 }));
        Assert.ThrowsException<ArgumentException>(() => ChineseRemainderTheorem.Solve(new long[] { 1 }, new long[] { 4, 5 }));
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => ChineseRemainderTheorem.Solve(new long[] { 1, 1 }, new long[] { 4, -5 }));
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => ChineseRemainderTheorem.TrySolve(new long[] { 1 }, new long[] { 0 }, out _, out _));
        Assert.ThrowsException<ArgumentNullException>(() => ChineseRemainderTheorem.Solve(null!, new long[] { 3 }));
        Assert.ThrowsException<OverflowException>(() =>
            ChineseRemainderTheorem.Solve(new long[] { 1, 2 }, new long[] { Mersenne61, 1000000007 }));
    }

    [TestMethod]
    public void Crt_TrySolve_NonCoprimeModuli()
    {
        Assert.IsTrue(ChineseRemainderTheorem.TrySolve(new long[] { 1, 3 }, new long[] { 4, 6 }, out long r, out long m));
        Assert.AreEqual((9L, 12L), (r, m));

        Assert.IsFalse(ChineseRemainderTheorem.TrySolve(new long[] { 1, 2 }, new long[] { 4, 6 }, out _, out _));
    }

    [TestMethod]
    public void Crt_TrySolve_MatchesBruteForce()
    {
        var rng = new Random(7);
        for (int trial = 0; trial < 1000; trial++)
        {
            int count = rng.Next(1, 4);
            long[] moduli = Enumerable.Range(0, count).Select(_ => (long)rng.Next(1, 13)).ToArray();
            long[] remainders = Enumerable.Range(0, count).Select(_ => (long)rng.Next(-50, 50)).ToArray();

            long lcm = EuclideanAlgorithm.Lcm(moduli);
            long expected = -1;
            for (long x = 0; x < lcm && expected < 0; x++)
                if (moduli.Select((mi, i) => ModularArithmetic.Mod(x - remainders[i], mi)).All(d => d == 0))
                    expected = x;

            bool found = ChineseRemainderTheorem.TrySolve(remainders, moduli, out long r, out long m);
            Assert.AreEqual(expected >= 0, found, $"{string.Join(",", remainders)} / {string.Join(",", moduli)}");
            if (found) Assert.AreEqual((expected, lcm), (r, m));
        }
    }

    [TestMethod]
    public void Crt_LargeModuliNearLongMaxValue()
    {
        long[] moduli = { 1000000007, 998244353, 7 }; // product ≈ 7.0·10¹⁸
        long x = 7_777_777_777_777_777;
        long[] remainders = moduli.Select(mi => x % mi).ToArray();
        var (r, m) = ChineseRemainderTheorem.Solve(remainders, moduli);
        Assert.AreEqual(1000000007L * 998244353 * 7, m);
        Assert.AreEqual(x % m, r);
    }

    [TestMethod]
    public void Crt_BigInteger()
    {
        var moduli = new[] { BigInteger.Pow(2, 127) - 1, BigInteger.Pow(2, 61) - 1, new BigInteger(1000000007) };
        var x = BigInteger.Pow(10, 60) + 12345;
        var remainders = moduli.Select(mi => x % mi).ToArray();

        var (r, m) = ChineseRemainderTheorem.Solve(remainders, moduli);
        Assert.AreEqual(moduli[0] * moduli[1] * moduli[2], m);
        Assert.AreEqual(x % m, r);

        Assert.IsTrue(ChineseRemainderTheorem.TrySolve(
            new[] { BigInteger.One, new BigInteger(3) }, new[] { new BigInteger(4), new BigInteger(6) }, out var r2, out var m2));
        Assert.AreEqual((new BigInteger(9), new BigInteger(12)), (r2, m2));
        Assert.IsFalse(ChineseRemainderTheorem.TrySolve(
            new[] { BigInteger.One, new BigInteger(2) }, new[] { new BigInteger(4), new BigInteger(6) }, out _, out _));
    }

    #endregion

    #region ModInt and extensions

    [TestMethod]
    public void ModInt_Arithmetic()
    {
        var a = new ModInt(5, 7);
        var b = new ModInt(-4, 7); // ≡ 3

        Assert.AreEqual(3, b.Value);
        Assert.AreEqual(new ModInt(1, 7), a + b);
        Assert.AreEqual(new ModInt(2, 7), a - b);
        Assert.AreEqual(new ModInt(1, 7), a * b);
        Assert.AreEqual(new ModInt(4, 7), a / b); // 5 · 3⁻¹ = 5 · 5 = 25 ≡ 4
        Assert.AreEqual(new ModInt(2, 7), -a);
        Assert.AreEqual(new ModInt(0, 7), -new ModInt(0, 7));
        Assert.AreEqual(a * a * a, a.Pow(3));
        Assert.AreEqual(a.Inverse(), a.Pow(-1));
        Assert.AreEqual("5 (mod 7)", a.ToString());
    }

    [TestMethod]
    public void ModInt_LargeModulus()
    {
        var a = new ModInt(long.MaxValue - 1, PLarge);
        var product = a * a.Inverse();
        Assert.AreEqual(new ModInt(1, PLarge), product);
    }

    [TestMethod]
    public void ModInt_EqualityAndUnits()
    {
        Assert.AreEqual(new ModInt(10, 7), new ModInt(3, 7));
        Assert.AreNotEqual(new ModInt(3, 7), new ModInt(3, 8));
        Assert.IsTrue(new ModInt(3, 7) == new ModInt(-4, 7));
        Assert.AreEqual(new ModInt(3, 7).GetHashCode(), new ModInt(10, 7).GetHashCode());
        Assert.IsTrue(new ModInt(5, 9).IsUnit);
        Assert.IsFalse(new ModInt(6, 9).IsUnit);
    }

    [TestMethod]
    public void ModInt_InvalidOperationsThrow()
    {
        Assert.ThrowsException<ArgumentException>(() => new ModInt(3, 7) + new ModInt(3, 8));
        Assert.ThrowsException<ArgumentException>(() => new ModInt(1, 9) / new ModInt(3, 9));
        Assert.ThrowsException<ArgumentException>(() => default(ModInt) * default(ModInt));
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => new ModInt(1, 0));
    }

    [TestMethod]
    public void Extensions_DelegateToStaticApi()
    {
        Assert.AreEqual(12, 84L.Gcd(36));
        Assert.AreEqual(36, 12L.Lcm(18));
        Assert.AreEqual(445, 4L.ModPow(13, 497));
        Assert.AreEqual(4, 3L.ModInverse(11));
        Assert.IsTrue(Mersenne61.IsPrime());
    }

    #endregion

    #region PrimalityTesting

    [TestMethod]
    public void IsPrime_MatchesSieveBelow200000()
    {
        var primes = PrimeSieve.PrimesUpTo(200_000).Select(p => (long)p).ToHashSet();
        for (long n = -5; n <= 200_000; n++)
            Assert.AreEqual(primes.Contains(n), PrimalityTesting.IsPrime(n), $"n = {n}");
    }

    [TestMethod]
    public void IsPrime_RejectsCarmichaelNumbersAndStrongPseudoprimes()
    {
        long[] carmichael = { 561, 1105, 1729, 2465, 2821, 6601, 8911, 41041, 825265, 321197185, 5394826801 };
        foreach (long n in carmichael)
            Assert.IsFalse(PrimalityTesting.IsPrime(n), $"Carmichael {n}");

        // Strong pseudoprimes to base 2: fool a single Miller–Rabin round, not the full test.
        foreach (long n in new long[] { 2047, 3277, 4033, 4681, 8321 })
        {
            Assert.IsTrue(PrimalityTesting.IsStrongProbablePrime(n, 2), $"spsp(2) {n}");
            Assert.IsFalse(PrimalityTesting.IsPrime(n), $"{n}");
        }

        // 149491 · 747451 · 34233211 is a strong pseudoprime to every prime base up to 23.
        const long spspTo23 = 3825123056546413051;
        foreach (long a in new long[] { 2, 3, 5, 7, 11, 13, 17, 19, 23 })
            Assert.IsTrue(PrimalityTesting.IsStrongProbablePrime(spspTo23, a));
        Assert.IsFalse(PrimalityTesting.IsPrime(spspTo23));

        Assert.ThrowsException<ArgumentException>(() => PrimalityTesting.IsStrongProbablePrime(10, 2));
    }

    [TestMethod]
    public void IsPrime_LargeKnownValues()
    {
        Assert.IsTrue(PrimalityTesting.IsPrime(Mersenne61));
        Assert.IsTrue(PrimalityTesting.IsPrime(PLarge));
        Assert.IsTrue(PrimalityTesting.IsPrime(P998));
        Assert.IsTrue(PrimalityTesting.IsPrime(PrimalityTesting.LargestLongPrime));
        Assert.IsFalse(PrimalityTesting.IsPrime(long.MaxValue));                    // 7²·73·127·337·92737·649657
        Assert.IsFalse(PrimalityTesting.IsPrime(2147483647L * 2147483629L));        // two 31-bit primes
        Assert.IsFalse(PrimalityTesting.IsPrime(1000000007L * 1000000007L));       // prime square
    }

    [TestMethod]
    public void IsProbablePrime_BigInteger()
    {
        Assert.IsTrue(PrimalityTesting.IsProbablePrime(BigInteger.Pow(2, 127) - 1, seed: 1));
        Assert.IsTrue(PrimalityTesting.IsProbablePrime(BigInteger.Pow(2, 521) - 1, rounds: 10, seed: 2));
        Assert.IsFalse(PrimalityTesting.IsProbablePrime(BigInteger.Pow(2, 128) + 1, seed: 3)); // Fermat F7
        Assert.IsFalse(PrimalityTesting.IsProbablePrime((BigInteger.Pow(2, 61) - 1) * (BigInteger.Pow(2, 89) - 1), seed: 4));
        Assert.IsFalse(PrimalityTesting.IsProbablePrime(BigInteger.Pow(2, 127) + 1));

        // Values that fit in a long use the deterministic test.
        Assert.IsTrue(PrimalityTesting.IsProbablePrime(Mersenne61));
        Assert.IsFalse(PrimalityTesting.IsProbablePrime(new BigInteger(3825123056546413051)));
        Assert.IsFalse(PrimalityTesting.IsProbablePrime(BigInteger.MinusOne));
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => PrimalityTesting.IsProbablePrime(7, rounds: 0));
    }

    [TestMethod]
    public void NextAndPreviousPrime()
    {
        Assert.AreEqual(2, PrimalityTesting.NextPrime(-10));
        Assert.AreEqual(3, PrimalityTesting.NextPrime(2));
        Assert.AreEqual(17, PrimalityTesting.NextPrime(13));
        Assert.AreEqual(1000000000039, PrimalityTesting.NextPrime(1_000_000_000_000));
        Assert.AreEqual(PrimalityTesting.LargestLongPrime, PrimalityTesting.NextPrime(PrimalityTesting.LargestLongPrime - 1));
        Assert.ThrowsException<OverflowException>(() => PrimalityTesting.NextPrime(PrimalityTesting.LargestLongPrime));

        Assert.AreEqual(2, PrimalityTesting.PreviousPrime(3));
        Assert.AreEqual(13, PrimalityTesting.PreviousPrime(17));
        Assert.AreEqual(PrimalityTesting.LargestLongPrime, PrimalityTesting.PreviousPrime(long.MaxValue));
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => PrimalityTesting.PreviousPrime(2));
    }

    #endregion

    #region PrimeSieve

    [TestMethod]
    public void PrimesUpTo_KnownCounts()
    {
        Assert.AreEqual(0, PrimeSieve.PrimesUpTo(1).Length);
        CollectionAssert.AreEqual(new[] { 2 }, PrimeSieve.PrimesUpTo(2));
        CollectionAssert.AreEqual(new[] { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 }, PrimeSieve.PrimesUpTo(30));
        Assert.AreEqual(25, PrimeSieve.PrimesUpTo(100).Length);
        Assert.AreEqual(78498, PrimeSieve.PrimesUpTo(1_000_000).Length);
        Assert.AreEqual(664579, PrimeSieve.PrimesUpTo(10_000_000).Length);
    }

    [TestMethod]
    public void PrimesInRange_SmallRangesMatchSieve()
    {
        var all = PrimeSieve.PrimesUpTo(10_000).Select(p => (long)p).ToArray();
        var rng = new Random(8);
        for (int i = 0; i < 500; i++)
        {
            long low = rng.Next(-20, 9000), high = low + rng.Next(0, 1000);
            var expected = all.Where(p => p >= low && p <= high).ToArray();
            CollectionAssert.AreEqual(expected, PrimeSieve.PrimesInRange(low, high), $"[{low}, {high}]");
        }
        CollectionAssert.AreEqual(Array.Empty<long>(), PrimeSieve.PrimesInRange(-10, 1));
        CollectionAssert.AreEqual(new long[] { 2 }, PrimeSieve.PrimesInRange(2, 2));
    }

    [TestMethod]
    public void PrimesInRange_SegmentedSieveFarFromZero()
    {
        // √high ≈ 10⁸, so the base primes are streamed block by block.
        const long low = 10_000_000_000_000_000, high = low + 200_000;
        var primes = PrimeSieve.PrimesInRange(low, high);
        Assert.AreEqual(5412, primes.Length);
        foreach (long p in primes) Assert.IsTrue(PrimalityTesting.IsPrime(p), $"{p}");
        Assert.AreEqual(primes[0], PrimalityTesting.NextPrime(low - 1));
        Assert.AreEqual(primes[^1], PrimalityTesting.PreviousPrime(high + 1));
    }

    [TestMethod]
    public void PrimesInRange_NearLongMaxValue()
    {
        var primes = PrimeSieve.PrimesInRange(long.MaxValue - 1000, long.MaxValue);
        var expected = Enumerable.Range(0, 1001).Select(i => long.MaxValue - 1000 + i).Where(PrimalityTesting.IsPrime).ToArray();
        CollectionAssert.AreEqual(expected, primes);
        Assert.AreEqual(PrimalityTesting.LargestLongPrime, primes[^1]);
    }

    [TestMethod]
    public void PrimesInRange_InvalidInput()
    {
        Assert.ThrowsException<ArgumentException>(() => PrimeSieve.PrimesInRange(10, 9));
        Assert.ThrowsException<ArgumentException>(() => PrimeSieve.PrimesInRange(0, PrimeSieve.MaxRangeLength + 1));
        Assert.ThrowsException<ArgumentException>(() => PrimeSieve.PrimesInRange(long.MinValue, long.MaxValue));
    }

    #endregion

    #region IntegerFactorization

    [TestMethod]
    public void Factor_KnownValues()
    {
        CollectionAssert.AreEqual(new (long, int)[] { (2, 3), (3, 2), (5, 1) }, IntegerFactorization.Factor(360).ToArray());
        Assert.AreEqual(0, IntegerFactorization.Factor(1).Count);
        CollectionAssert.AreEqual(new (long, int)[] { (2, 62) }, IntegerFactorization.Factor(1L << 62).ToArray());
        CollectionAssert.AreEqual(new long[] { 71, 839, 1471, 6857 }, IntegerFactorization.PrimeFactors(600851475143));
        CollectionAssert.AreEqual(new (long, int)[] { (7, 2), (73, 1), (127, 1), (337, 1), (92737, 1), (649657, 1) },
            IntegerFactorization.Factor(long.MaxValue).ToArray());
        CollectionAssert.AreEqual(new long[] { 2, 3, 715827883, 2147483647 }, IntegerFactorization.PrimeFactors(long.MaxValue - 1));
        CollectionAssert.AreEqual(new long[] { 2147483629, 2147483647 }, IntegerFactorization.PrimeFactors(2147483647L * 2147483629L));
        CollectionAssert.AreEqual(new long[] { 1000000007, 1000000007 }, IntegerFactorization.PrimeFactors(1000000007L * 1000000007L));
        CollectionAssert.AreEqual(new long[] { 149491, 747451, 34233211 }, IntegerFactorization.PrimeFactors(3825123056546413051));
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => IntegerFactorization.Factor(0));
    }

    [TestMethod]
    public void Factor_RandomValuesMultiplyBack()
    {
        var rng = new Random(9);
        for (int i = 0; i < 1000; i++)
        {
            long n = rng.NextInt64(1, long.MaxValue);
            var factors = IntegerFactorization.Factor(n);

            BigInteger product = BigInteger.One;
            long previous = 1;
            foreach (var (p, k) in factors)
            {
                Assert.IsTrue(p > previous, $"ascending for {n}");
                Assert.IsTrue(PrimalityTesting.IsPrime(p), $"{p} | {n}");
                product *= BigInteger.Pow(p, k);
                previous = p;
            }
            Assert.AreEqual(new BigInteger(n), product, $"n = {n}");
        }
    }

    [TestMethod]
    public void DivisorFunctions_MatchBruteForce()
    {
        for (long n = 1; n <= 500; n++)
        {
            var divisors = Enumerable.Range(1, (int)n).Where(d => n % d == 0).Select(d => (long)d).ToArray();
            CollectionAssert.AreEqual(divisors, IntegerFactorization.Divisors(n), $"divisors({n})");
            Assert.AreEqual(divisors.Length, IntegerFactorization.DivisorCount(n), $"τ({n})");
            Assert.AreEqual(divisors.Sum(), IntegerFactorization.DivisorSum(n), $"σ({n})");

            var factors = IntegerFactorization.Factor(n);
            int expectedMu = factors.Any(f => f.Exponent > 1) ? 0 : (factors.Count % 2 == 0 ? 1 : -1);
            Assert.AreEqual(expectedMu, IntegerFactorization.Mobius(n), $"μ({n})");
        }
    }

    [TestMethod]
    public void DivisorFunctions_LargeValues()
    {
        Assert.AreEqual(long.MaxValue, IntegerFactorization.DivisorSum(1L << 62)); // 2⁶³ − 1
        Assert.ThrowsException<OverflowException>(() => IntegerFactorization.DivisorSum(3L << 61));
        Assert.AreEqual(63, IntegerFactorization.DivisorCount(1L << 62));
        Assert.AreEqual(96, IntegerFactorization.DivisorCount(long.MaxValue)); // 3·2⁵
        Assert.AreEqual(0, IntegerFactorization.Mobius(long.MaxValue));      // 7² divides it
        Assert.AreEqual(-1, IntegerFactorization.Mobius(Mersenne61));
        Assert.AreEqual(1, IntegerFactorization.Mobius(1));

        var divisors = IntegerFactorization.Divisors(long.MaxValue);
        Assert.AreEqual(96, divisors.Length);
        Assert.AreEqual(long.MaxValue, divisors[^1]);
        Assert.IsTrue(divisors.All(d => long.MaxValue % d == 0));
    }

    [TestMethod]
    public void MultiplicativeGroup_UsesFastFactorization()
    {
        Assert.AreEqual(Mersenne61 - 1, ModularArithmetic.Totient(Mersenne61));
        Assert.AreEqual(37, ModularArithmetic.PrimitiveRoot(Mersenne61));
        Assert.AreEqual(Mersenne61 - 1, ModularArithmetic.MultiplicativeOrder(37, Mersenne61));
        Assert.AreEqual(2147483646L * 2147483628L, ModularArithmetic.Totient(2147483647L * 2147483629L));
    }

    #endregion

    #region Backward compatibility

    [TestMethod]
    public void NumericExtensions_PrimeHelpersKeepTheirBehaviour()
    {
        // Both namespaces are in scope here: an int receiver still binds to the int overload.
        Assert.IsTrue(CSharpNumerics.Numerics.NumericExtensions.IsPrime(79));
        Assert.IsTrue(79.IsPrime());
        Assert.IsTrue(int.MaxValue.IsPrime());
        Assert.IsFalse(1.IsPrime());
        Assert.IsFalse((-7).IsPrime());

        CollectionAssert.AreEqual(new[] { 2, 3, 13 }, 78.GetPrimeFactors());
        CollectionAssert.AreEqual(new[] { 2, 2, 3 }, 12.GetPrimeFactors());
        CollectionAssert.AreEqual(new[] { int.MaxValue }, int.MaxValue.GetPrimeFactors());
        Assert.AreEqual(0, 1.GetPrimeFactors().Count);
        Assert.AreEqual(0, 0.GetPrimeFactors().Count);
        Assert.AreEqual(0, (-12).GetPrimeFactors().Count);
    }

    #endregion

    #region Helpers

    private static BigInteger Ref(BigInteger a, BigInteger m)
    {
        var r = a % m;
        return r.Sign < 0 ? r + m : r;
    }

    private static long BruteForceOrder(long a, long m)
    {
        if (m == 1) return 1;
        long x = a % m;
        for (long k = 1; ; k++)
        {
            if (x == 1) return k;
            x = x * a % m;
        }
    }

    private static long BruteForceExponent(long n)
    {
        long exponent = 1;
        for (long a = 1; a < Math.Max(n, 2); a++)
            if (EuclideanAlgorithm.Gcd(a, n) == 1)
                exponent = EuclideanAlgorithm.Lcm(exponent, BruteForceOrder(a, n));
        return exponent;
    }

    #endregion
}
