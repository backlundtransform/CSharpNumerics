using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using CSharpNumerics.Numerics.NumberTheory;
using CSharpNumerics.Numerics.NumberTheory.Interfaces;

namespace NumericsTests;

[TestClass]
public class FiniteFieldTests
{
    private const long Mersenne61 = 2305843009213693951; // 2⁶¹ − 1

    #region PrimeField

    [TestMethod]
    public void PrimeField_Basics()
    {
        var gf7 = new PrimeField(7);
        Assert.AreEqual(new ModInt(5, 7), gf7.Inverse(gf7.Element(3)));       // 3·5 = 15 ≡ 1
        Assert.AreEqual(new ModInt(4, 7), gf7.Divide(gf7.Element(1), gf7.Element(2)));
        Assert.AreEqual(new ModInt(4, 7), gf7.Element(-3));
        Assert.AreEqual(new ModInt(3, 7), gf7.PrimitiveElement);
        Assert.AreEqual(3, gf7.MultiplicativeOrder(gf7.Element(2)));          // 2³ = 8 ≡ 1
        Assert.AreEqual(new BigInteger(7), gf7.Order);
        Assert.AreEqual("GF(7)", gf7.ToString());
        Assert.AreEqual(new PrimeField(7), gf7);

        Assert.ThrowsException<ArgumentException>(() => new PrimeField(91));
        Assert.ThrowsException<DivideByZeroException>(() => gf7.Inverse(gf7.Zero));
        Assert.ThrowsException<DivideByZeroException>(() => gf7.Pow(gf7.Zero, -1));
        Assert.ThrowsException<ArgumentException>(() => gf7.Add(gf7.One, new ModInt(1, 11)));
        Assert.AreEqual(gf7.One, gf7.Pow(gf7.Zero, 0));
    }

    [TestMethod]
    public void PrimeField_SatisfiesFieldAxioms()
    {
        AssertFieldAxioms(new PrimeField(2));
        AssertFieldAxioms(new PrimeField(13));
    }

    #endregion

    #region PrimeFieldPolynomial

    [TestMethod]
    public void Polynomial_Arithmetic()
    {
        var a = new PrimeFieldPolynomial(new long[] { 1, 2, 1 }, 3);  // x² + 2x + 1 = (x + 1)²
        var b = new PrimeFieldPolynomial(new long[] { 1, 1 }, 3);     // x + 1
        Assert.AreEqual("x^2 + 2x + 1 (mod 3)", a.ToString());
        Assert.AreEqual(2, a.Degree);
        Assert.AreEqual(b * b, a);
        Assert.AreEqual(b, a / b);
        Assert.IsTrue((a % b).IsZero);
        Assert.AreEqual(Poly(3, 2, 0, 1), a + b);                     // x² + 3x + 2 ≡ x² + 2
        Assert.AreEqual(Poly(3, 0, 1, 1), a - b);
        Assert.AreEqual(Poly(3, 2, 1, 2), 2 * a);
        Assert.AreEqual(Poly(3, 2, 2), a.Derivative());               // 2x + 2
        Assert.AreEqual(1, a.Evaluate(0));
        Assert.AreEqual(0, a.Evaluate(2));                            // root x = −1
        Assert.AreEqual(new ModInt(0, 3), a.Evaluate(new ModInt(-1, 3)));
        Assert.AreEqual(Poly(3, 0, 0, 0, 1).Derivative(), Poly(3));   // 3x² ≡ 0

        var zero = PrimeFieldPolynomial.Zero(3);
        Assert.AreEqual(-1, zero.Degree);
        Assert.AreEqual("0 (mod 3)", zero.ToString());
        Assert.AreEqual(zero, -zero);
        Assert.AreEqual(zero, a + (-a));
        Assert.AreEqual(0, a[7]);

        Assert.AreEqual("x^8 + x^4 + x^3 + x + 1 (mod 2)", PrimeFieldPolynomial.FromInteger(0x11B, 2).ToString());
        Assert.AreEqual(0x11B, PrimeFieldPolynomial.FromInteger(0x11B, 2).ToInteger());
        Assert.AreEqual(Poly(5, 0, 0, 0, 4), PrimeFieldPolynomial.Monomial(-1, 3, 5));

        Assert.ThrowsException<ArgumentException>(() => new PrimeFieldPolynomial(new long[] { 1 }, 4));
        Assert.ThrowsException<ArgumentException>(() => a + Poly(5, 1));
        Assert.ThrowsException<DivideByZeroException>(() => a / zero);
        Assert.ThrowsException<InvalidOperationException>(() => zero.Monic());
    }

    [TestMethod]
    public void Polynomial_DivRemAndExtendedGcd_RandomIdentities()
    {
        var rng = new Random(5);
        foreach (long p in new long[] { 2, 3, 7, 998244353, Mersenne61 })
        {
            for (int trial = 0; trial < 50; trial++)
            {
                var a = RandomPolynomial(p, rng.Next(0, 9), rng);
                var b = RandomPolynomial(p, rng.Next(0, 6), rng);
                if (b.IsZero) continue;

                var (q, r) = PrimeFieldPolynomial.DivRem(a, b);
                Assert.AreEqual(a, q * b + r);
                Assert.IsTrue(r.Degree < b.Degree);

                var (g, s, t) = PrimeFieldPolynomial.ExtendedGcd(a, b);
                Assert.AreEqual(g, a * s + b * t, $"p = {p}");
                Assert.IsTrue(g.IsMonic);
                Assert.IsTrue((a % g).IsZero && (b % g).IsZero);
                Assert.AreEqual(g, PrimeFieldPolynomial.Gcd(a, b));
            }
        }

        // A common factor survives: gcd((x + 1)(x + 2), (x + 1)(x + 3)) = x + 1 over GF(5).
        var common = PrimeFieldPolynomial.Gcd(Poly(5, 1, 1) * Poly(5, 2, 1), Poly(5, 1, 1) * Poly(5, 3, 1));
        Assert.AreEqual(Poly(5, 1, 1), common);
        Assert.IsTrue(PrimeFieldPolynomial.Gcd(PrimeFieldPolynomial.Zero(5), PrimeFieldPolynomial.Zero(5)).IsZero);
    }

    [TestMethod]
    public void Polynomial_PowMod_MatchesRepeatedMultiplication()
    {
        var f = Poly(3, 2, 1, 0, 1);  // x³ + x + 2
        var a = Poly(3, 1, 2, 1, 2);
        var expected = PrimeFieldPolynomial.One(3);
        for (int e = 0; e <= 30; e++)
        {
            Assert.AreEqual(expected % f, a.PowMod(e, f), $"e = {e}");
            expected = expected * a % f;
        }
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => a.PowMod(-1, f));
    }

    [TestMethod]
    public void Polynomial_IrreducibleAndPrimitiveCounts()
    {
        // Gauss: the number of monic irreducibles of degree n over GF(p) is (1/n) Σ μ(d)·p^(n/d) over d | n.
        // A primitive polynomial is the minimal polynomial of a generator: φ(pⁿ − 1)/n of them.
        foreach (var (p, maxDegree) in new[] { (2L, 10), (3L, 5), (5L, 3), (7L, 2) })
        {
            for (int n = 1; n <= maxDegree; n++)
            {
                long irreducible = 0, primitive = 0;
                foreach (var f in MonicPolynomials(p, n))
                {
                    if (f.IsIrreducible()) irreducible++;
                    if (f.IsPrimitive()) primitive++;
                }

                long gauss = IntegerFactorization.Divisors(n)
                    .Sum(d => IntegerFactorization.Mobius(d) * (long)Math.Pow(p, n / d)) / n;
                long order = (long)Math.Pow(p, n) - 1;
                Assert.AreEqual(gauss, irreducible, $"irreducible, p = {p}, n = {n}");
                Assert.AreEqual(ModularArithmetic.Totient(order) / n, primitive, $"primitive, p = {p}, n = {n}");
            }
        }
    }

    [TestMethod]
    public void Polynomial_KnownIrreducibleAndPrimitive()
    {
        var aes = PrimeFieldPolynomial.FromInteger(0x11B, 2);
        Assert.IsTrue(aes.IsIrreducible());
        Assert.IsFalse(aes.IsPrimitive());                                       // x has order 51
        Assert.AreEqual(aes, PrimeFieldPolynomial.FindIrreducible(2, 8));
        Assert.AreEqual(0x11D, PrimeFieldPolynomial.FindPrimitive(2, 8).ToInteger());

        Assert.IsTrue(PrimeFieldPolynomial.FromInteger(0b11111, 2).IsIrreducible()); // x⁴ + x³ + x² + x + 1
        Assert.IsFalse(PrimeFieldPolynomial.FromInteger(0b11111, 2).IsPrimitive());  // x⁵ = 1
        Assert.IsFalse(PrimeFieldPolynomial.FromInteger(0b101, 2).IsIrreducible());   // x² + 1 = (x + 1)²
        Assert.IsFalse((Poly(2, 1, 1, 1) * Poly(2, 1, 1, 0, 0, 1)).IsIrreducible()); // (x² + x + 1)(x⁴ + x + 1): no roots, still reducible
        Assert.IsTrue(Poly(2, 1, 1, 0, 0, 0, 0, 1).IsPrimitive());                  // x⁶ + x + 1

        // Non-monic polynomials are judged up to the unit factor.
        Assert.IsTrue((3 * Poly(5, 2, 0, 1)).IsIrreducible());                     // 3(x² + 2): x² = 3 has no root mod 5
        Assert.IsFalse(PrimeFieldPolynomial.One(5).IsIrreducible());
        Assert.IsFalse(PrimeFieldPolynomial.Zero(5).IsIrreducible());

        // Degree 1: x − a is primitive exactly when a is a primitive root.
        Assert.AreEqual(Poly(7, 2, 1), PrimeFieldPolynomial.FindPrimitive(7, 1)); // x + 2 = x − 5, and 5 is a primitive root
        Assert.AreEqual(Poly(7, 0, 1), PrimeFieldPolynomial.FindIrreducible(7, 1));

        Assert.ThrowsException<ArgumentOutOfRangeException>(() => PrimeFieldPolynomial.FindIrreducible(2, 0));
        Assert.ThrowsException<OverflowException>(() => PrimeFieldPolynomial.FindPrimitive(2, 64));
        Assert.AreEqual(64, PrimeFieldPolynomial.FindIrreducible(2, 64).Degree);
    }

    #endregion

    #region ExtensionField

    [TestMethod]
    public void ExtensionField_Aes_KnownProducts()
    {
        var aes = new ExtensionField(PrimeFieldPolynomial.FromInteger(0x11B, 2));
        Assert.IsTrue(aes.UsesTables);
        Assert.AreEqual("GF(2^8)", aes.ToString());

        // FIPS 197, section 4.2.
        Assert.AreEqual(0xC1, (aes.Element(0x57) * aes.Element(0x83)).ToInt64());
        Assert.AreEqual(0xFE, (aes.Element(0x57) * aes.Element(0x13)).ToInt64());
        Assert.AreEqual(0xCA, aes.Element(0x53).Inverse().ToInt64());
        Assert.AreEqual(0xD4, (aes.Element(0x57) + aes.Element(0x83)).ToInt64()); // XOR

        Assert.AreEqual("x^6 + x^4 + x^2 + x + 1 in GF(2^8)", aes.Element(0x57).ToString());
        Assert.AreEqual(51, aes.MultiplicativeOrder(aes.Root));
        Assert.AreEqual(3, aes.PrimitiveElement.ToInt64());                      // x + 1

        // Every product agrees with shift-and-xor multiplication modulo 0x11B.
        for (int a = 0; a < 256; a++)
            for (int b = 0; b < 256; b++)
                Assert.AreEqual(AesMultiply(a, b), (aes.Element(a) * aes.Element(b)).ToInt64(), $"{a:X2}·{b:X2}");
    }

    [TestMethod]
    public void ExtensionField_SatisfiesFieldAxioms()
    {
        AssertFieldAxioms(new ExtensionField(2, 1));
        AssertFieldAxioms(new ExtensionField(2, 3));
        AssertFieldAxioms(new ExtensionField(3, 2));
        AssertFieldAxioms(new ExtensionField(5, 2));
        AssertFieldAxioms(new ExtensionField(PrimeFieldPolynomial.FromInteger(0x11B, 2))); // root not primitive
    }

    [TestMethod]
    public void ExtensionField_MultiplicativeGroupIsCyclicOfOrderQMinusOne()
    {
        foreach (var field in new[] { new ExtensionField(2, 4), new ExtensionField(3, 3), new ExtensionField(7, 2) })
        {
            long q = (long)field.Order;
            var g = field.PrimitiveElement;
            Assert.AreEqual(field.Root, g, "the default polynomial is primitive");
            Assert.AreEqual(q - 1, field.MultiplicativeOrder(g));

            // The powers of g run through every nonzero element once.
            var powers = Enumerable.Range(0, (int)q - 1).Select(k => field.Pow(g, k).ToInt64()).ToList();
            CollectionAssert.AreEquivalent(Enumerable.Range(1, (int)q - 1).Select(i => (long)i).ToList(), powers);

            // Lagrange: every order divides q − 1, and the number of generators is φ(q − 1).
            var nonzero = field.Elements().Skip(1).ToList();
            Assert.IsTrue(nonzero.All(a => (q - 1) % field.MultiplicativeOrder(a) == 0));
            Assert.AreEqual(ModularArithmetic.Totient(q - 1), nonzero.Count(a => field.MultiplicativeOrder(a) == q - 1));
        }
    }

    [TestMethod]
    public void ExtensionField_TablesMatchPolynomialArithmetic()
    {
        var rng = new Random(8);
        foreach (var field in new[] { new ExtensionField(2, 16), new ExtensionField(3, 7), new ExtensionField(251, 2) })
        {
            Assert.IsTrue(field.UsesTables, field.ToString());
            var f = field.DefiningPolynomial;
            for (int trial = 0; trial < 500; trial++)
            {
                var a = RandomElement(field, rng);
                var b = RandomElement(field, rng);
                Assert.AreEqual(a.ToPolynomial() * b.ToPolynomial() % f, (a * b).ToPolynomial());
                if (!a.IsZero)
                {
                    Assert.AreEqual(field.One, a * a.Inverse());
                    long e = rng.Next(-1000, 1000);
                    var x = e >= 0 ? a : a.Inverse();
                    Assert.AreEqual(x.ToPolynomial().PowMod(Math.Abs(e), f), a.Pow(e).ToPolynomial(), $"e = {e}");
                }
            }
        }
    }

    [TestMethod]
    public void ExtensionField_WithoutTables()
    {
        var rng = new Random(9);
        var gf2_17 = new ExtensionField(2, 17);
        var gf3_11 = new ExtensionField(3, 11);
        foreach (var field in new[] { gf2_17, gf3_11 })
        {
            Assert.IsFalse(field.UsesTables);
            long q = (long)field.Order;
            for (int trial = 0; trial < 100; trial++)
            {
                var a = RandomElement(field, rng);
                var b = RandomElement(field, rng);
                var c = RandomElement(field, rng);
                Assert.AreEqual(a * (b + c), a * b + a * c);
                Assert.AreEqual(a * b, b * a);
                if (a.IsZero) continue;
                Assert.AreEqual(field.One, a * a.Inverse());
                Assert.AreEqual(field.One, a.Pow(q - 1));
                Assert.AreEqual(a, a.Pow(q));                                      // aᵠ = a
                Assert.AreEqual(a.Inverse(), a.Pow(long.MinValue) * a.Pow(long.MaxValue) * a.Pow(-1) * a);
                Assert.AreEqual(b, b * a / a);
            }
        }
        // 2¹⁷ − 1 is a Mersenne prime, so every element except 0 and 1 generates the group.
        Assert.AreEqual((1L << 17) - 1, gf2_17.MultiplicativeOrder(gf2_17.Element(12345)));
    }

    [TestMethod]
    public void ExtensionField_LargeCharacteristic()
    {
        // GF(p²) for p = 2⁶¹ − 1: p² − 1 overflows a long, so it cannot be factored. The field uses an
        // irreducible (not necessarily primitive) polynomial, and PrimitiveElement is unavailable.
        var field = new ExtensionField(Mersenne61, 2);
        Assert.AreEqual(BigInteger.Pow(Mersenne61, 2), field.Order);
        Assert.ThrowsException<OverflowException>(() => field.PrimitiveElement);

        var rng = new Random(10);
        for (int trial = 0; trial < 20; trial++)
        {
            var a = RandomElement(field, rng);
            var b = RandomElement(field, rng);
            Assert.AreEqual((a + b) * (a - b), a * a - b * b);
            if (a.IsZero) continue;
            Assert.AreEqual(field.One, a / a);

            // Frobenius: a ↦ aᵖ is an automorphism of order 2, and it fixes exactly GF(p).
            var frobenius = a.Pow(Mersenne61);
            Assert.AreEqual(a, frobenius.Pow(Mersenne61));
            Assert.AreEqual(frobenius + b.Pow(Mersenne61), (a + b).Pow(Mersenne61));
        }
        var constant = field.FromCoefficients(new long[] { 12345 });
        Assert.AreEqual(constant, constant.Pow(Mersenne61));
    }

    [TestMethod]
    public void ExtensionField_ConstructionAndValidation()
    {
        var gf9 = new ExtensionField(3, 2);
        Assert.AreEqual(9, gf9.Elements().Count());
        Assert.AreEqual(9, gf9.Elements().Select(e => e.ToInt64()).Distinct().Count());
        Assert.AreEqual("GF(3^2)", gf9.ToString());
        Assert.AreEqual(3, gf9.Characteristic);
        Assert.AreEqual(2, gf9.Degree);

        // α is a root of the defining polynomial: f(α) = 0, evaluated with Horner in the field.
        var alpha = gf9.Root;
        var value = gf9.Zero;
        foreach (long c in gf9.DefiningPolynomial.Coefficients.Reverse())
            value = value * alpha + gf9.FromCoefficients(new[] { c });
        Assert.IsTrue(value.IsZero);

        // Coefficients and polynomials are reduced modulo p and modulo f.
        Assert.AreEqual(gf9.Element(5), gf9.FromCoefficients(new long[] { -1, 4 }));  // 2 + x → index 5
        Assert.AreEqual(alpha * alpha, gf9.FromPolynomial(PrimeFieldPolynomial.Monomial(1, 2, 3)));
        CollectionAssert.AreEqual(new long[] { 2, 1 }, gf9.Element(5).Coefficients.ToArray());

        // A field equals another with the same polynomial, and elements of equal fields mix.
        var same = new ExtensionField(gf9.DefiningPolynomial);
        Assert.AreEqual(gf9, same);
        Assert.AreEqual(gf9.Element(4), gf9.Element(1) + same.Element(3));

        Assert.ThrowsException<ArgumentException>(() => new ExtensionField(Poly(3, 1, 0, 1) * Poly(3, 1, 1)));
        Assert.ThrowsException<ArgumentException>(() => new ExtensionField(Poly(2, 1, 0, 1)));   // x² + 1 = (x + 1)²
        Assert.ThrowsException<ArgumentException>(() => gf9.Element(1) + new ExtensionField(2, 3).One);
        Assert.ThrowsException<ArgumentException>(() => gf9.Element(1) + default(ExtensionFieldElement));
        Assert.ThrowsException<ArgumentException>(() => default(ExtensionFieldElement).Inverse());
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => gf9.Element(9));
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => new ExtensionField(2, 0));
        Assert.ThrowsException<ArgumentException>(() => new ExtensionField(4, 2));
        Assert.ThrowsException<DivideByZeroException>(() => gf9.One / gf9.Zero);
        Assert.ThrowsException<DivideByZeroException>(() => gf9.Zero.Pow(-2));
        Assert.AreEqual(gf9.One, gf9.Zero.Pow(0));
        Assert.AreEqual(gf9.Zero, gf9.Zero.Pow(5));
    }

    #endregion

    #region Helpers

    /// <summary>Exhaustive field axioms over every pair, and a sample of triples.</summary>
    private static void AssertFieldAxioms<T>(IFiniteField<T> field)
    {
        var elements = field.Elements().ToList();
        Assert.AreEqual(field.Order, new BigInteger(elements.Count), "Elements() count");
        Assert.AreEqual(elements.Count, elements.Distinct().Count(), "Elements() distinct");
        Assert.AreEqual(field.Zero, elements[0]);

        var rng = new Random(4);
        var zero = field.Zero;
        var one = field.One;
        long q = elements.Count;
        int stride = Math.Max(1, elements.Count / 64); // keeps GF(2⁸) pairs quick

        for (int i = 0; i < elements.Count; i++)
        {
            var a = elements[i];
            Assert.AreEqual(a, field.Add(a, zero));
            Assert.AreEqual(a, field.Multiply(a, one));
            Assert.AreEqual(zero, field.Add(a, field.Negate(a)));
            Assert.AreEqual(zero, field.Multiply(a, zero));

            if (!a.Equals(zero))
            {
                Assert.AreEqual(one, field.Multiply(a, field.Inverse(a)), $"inverse of {a}");
                Assert.AreEqual(one, field.Pow(a, q - 1), $"Fermat for {a}");
                Assert.AreEqual(field.Inverse(a), field.Pow(a, -1));
            }

            for (int j = 0; j < elements.Count; j += stride)
            {
                var b = elements[j];
                Assert.AreEqual(field.Add(a, b), field.Add(b, a));
                Assert.AreEqual(field.Multiply(a, b), field.Multiply(b, a));
                Assert.AreEqual(a, field.Add(field.Subtract(a, b), b));
                if (!b.Equals(zero)) Assert.AreEqual(a, field.Multiply(field.Divide(a, b), b));

                var c = elements[rng.Next(elements.Count)];
                Assert.AreEqual(field.Multiply(field.Multiply(a, b), c), field.Multiply(a, field.Multiply(b, c)));
                Assert.AreEqual(field.Add(field.Add(a, b), c), field.Add(a, field.Add(b, c)));
                Assert.AreEqual(field.Multiply(a, field.Add(b, c)), field.Add(field.Multiply(a, b), field.Multiply(a, c)));
            }
        }

        // The primitive element generates every nonzero element.
        var g = field.PrimitiveElement;
        var powers = new HashSet<T>();
        var power = one;
        for (long k = 0; k < q - 1; k++, power = field.Multiply(power, g))
            powers.Add(power);
        Assert.AreEqual(q - 1, powers.Count, "primitive element");
    }

    private static PrimeFieldPolynomial Poly(long p, params long[] coefficients) => new PrimeFieldPolynomial(coefficients, p);

    private static PrimeFieldPolynomial RandomPolynomial(long p, int degree, Random rng)
    {
        var c = new long[degree + 1];
        for (int i = 0; i <= degree; i++) c[i] = rng.NextInt64(p);
        return new PrimeFieldPolynomial(c, p);
    }

    private static ExtensionFieldElement RandomElement(ExtensionField field, Random rng)
    {
        var c = new long[field.Degree];
        for (int i = 0; i < c.Length; i++) c[i] = rng.NextInt64(field.Characteristic);
        return field.FromCoefficients(c);
    }

    private static IEnumerable<PrimeFieldPolynomial> MonicPolynomials(long p, int degree)
    {
        long count = (long)Math.Pow(p, degree);
        for (long lower = 0; lower < count; lower++)
        {
            var c = new long[degree + 1];
            c[degree] = 1;
            long v = lower;
            for (int i = 0; i < degree; i++, v /= p) c[i] = v % p;
            yield return new PrimeFieldPolynomial(c, p);
        }
    }

    /// <summary>GF(2⁸) multiplication by shift-and-xor modulo 0x11B, as in FIPS 197.</summary>
    private static long AesMultiply(int a, int b)
    {
        int product = 0;
        for (; b > 0; b >>= 1)
        {
            if ((b & 1) != 0) product ^= a;
            a <<= 1;
            if ((a & 0x100) != 0) a ^= 0x11B;
        }
        return product;
    }

    #endregion
}
