using System;
using System.Linq;
using System.Numerics;
using CSharpNumerics.Numerics;
using CSharpNumerics.Numerics.NumberTheory;

namespace NumericsTests;

[TestClass]
public class NumberTheoreticTransformTests
{
    private const long P998 = 998244353;               // 119·2²³ + 1
    private const long P7681 = 7681;                   // 15·2⁹ + 1
    private const long PLarge = 4179340454199820289;   // 29·2⁵⁷ + 1, above the direct-multiply limit

    [TestMethod]
    public void Forward_MatchesDefinition()
    {
        var rng = new Random(1);
        foreach (long p in new[] { P998, P7681, PLarge })
        {
            foreach (int n in new[] { 1, 2, 8, 32 })
            {
                var x = Enumerable.Range(0, n).Select(_ => rng.NextInt64(p)).ToArray();
                long w = ModularArithmetic.Pow(ModularArithmetic.PrimitiveRoot(p), (p - 1) / n, p);

                var expected = new long[n];
                for (int k = 0; k < n; k++)
                    for (int j = 0; j < n; j++)
                        expected[k] = ModularArithmetic.AddMod(expected[k],
                            ModularArithmetic.MulMod(x[j], ModularArithmetic.Pow(w, (long)j * k, p), p), p);

                CollectionAssert.AreEqual(expected, NumberTheoreticTransform.Forward(x, p), $"p = {p}, n = {n}");
            }
        }
    }

    [TestMethod]
    public void InverseUndoesForward()
    {
        var rng = new Random(2);
        foreach (long p in new[] { P998, P7681, PLarge })
        {
            for (int n = 1; n <= 512; n <<= 1)
            {
                var x = Enumerable.Range(0, n).Select(_ => rng.NextInt64(p)).ToArray();
                CollectionAssert.AreEqual(x, NumberTheoreticTransform.Inverse(NumberTheoreticTransform.Forward(x, p), p), $"p = {p}, n = {n}");
            }
        }

        // Negative inputs are reduced modulo p first.
        CollectionAssert.AreEqual(new long[] { P998 - 1, 2 }, NumberTheoreticTransform.Inverse(NumberTheoreticTransform.Forward(new long[] { -1, 2 })));
    }

    [TestMethod]
    public void Convolution_MatchesNaive()
    {
        var rng = new Random(3);
        foreach (long p in new[] { P998, P7681, PLarge })
        {
            for (int trial = 0; trial < 20; trial++)
            {
                var a = Enumerable.Range(0, rng.Next(1, 60)).Select(_ => rng.NextInt64(-p, p)).ToArray();
                var b = Enumerable.Range(0, rng.Next(1, 60)).Select(_ => rng.NextInt64(-p, p)).ToArray();
                var expected = NaiveConvolution(a, b).Select(c => (long)ModularArithmetic.Mod(c, p)).ToArray();
                CollectionAssert.AreEqual(expected, NumberTheoreticTransform.Convolution(a, b, p), $"p = {p}");
            }
        }
        Assert.AreEqual(0, NumberTheoreticTransform.Convolution(new long[0], new long[] { 1 }).Length);
    }

    [TestMethod]
    public void ExactConvolution_MatchesNaive()
    {
        // (1 + 2x + 3x²)(4 + 5x) = 4 + 13x + 22x² + 15x³
        CollectionAssert.AreEqual(new long[] { 4, 13, 22, 15 }, NumberTheoreticTransform.ExactConvolution(new long[] { 1, 2, 3 }, new long[] { 4, 5 }));

        var rng = new Random(4);
        foreach (long max in new[] { 10L, 1L << 20, 1L << 28, 3037000499 })
        {
            for (int trial = 0; trial < 10; trial++)
            {
                // Keep every coefficient within a long: length · max² < 2⁶³.
                int maxLength = (int)Math.Min(300, long.MaxValue / max / max);
                var a = Enumerable.Range(0, rng.Next(1, maxLength + 1)).Select(_ => rng.NextInt64(-max, max + 1)).ToArray();
                var b = Enumerable.Range(0, rng.Next(1, maxLength + 1)).Select(_ => rng.NextInt64(-max, max + 1)).ToArray();
                var expected = NaiveConvolution(a, b).Select(c => (long)c).ToArray();
                CollectionAssert.AreEqual(expected, NumberTheoreticTransform.ExactConvolution(a, b), $"max = {max}");
            }
        }
    }

    [TestMethod]
    public void ExactConvolution_LargeValuesAndOverflow()
    {
        // The intermediate bound is 2·2⁶²·1, far above any single prime, and the terms cancel.
        CollectionAssert.AreEqual(new[] { 1L << 62, 0, -(1L << 62) },
            NumberTheoreticTransform.ExactConvolution(new[] { 1L << 62, 1L << 62 }, new long[] { 1, -1 }));
        CollectionAssert.AreEqual(new[] { long.MinValue, long.MaxValue },
            NumberTheoreticTransform.ExactConvolution(new[] { long.MinValue, long.MaxValue }, new long[] { 1 }));

        // Products near 2¹²⁶ need several primes; their sum still fits after cancellation.
        long big = long.MaxValue;
        CollectionAssert.AreEqual(new[] { big, 0L, -big },
            NumberTheoreticTransform.ExactConvolution(new[] { big, big }, new long[] { 1, -1 }));
        var squares = NumberTheoreticTransform.ExactConvolution(new[] { 3037000499L, 1 }, new[] { 3037000499L, -1 });
        CollectionAssert.AreEqual(new[] { 3037000499L * 3037000499L, 0, -1 }, squares);

        Assert.ThrowsException<OverflowException>(() => NumberTheoreticTransform.ExactConvolution(new[] { long.MaxValue }, new long[] { 2 }));
        Assert.ThrowsException<OverflowException>(() => NumberTheoreticTransform.ExactConvolution(new[] { long.MinValue }, new long[] { -1 }));
        CollectionAssert.AreEqual(new long[3], NumberTheoreticTransform.ExactConvolution(new long[2], new long[] { 0, 0 }));
    }

    [TestMethod]
    public void ExactConvolution_LongInputs()
    {
        // All-ones sequences: cₖ = number of pairs i + j = k, a triangle.
        int n = 100_000, m = 70_000;
        var c = NumberTheoreticTransform.ExactConvolution(Enumerable.Repeat(1L, n).ToArray(), Enumerable.Repeat(1L, m).ToArray());
        Assert.AreEqual(n + m - 1, c.Length);
        for (int k = 0; k < c.Length; k++)
            Assert.AreEqual(Math.Min(Math.Min(k + 1, m), n + m - 1 - k), c[k], $"k = {k}");
    }

    [TestMethod]
    public void ExtensionsAndValidation()
    {
        long[] x = { 1, 2, 3, 4 };
        CollectionAssert.AreEqual(NumberTheoreticTransform.Forward(x), x.NumberTheoreticTransform());
        CollectionAssert.AreEqual(x, x.NumberTheoreticTransform().InverseNumberTheoreticTransform());
        CollectionAssert.AreEqual(new long[] { 4, 13, 22, 15 }, new long[] { 1, 2, 3 }.ExactConvolution(new long[] { 4, 5 }));

        Assert.ThrowsException<ArgumentException>(() => NumberTheoreticTransform.Forward(new long[3]));            // not a power of two
        Assert.ThrowsException<ArgumentException>(() => NumberTheoreticTransform.Forward(new long[0]));
        Assert.ThrowsException<ArgumentException>(() => NumberTheoreticTransform.Forward(new long[4], 15));        // not prime
        Assert.ThrowsException<ArgumentException>(() => NumberTheoreticTransform.Forward(new long[1024], P7681));  // 1024 ∤ 7680
        Assert.ThrowsException<ArgumentException>(() => NumberTheoreticTransform.Convolution(new long[300], new long[300], P7681));
        Assert.ThrowsException<ArgumentNullException>(() => NumberTheoreticTransform.ExactConvolution(null!, x));
    }

    private static BigInteger[] NaiveConvolution(long[] a, long[] b)
    {
        var c = new BigInteger[a.Length + b.Length - 1];
        for (int i = 0; i < a.Length; i++)
            for (int j = 0; j < b.Length; j++)
                c[i + j] += (BigInteger)a[i] * b[j];
        return c;
    }
}
