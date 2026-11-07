using System;
using System.Collections.Generic;
using System.Numerics;

namespace CSharpNumerics.Numerics.NumberTheory;

/// <summary>
/// The number theoretic transform (NTT): the discrete Fourier transform over GF(p) instead of ℂ,
/// and exact integer convolution built on it.
/// </summary>
/// <remarks>
/// <para>
/// For a prime p and a length n = 2ᵏ dividing p − 1, the forward transform is
/// Xₖ = Σⱼ xⱼ·ωʲᵏ (mod p), with ω = g^((p − 1)/n) for the smallest primitive root g modulo p. The
/// inverse uses ω⁻¹ and multiplies by n⁻¹. Primes of the form c·2ᵏ + 1, such as
/// 998244353 = 119·2²³ + 1, allow long power-of-two lengths.
/// </para>
/// <para>
/// Unlike the floating-point FFT, the arithmetic is exact, so a convolution has no rounding error.
/// <see cref="ExactConvolution"/> runs the NTT modulo several primes below 2³¹ and joins the results
/// with the Chinese remainder theorem (Garner's algorithm). It uses as many primes as the size of the
/// inputs requires.
/// </para>
/// </remarks>
/// <example>
/// <code>
/// // (1 + 2x + 3x²)(4 + 5x) = 4 + 13x + 22x² + 15x³
/// long[] c = NumberTheoreticTransform.ExactConvolution(new long[] { 1, 2, 3 }, new long[] { 4, 5 });
/// long[] spectrum = NumberTheoreticTransform.Forward(new long[] { 1, 2, 3, 4 }, 998244353);
/// long[] back = NumberTheoreticTransform.Inverse(spectrum, 998244353); // 1, 2, 3, 4
/// </code>
/// </example>
public static class NumberTheoreticTransform
{
    /// <summary>998244353 = 119·2²³ + 1, the usual NTT prime: transform lengths up to 2²³.</summary>
    public const long DefaultModulus = 998244353;

    // NTT-friendly primes c·2ᵏ + 1 below 2³¹ (so MulMod multiplies directly), by decreasing k.
    // Their product is about 2²³⁸.
    private static readonly long[] ConvolutionPrimes =
    {
        2013265921, // 15·2²⁷ + 1
        469762049,  // 7·2²⁶ + 1
        1811939329, // 27·2²⁶ + 1
        167772161,  // 5·2²⁵ + 1
        2113929217, // 63·2²⁵ + 1
        754974721,  // 45·2²⁴ + 1
        998244353,  // 119·2²³ + 1
        1004535809, // 479·2²¹ + 1
    };

    /// <summary>
    /// The forward transform of <paramref name="values"/> modulo <paramref name="modulus"/>. The
    /// values are reduced modulo p first.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="values"/> is null.</exception>
    /// <exception cref="ArgumentException">The modulus is not prime, the length is not a power of two, or the length does not divide p − 1.</exception>
    public static long[] Forward(IReadOnlyList<long> values, long modulus = DefaultModulus) => Transform(values, modulus, inverse: false);

    /// <summary>The inverse of <see cref="Forward"/>: Inverse(Forward(x)) = x mod p.</summary>
    /// <exception cref="ArgumentNullException"><paramref name="values"/> is null.</exception>
    /// <exception cref="ArgumentException">The modulus is not prime, the length is not a power of two, or the length does not divide p − 1.</exception>
    public static long[] Inverse(IReadOnlyList<long> values, long modulus = DefaultModulus) => Transform(values, modulus, inverse: true);

    /// <summary>
    /// The linear convolution cₖ = Σ aᵢ·bₖ₋ᵢ reduced modulo the prime <paramref name="modulus"/>,
    /// of length |a| + |b| − 1 (empty if either input is empty).
    /// </summary>
    /// <exception cref="ArgumentException">The modulus is not prime, or the padded length 2ᵏ ≥ |a| + |b| − 1 does not divide p − 1.</exception>
    public static long[] Convolution(IReadOnlyList<long> a, IReadOnlyList<long> b, long modulus = DefaultModulus)
    {
        if (a == null) throw new ArgumentNullException(nameof(a));
        if (b == null) throw new ArgumentNullException(nameof(b));
        ValidateModulus(modulus);
        if (a.Count == 0 || b.Count == 0) return new long[0];

        int resultLength = a.Count + b.Count - 1;
        int n = PaddedLength(resultLength);
        ValidateLength(n, modulus);
        return ConvolveModulo(a, b, modulus, n, resultLength);
    }

    /// <summary>
    /// The exact linear convolution cₖ = Σ aᵢ·bₖ₋ᵢ of two integer sequences, with no rounding error
    /// and no overflow in the intermediate sums. Equivalent to multiplying the polynomials with
    /// coefficients a and b.
    /// </summary>
    /// <exception cref="ArgumentException">The inputs are too long: the result may have at most 2²¹ coefficients when the values are large, and up to 2²⁷ when they are small.</exception>
    /// <exception cref="OverflowException">A coefficient of the result does not fit in a <see cref="long"/>.</exception>
    public static long[] ExactConvolution(IReadOnlyList<long> a, IReadOnlyList<long> b)
    {
        if (a == null) throw new ArgumentNullException(nameof(a));
        if (b == null) throw new ArgumentNullException(nameof(b));
        if (a.Count == 0 || b.Count == 0) return new long[0];

        int resultLength = a.Count + b.Count - 1;
        int n = PaddedLength(resultLength);

        // |cₖ| ≤ min(|a|, |b|)·max|aᵢ|·max|bⱼ|. The primes must cover twice that, for the sign.
        BigInteger bound = Math.Min(a.Count, b.Count) * MaxAbs(a) * MaxAbs(b);
        if (bound.IsZero) return new long[resultLength];

        int count = 0;
        BigInteger product = BigInteger.One;
        while (product <= 2 * bound)
            product *= ConvolutionPrimes[count++]; // the product of all eight exceeds 2·2²¹·2⁶³·2⁶³

        long[] primes = new long[count];
        Array.Copy(ConvolutionPrimes, primes, count);
        foreach (long p in primes)
            if ((p - 1) % n != 0)
                throw new ArgumentException($"The result has {resultLength} coefficients; with values this large at most {MaxLength(primes)} are supported.");

        var residues = new long[count][];
        for (int i = 0; i < count; i++)
            residues[i] = ConvolveModulo(a, b, primes[i], n, resultLength);

        return Garner(residues, primes, product, resultLength);
    }

    private static long[] ConvolveModulo(IReadOnlyList<long> a, IReadOnlyList<long> b, long p, int n, int resultLength)
    {
        var fa = Pad(a, n, p);
        var fb = Pad(b, n, p);
        TransformInPlace(fa, p, inverse: false);
        TransformInPlace(fb, p, inverse: false);
        for (int i = 0; i < n; i++)
            fa[i] = ModularArithmetic.MulMod(fa[i], fb[i], p);
        TransformInPlace(fa, p, inverse: true);

        var result = new long[resultLength];
        Array.Copy(fa, result, resultLength);
        return result;
    }

    /// <summary>
    /// Garner's algorithm: from cₖ mod pᵢ, the value in (−M/2, M/2] for M = Πpᵢ, built as
    /// v₀ + v₁·p₀ + v₂·p₀p₁ + … with each digit vᵢ in [0, pᵢ).
    /// </summary>
    private static long[] Garner(long[][] residues, long[] primes, BigInteger product, int length)
    {
        int count = primes.Length;

        // inverses[i][j] = pⱼ⁻¹ mod pᵢ for j < i.
        var inverses = new long[count][];
        for (int i = 0; i < count; i++)
        {
            inverses[i] = new long[i];
            for (int j = 0; j < i; j++)
                inverses[i][j] = ModularArithmetic.Inverse(primes[j], primes[i]);
        }

        BigInteger half = product / 2;
        var result = new long[length];
        var digits = new long[count];
        for (int k = 0; k < length; k++)
        {
            for (int i = 0; i < count; i++)
            {
                long v = residues[i][k];
                for (int j = 0; j < i; j++)
                    v = ModularArithmetic.MulMod(ModularArithmetic.SubtractMod(v, digits[j] % primes[i], primes[i]), inverses[i][j], primes[i]);
                digits[i] = v;
            }

            BigInteger value = 0, radix = 1;
            for (int i = 0; i < count; i++)
            {
                value += digits[i] * radix;
                radix *= primes[i];
            }
            if (value > half) value -= product;

            if (value < long.MinValue || value > long.MaxValue)
                throw new OverflowException($"Coefficient {k} of the convolution is {value}, which does not fit in a long.");
            result[k] = (long)value;
        }
        return result;
    }

    private static long[] Transform(IReadOnlyList<long> values, long modulus, bool inverse)
    {
        if (values == null) throw new ArgumentNullException(nameof(values));
        ValidateModulus(modulus);
        int n = values.Count;
        if (n == 0 || (n & (n - 1)) != 0)
            throw new ArgumentException($"The length must be a power of two, not {n}.", nameof(values));
        ValidateLength(n, modulus);

        var a = Pad(values, n, modulus);
        TransformInPlace(a, modulus, inverse);
        return a;
    }

    /// <summary>Iterative radix-2 Cooley–Tukey on a power-of-two array of reduced values.</summary>
    private static void TransformInPlace(long[] a, long p, bool inverse)
    {
        int n = a.Length;
        if (n == 1) return;

        // Bit-reversal permutation.
        for (int i = 1, j = 0; i < n; i++)
        {
            int bit = n >> 1;
            for (; (j & bit) != 0; bit >>= 1) j ^= bit;
            j ^= bit;
            if (i < j) (a[i], a[j]) = (a[j], a[i]);
        }

        long g = ModularArithmetic.PrimitiveRoot(p);
        for (int length = 2; length <= n; length <<= 1)
        {
            long w = ModularArithmetic.Pow(g, (p - 1) / length, p);
            if (inverse) w = ModularArithmetic.Inverse(w, p);
            int half = length >> 1;

            // Twiddle factors wᵏ for this stage.
            var twiddles = new long[half];
            twiddles[0] = 1;
            for (int k = 1; k < half; k++)
                twiddles[k] = ModularArithmetic.MulMod(twiddles[k - 1], w, p);

            for (int start = 0; start < n; start += length)
            {
                for (int k = 0; k < half; k++)
                {
                    long u = a[start + k];
                    long v = ModularArithmetic.MulMod(a[start + k + half], twiddles[k], p);
                    a[start + k] = ModularArithmetic.AddMod(u, v, p);
                    a[start + k + half] = ModularArithmetic.SubtractMod(u, v, p);
                }
            }
        }

        if (inverse)
        {
            long nInverse = ModularArithmetic.Inverse(n, p);
            for (int i = 0; i < n; i++)
                a[i] = ModularArithmetic.MulMod(a[i], nInverse, p);
        }
    }

    private static long[] Pad(IReadOnlyList<long> values, int n, long p)
    {
        var a = new long[n];
        for (int i = 0; i < values.Count; i++)
            a[i] = ModularArithmetic.Mod(values[i], p);
        return a;
    }

    private static int PaddedLength(int length)
    {
        int n = 1;
        while (n < length)
        {
            if (n >= 1 << 30) throw new ArgumentException($"A result of {length} coefficients is too long.");
            n <<= 1;
        }
        return n;
    }

    private static BigInteger MaxAbs(IReadOnlyList<long> values)
    {
        BigInteger max = 0;
        foreach (long v in values)
            max = BigInteger.Max(max, BigInteger.Abs(v));
        return max;
    }

    private static long MaxLength(long[] primes)
    {
        long max = long.MaxValue;
        foreach (long p in primes)
            max = Math.Min(max, (p - 1) & -(p - 1)); // largest power of two dividing p − 1
        return max;
    }

    private static void ValidateModulus(long modulus)
    {
        if (!PrimalityTesting.IsPrime(modulus))
            throw new ArgumentException($"The modulus {modulus} is not prime.", nameof(modulus));
    }

    private static void ValidateLength(int n, long modulus)
    {
        if ((modulus - 1) % n != 0)
            throw new ArgumentException($"The transform length {n} does not divide {modulus} − 1, so GF({modulus}) has no root of unity of order {n}.");
    }
}
