using CSharpNumerics.Numerics.NumberTheory.Interfaces;
using System;
using System.Collections.Generic;
using System.Numerics;

namespace CSharpNumerics.Numerics.NumberTheory;

/// <summary>
/// The finite field GF(pⁿ) = GF(p)[x] / (f), where f is a monic irreducible polynomial of degree n.
/// Elements are polynomials of degree &lt; n, see <see cref="ExtensionFieldElement"/>.
/// </summary>
/// <remarks>
/// <para>
/// Addition is coefficient-wise modulo p. Multiplication multiplies the polynomials and reduces
/// modulo f, and the inverse comes from the extended Euclidean algorithm.
/// </para>
/// <para>
/// When the field has at most <see cref="MaxTableOrder"/> elements (for example GF(2⁸)),
/// multiplication, inversion and powers use precomputed log/antilog tables over a primitive
/// element instead, so each costs a few table lookups. The results are the same either way.
/// </para>
/// <para>
/// Element indices follow the usual bit notation: the coefficients are the base-p digits of the
/// index, constant term first. In GF(2⁸) with the AES polynomial, <c>Element(0x57)</c> is
/// x⁶ + x⁴ + x² + x + 1.
/// </para>
/// <para>
/// Not constant-time: the table lookups and branches depend on the values, so this is not
/// suitable for cryptographic use.
/// </para>
/// </remarks>
/// <example>
/// <code>
/// var aes = new ExtensionField(PrimeFieldPolynomial.FromInteger(0x11B, 2)); // x⁸ + x⁴ + x³ + x + 1
/// var product = aes.Element(0x57) * aes.Element(0x83);   // 0xC1
/// var inverse = aes.Element(0x53).Inverse();             // 0xCA
/// var gf9 = new ExtensionField(3, 2);                    // GF(3²), polynomial found automatically
/// </code>
/// </example>
public sealed class ExtensionField : IFiniteField<ExtensionFieldElement>, IEquatable<ExtensionField>
{
    /// <summary>Fields with at most this many elements use log/antilog tables.</summary>
    public const int MaxTableOrder = 1 << 16;

    private readonly long[] _modulus;     // monic defining polynomial, constant term first
    private readonly long _groupOrder;    // pⁿ − 1, or 0 when it does not fit in a long
    private readonly int[] _exp;          // _exp[i] = index of gⁱ for 0 ≤ i < 2(q − 1); null without tables
    private readonly int[] _log;          // _log[index] = i with gⁱ = element; entry 0 unused
    private IReadOnlyList<(long Prime, int Exponent)> _groupOrderFactors;
    private ExtensionFieldElement? _primitiveElement;

    /// <summary>
    /// Creates GF(<paramref name="prime"/>^<paramref name="degree"/>) with the first primitive
    /// polynomial of <see cref="PrimeFieldPolynomial.FindPrimitive"/>, so <see cref="Root"/> is a
    /// primitive element. When pⁿ − 1 does not fit in a <see cref="long"/> it uses the first
    /// irreducible polynomial instead.
    /// </summary>
    /// <exception cref="ArgumentException"><paramref name="prime"/> is not prime.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="degree"/> &lt; 1.</exception>
    public ExtensionField(long prime, int degree) : this(DefaultPolynomial(prime, degree)) { }

    /// <summary>Creates GF(p)[x] / (f) for an irreducible polynomial f, made monic.</summary>
    /// <exception cref="ArgumentNullException"><paramref name="definingPolynomial"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="definingPolynomial"/> is not irreducible.</exception>
    public ExtensionField(PrimeFieldPolynomial definingPolynomial)
    {
        if (definingPolynomial is null) throw new ArgumentNullException(nameof(definingPolynomial));
        if (!definingPolynomial.IsIrreducible())
            throw new ArgumentException($"{definingPolynomial} is not irreducible.", nameof(definingPolynomial));

        DefiningPolynomial = definingPolynomial.Monic();
        _modulus = DefiningPolynomial.RawCoefficients;
        Characteristic = definingPolynomial.Prime;
        Degree = definingPolynomial.Degree;
        Order = BigInteger.Pow(Characteristic, Degree);

        try { _groupOrder = PrimeFieldPolynomial.GroupOrder(Characteristic, Degree); }
        catch (OverflowException) { _groupOrder = 0; }

        if (Order <= MaxTableOrder) (_exp, _log) = BuildTables();
    }

    /// <summary>The monic irreducible polynomial f that defines the field.</summary>
    public PrimeFieldPolynomial DefiningPolynomial { get; }

    public long Characteristic { get; }

    public int Degree { get; }

    public BigInteger Order { get; }

    public ExtensionFieldElement Zero => new ExtensionFieldElement(this, new long[Degree]);

    public ExtensionFieldElement One => new ExtensionFieldElement(this, Reduce(new long[] { 1 }));

    /// <summary>α, the class of x: a root of <see cref="DefiningPolynomial"/>.</summary>
    public ExtensionFieldElement Root => new ExtensionFieldElement(this, Reduce(new long[] { 0, 1 }));

    /// <summary>
    /// A generator of the multiplicative group: <see cref="Root"/> when the defining polynomial is
    /// primitive, otherwise the nonzero element with the smallest index whose order is pⁿ − 1.
    /// </summary>
    /// <exception cref="OverflowException">pⁿ − 1 does not fit in a <see cref="long"/>.</exception>
    public ExtensionFieldElement PrimitiveElement => _primitiveElement ??= FindPrimitiveElement();

    /// <summary>True when the field uses log/antilog tables (at most <see cref="MaxTableOrder"/> elements).</summary>
    public bool UsesTables => _exp != null;

    #region Elements

    /// <summary>The element whose coefficients are the base-p digits of <paramref name="index"/>, constant term first.</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="index"/> is outside [0, pⁿ).</exception>
    public ExtensionFieldElement Element(long index)
    {
        if (index < 0 || index >= Order)
            throw new ArgumentOutOfRangeException(nameof(index), index, $"The index must be in [0, {Order}).");
        return new ExtensionFieldElement(this, FromIndex(index));
    }

    /// <summary>The element with the given coefficients (constant term first), reduced modulo p and modulo f.</summary>
    /// <exception cref="ArgumentNullException"><paramref name="coefficients"/> is null.</exception>
    public ExtensionFieldElement FromCoefficients(IReadOnlyList<long> coefficients)
    {
        if (coefficients == null) throw new ArgumentNullException(nameof(coefficients));
        var c = new long[coefficients.Count];
        for (int i = 0; i < c.Length; i++)
            c[i] = ModularArithmetic.Mod(coefficients[i], Characteristic);
        return new ExtensionFieldElement(this, Reduce(c));
    }

    /// <summary>The class of <paramref name="polynomial"/> modulo f.</summary>
    /// <exception cref="ArgumentException">The polynomial is over another prime.</exception>
    public ExtensionFieldElement FromPolynomial(PrimeFieldPolynomial polynomial)
    {
        if (polynomial is null) throw new ArgumentNullException(nameof(polynomial));
        if (polynomial.Prime != Characteristic)
            throw new ArgumentException($"{polynomial} is not over GF({Characteristic}).", nameof(polynomial));
        return new ExtensionFieldElement(this, Reduce(polynomial.RawCoefficients));
    }

    public IEnumerable<ExtensionFieldElement> Elements()
    {
        // Base-p counter over the coefficients, constant term as the least significant digit.
        var c = new long[Degree];
        while (true)
        {
            yield return new ExtensionFieldElement(this, (long[])c.Clone());
            int i = 0;
            while (i < Degree && ++c[i] == Characteristic) c[i++] = 0;
            if (i == Degree) yield break;
        }
    }

    #endregion

    #region Arithmetic

    public ExtensionFieldElement Add(ExtensionFieldElement a, ExtensionFieldElement b)
    {
        long[] x = Check(a), y = Check(b);
        var c = new long[Degree];
        for (int i = 0; i < c.Length; i++)
            c[i] = ModularArithmetic.AddMod(x[i], y[i], Characteristic);
        return new ExtensionFieldElement(this, c);
    }

    public ExtensionFieldElement Subtract(ExtensionFieldElement a, ExtensionFieldElement b)
    {
        long[] x = Check(a), y = Check(b);
        var c = new long[Degree];
        for (int i = 0; i < c.Length; i++)
            c[i] = ModularArithmetic.SubtractMod(x[i], y[i], Characteristic);
        return new ExtensionFieldElement(this, c);
    }

    public ExtensionFieldElement Negate(ExtensionFieldElement a)
    {
        long[] x = Check(a);
        var c = new long[Degree];
        for (int i = 0; i < c.Length; i++)
            c[i] = x[i] == 0 ? 0 : Characteristic - x[i];
        return new ExtensionFieldElement(this, c);
    }

    public ExtensionFieldElement Multiply(ExtensionFieldElement a, ExtensionFieldElement b)
    {
        long[] x = Check(a), y = Check(b);
        if (_exp == null) return new ExtensionFieldElement(this, MultiplyRaw(x, y));

        int i = (int)Index(x), j = (int)Index(y);
        if (i == 0 || j == 0) return Zero;
        return new ExtensionFieldElement(this, FromIndex(_exp[_log[i] + _log[j]]));
    }

    public ExtensionFieldElement Divide(ExtensionFieldElement a, ExtensionFieldElement b) => Multiply(a, Inverse(b));

    public ExtensionFieldElement Inverse(ExtensionFieldElement a)
    {
        long[] x = Check(a);
        if (IsZero(x)) throw new DivideByZeroException("Zero has no inverse.");
        return new ExtensionFieldElement(this, InverseRaw(x));
    }

    public ExtensionFieldElement Pow(ExtensionFieldElement a, long exponent)
    {
        long[] x = Check(a);
        if (exponent == 0) return One;
        if (IsZero(x))
        {
            if (exponent < 0) throw new DivideByZeroException("Zero has no inverse.");
            return Zero;
        }

        if (_exp != null)
        {
            int groupOrder = (int)_groupOrder;
            long k = ModularArithmetic.Mod(exponent, groupOrder);
            return new ExtensionFieldElement(this, FromIndex(_exp[_log[Index(x)] * k % groupOrder]));
        }

        if (exponent < 0)
        {
            x = InverseRaw(x);
            if (exponent == long.MinValue)
                return new ExtensionFieldElement(this, MultiplyRaw(PrimeFieldPolynomial.PowModRaw(x, long.MaxValue, _modulus, Characteristic), x));
            exponent = -exponent;
        }
        // a^(pⁿ − 1) = 1 for every nonzero a.
        if (_groupOrder != 0) exponent %= _groupOrder;
        return new ExtensionFieldElement(this, PrimeFieldPolynomial.PowModRaw(x, exponent, _modulus, Characteristic));
    }

    /// <summary>The order of <paramref name="a"/> in the multiplicative group: the smallest k ≥ 1 with aᵏ = 1. It divides pⁿ − 1.</summary>
    /// <exception cref="ArgumentException"><paramref name="a"/> is zero.</exception>
    /// <exception cref="OverflowException">pⁿ − 1 does not fit in a <see cref="long"/>.</exception>
    public long MultiplicativeOrder(ExtensionFieldElement a)
    {
        long[] x = Check(a);
        if (IsZero(x)) throw new ArgumentException("Zero has no multiplicative order.", nameof(a));
        long order = RequireGroupOrder();

        if (_exp != null) return order / EuclideanAlgorithm.Gcd(_log[Index(x)], order);

        _groupOrderFactors ??= IntegerFactorization.Factor(order);
        foreach (var (q, e) in _groupOrderFactors)
        {
            for (int i = 0; i < e && Pow(a, order / q) == One; i++)
                order /= q;
        }
        return order;
    }

    #endregion

    #region Equality and formatting

    public bool Equals(ExtensionField other) =>
        ReferenceEquals(this, other) || (other is not null && DefiningPolynomial == other.DefiningPolynomial);

    public override bool Equals(object obj) => obj is ExtensionField other && Equals(other);

    public override int GetHashCode() => DefiningPolynomial.GetHashCode();

    /// <summary>Formats as "GF(2^8)", or "GF(7)" for degree 1.</summary>
    public override string ToString() => Degree == 1 ? $"GF({Characteristic})" : $"GF({Characteristic}^{Degree})";

    #endregion

    #region Internals

    /// <summary>The base-p value of a coefficient array, constant term as the least significant digit.</summary>
    internal long Index(long[] c)
    {
        long index = 0;
        for (int i = c.Length - 1; i >= 0; i--)
            index = checked(index * Characteristic + c[i]);
        return index;
    }

    internal string Format(long[] c) => $"{PrimeFieldPolynomial.FormatTerms(PrimeFieldPolynomial.Trim(c), "x")} in {this}";

    private long[] FromIndex(long index)
    {
        var c = new long[Degree];
        for (int i = 0; i < c.Length && index > 0; i++, index /= Characteristic)
            c[i] = index % Characteristic;
        return c;
    }

    private long[] Reduce(long[] c) => PrimeFieldPolynomial.Reduce(c, _modulus, Characteristic);

    private long[] MultiplyRaw(long[] x, long[] y) => Reduce(PrimeFieldPolynomial.Multiply(x, y, Characteristic));

    private long[] InverseRaw(long[] x)
    {
        if (_exp != null)
        {
            int groupOrder = (int)_groupOrder;
            return FromIndex(_exp[(groupOrder - _log[Index(x)]) % groupOrder]);
        }

        // f is irreducible and x ≠ 0, so gcd(x, f) = 1 and x·s + f·t = 1 gives x⁻¹ = s.
        var (_, s, _) = PrimeFieldPolynomial.ExtendedGcd(PrimeFieldPolynomial.FromTrusted(x, Characteristic), DefiningPolynomial);
        return Reduce(s.RawCoefficients);
    }

    private static bool IsZero(long[] c)
    {
        foreach (long v in c)
            if (v != 0) return false;
        return true;
    }

    private long[] Check(ExtensionFieldElement a)
    {
        if (a.Field is null)
            throw new ArgumentException("default(ExtensionFieldElement) is not a valid element.");
        if (!Equals(a.Field))
            throw new ArgumentException($"{a} is not an element of {this}.");
        return a.Raw;
    }

    private long RequireGroupOrder()
    {
        if (_groupOrder == 0)
            throw new OverflowException($"{Characteristic}^{Degree} − 1 does not fit in a long.");
        return _groupOrder;
    }

    private ExtensionFieldElement FindPrimitiveElement()
    {
        long order = RequireGroupOrder();
        var root = Root;
        if (!root.IsZero && MultiplicativeOrder(root) == order) return root;

        for (long index = 1; ; index++)
        {
            var candidate = new ExtensionFieldElement(this, FromIndex(index));
            if (MultiplicativeOrder(candidate) == order) return candidate;
        }
    }

    private (int[] Exp, int[] Log) BuildTables()
    {
        // Runs from the constructor before the tables exist, so PrimitiveElement uses the polynomial arithmetic.
        int groupOrder = (int)_groupOrder;
        long[] generator = PrimitiveElement.Raw;

        var exp = new int[2 * groupOrder];
        var log = new int[groupOrder + 1];
        long[] power = Reduce(new long[] { 1 });
        for (int i = 0; i < groupOrder; i++)
        {
            int index = (int)Index(power);
            exp[i] = exp[i + groupOrder] = index;
            log[index] = i;
            power = MultiplyRaw(power, generator);
        }
        return (exp, log);
    }

    private static PrimeFieldPolynomial DefaultPolynomial(long prime, int degree)
    {
        if (degree < 1) throw new ArgumentOutOfRangeException(nameof(degree), degree, "The degree must be at least 1.");
        if (!PrimalityTesting.IsPrime(prime)) throw new ArgumentException($"{prime} is not prime.", nameof(prime));

        try { PrimeFieldPolynomial.GroupOrder(prime, degree); }
        catch (OverflowException) { return PrimeFieldPolynomial.FindIrreducible(prime, degree); }
        return PrimeFieldPolynomial.FindPrimitive(prime, degree);
    }

    #endregion
}
