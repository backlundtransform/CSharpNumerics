using System;
using System.Collections.Generic;
using System.Text;

namespace CSharpNumerics.Numerics.NumberTheory;

/// <summary>
/// An immutable polynomial c₀ + c₁x + … + cₙxⁿ with coefficients in the prime field GF(p).
/// </summary>
/// <remarks>
/// <para>
/// Coefficients are stored constant term first and always reduced into [0, p). The zero
/// polynomial has <see cref="Degree"/> −1. Operations on polynomials over different primes throw
/// <see cref="ArgumentException"/>.
/// </para>
/// <para>
/// Besides ring arithmetic (+, −, ×, division with remainder, gcd), the class tests
/// irreducibility (Rabin's test) and primitivity, and searches for irreducible and primitive
/// polynomials. These define the extension fields GF(pⁿ) in <see cref="ExtensionField"/>.
/// </para>
/// <para>
/// Not constant-time: do not use for cryptographic purposes.
/// </para>
/// </remarks>
/// <example>
/// <code>
/// var f = new PrimeFieldPolynomial(new long[] { 1, 1, 0, 1, 1, 0, 0, 0, 1 }, 2); // x⁸ + x⁴ + x³ + x + 1
/// bool irreducible = f.IsIrreducible();  // true: the AES polynomial
/// bool primitive = f.IsPrimitive();      // false: x has order 51, not 255
/// var g = PrimeFieldPolynomial.FindPrimitive(2, 8); // x⁸ + x⁴ + x³ + x² + 1
/// </code>
/// </example>
public sealed class PrimeFieldPolynomial : IEquatable<PrimeFieldPolynomial>
{
    private readonly long[] _coefficients; // constant term first, no trailing zeros; empty for the zero polynomial

    /// <summary>Creates a polynomial from its coefficients, constant term first. Each is reduced modulo p.</summary>
    /// <exception cref="ArgumentNullException"><paramref name="coefficients"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="prime"/> is not prime.</exception>
    public PrimeFieldPolynomial(IReadOnlyList<long> coefficients, long prime)
    {
        if (coefficients == null) throw new ArgumentNullException(nameof(coefficients));
        ValidatePrime(prime);

        var c = new long[coefficients.Count];
        for (int i = 0; i < c.Length; i++)
            c[i] = ModularArithmetic.Mod(coefficients[i], prime);
        _coefficients = Trim(c);
        Prime = prime;
    }

    // Trusted constructor: coefficients already reduced and trimmed, prime already validated.
    private PrimeFieldPolynomial(long[] trimmed, long prime)
    {
        _coefficients = trimmed;
        Prime = prime;
    }

    /// <summary>The prime p of the coefficient field.</summary>
    public long Prime { get; }

    /// <summary>The degree, or −1 for the zero polynomial.</summary>
    public int Degree => _coefficients.Length - 1;

    /// <summary>True for the zero polynomial.</summary>
    public bool IsZero => _coefficients.Length == 0;

    /// <summary>The coefficient of the highest power, or 0 for the zero polynomial.</summary>
    public long LeadingCoefficient => IsZero ? 0 : _coefficients[_coefficients.Length - 1];

    /// <summary>True when the leading coefficient is 1.</summary>
    public bool IsMonic => LeadingCoefficient == 1;

    /// <summary>The coefficients c₀ … cₙ, constant term first. Empty for the zero polynomial.</summary>
    public IReadOnlyList<long> Coefficients => Array.AsReadOnly(_coefficients);

    /// <summary>The coefficient of xⁱ: zero above the degree.</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="power"/> is negative.</exception>
    public long this[int power]
    {
        get
        {
            if (power < 0) throw new ArgumentOutOfRangeException(nameof(power), power, "The power must be non-negative.");
            return power < _coefficients.Length ? _coefficients[power] : 0;
        }
    }

    #region Construction

    /// <summary>The zero polynomial over GF(<paramref name="prime"/>).</summary>
    public static PrimeFieldPolynomial Zero(long prime)
    {
        ValidatePrime(prime);
        return new PrimeFieldPolynomial(new long[0], prime);
    }

    /// <summary>The constant polynomial 1 over GF(<paramref name="prime"/>).</summary>
    public static PrimeFieldPolynomial One(long prime)
    {
        ValidatePrime(prime);
        return new PrimeFieldPolynomial(new long[] { 1 }, prime);
    }

    /// <summary>The monomial c·xᵈ over GF(<paramref name="prime"/>).</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="degree"/> is negative.</exception>
    public static PrimeFieldPolynomial Monomial(long coefficient, int degree, long prime)
    {
        if (degree < 0) throw new ArgumentOutOfRangeException(nameof(degree), degree, "The degree must be non-negative.");
        ValidatePrime(prime);
        var c = new long[degree + 1];
        c[degree] = ModularArithmetic.Mod(coefficient, prime);
        return new PrimeFieldPolynomial(Trim(c), prime);
    }

    /// <summary>
    /// The polynomial whose coefficients are the base-p digits of <paramref name="value"/>, least
    /// significant digit as the constant term. Over GF(2) this is the usual bit notation:
    /// 0x11B is x⁸ + x⁴ + x³ + x + 1.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="value"/> is negative.</exception>
    /// <exception cref="ArgumentException"><paramref name="prime"/> is not prime.</exception>
    public static PrimeFieldPolynomial FromInteger(long value, long prime)
    {
        if (value < 0) throw new ArgumentOutOfRangeException(nameof(value), value, "The value must be non-negative.");
        ValidatePrime(prime);
        var c = new List<long>();
        for (; value > 0; value /= prime)
            c.Add(value % prime);
        return new PrimeFieldPolynomial(c.ToArray(), prime);
    }

    /// <summary>
    /// The base-p value of the coefficients, the inverse of <see cref="FromInteger"/>.
    /// </summary>
    /// <exception cref="OverflowException">The value does not fit in a <see cref="long"/>.</exception>
    public long ToInteger()
    {
        long value = 0;
        for (int i = _coefficients.Length - 1; i >= 0; i--)
            value = checked(value * Prime + _coefficients[i]);
        return value;
    }

    #endregion

    #region Arithmetic

    public static PrimeFieldPolynomial operator +(PrimeFieldPolynomial a, PrimeFieldPolynomial b)
    {
        long p = CommonPrime(a, b);
        var c = new long[Math.Max(a._coefficients.Length, b._coefficients.Length)];
        for (int i = 0; i < c.Length; i++)
            c[i] = ModularArithmetic.AddMod(a.Coefficient(i), b.Coefficient(i), p);
        return new PrimeFieldPolynomial(Trim(c), p);
    }

    public static PrimeFieldPolynomial operator -(PrimeFieldPolynomial a, PrimeFieldPolynomial b)
    {
        long p = CommonPrime(a, b);
        var c = new long[Math.Max(a._coefficients.Length, b._coefficients.Length)];
        for (int i = 0; i < c.Length; i++)
            c[i] = ModularArithmetic.SubtractMod(a.Coefficient(i), b.Coefficient(i), p);
        return new PrimeFieldPolynomial(Trim(c), p);
    }

    public static PrimeFieldPolynomial operator -(PrimeFieldPolynomial a)
    {
        if (a == null) throw new ArgumentNullException(nameof(a));
        var c = new long[a._coefficients.Length];
        for (int i = 0; i < c.Length; i++)
            c[i] = a._coefficients[i] == 0 ? 0 : a.Prime - a._coefficients[i];
        return new PrimeFieldPolynomial(c, a.Prime);
    }

    public static PrimeFieldPolynomial operator *(PrimeFieldPolynomial a, PrimeFieldPolynomial b)
    {
        long p = CommonPrime(a, b);
        return new PrimeFieldPolynomial(Trim(Multiply(a._coefficients, b._coefficients, p)), p);
    }

    /// <summary>Multiplies every coefficient by <paramref name="scalar"/> (reduced modulo p).</summary>
    public static PrimeFieldPolynomial operator *(long scalar, PrimeFieldPolynomial a)
    {
        if (a == null) throw new ArgumentNullException(nameof(a));
        return a.Scale(ModularArithmetic.Mod(scalar, a.Prime));
    }

    /// <summary>Multiplies every coefficient by <paramref name="scalar"/> (reduced modulo p).</summary>
    public static PrimeFieldPolynomial operator *(PrimeFieldPolynomial a, long scalar) => scalar * a;

    /// <summary>The quotient of <see cref="DivRem"/>.</summary>
    /// <exception cref="DivideByZeroException"><paramref name="b"/> is the zero polynomial.</exception>
    public static PrimeFieldPolynomial operator /(PrimeFieldPolynomial a, PrimeFieldPolynomial b) => DivRem(a, b).Quotient;

    /// <summary>The remainder of <see cref="DivRem"/>.</summary>
    /// <exception cref="DivideByZeroException"><paramref name="b"/> is the zero polynomial.</exception>
    public static PrimeFieldPolynomial operator %(PrimeFieldPolynomial a, PrimeFieldPolynomial b) => DivRem(a, b).Remainder;

    /// <summary>
    /// Division with remainder: a = q·b + r with deg r &lt; deg b. Exact over a field, since every
    /// nonzero leading coefficient is invertible.
    /// </summary>
    /// <exception cref="DivideByZeroException"><paramref name="divisor"/> is the zero polynomial.</exception>
    public static (PrimeFieldPolynomial Quotient, PrimeFieldPolynomial Remainder) DivRem(PrimeFieldPolynomial dividend, PrimeFieldPolynomial divisor)
    {
        long p = CommonPrime(dividend, divisor);
        if (divisor.IsZero) throw new DivideByZeroException("Division by the zero polynomial.");

        var quotient = new long[Math.Max(0, dividend.Degree - divisor.Degree + 1)];
        var remainder = (long[])dividend._coefficients.Clone();
        DivRemInPlace(remainder, divisor._coefficients, p, quotient);
        return (new PrimeFieldPolynomial(Trim(quotient), p), new PrimeFieldPolynomial(Trim(remainder), p));
    }

    /// <summary>The monic greatest common divisor. Gcd(0, 0) = 0.</summary>
    public static PrimeFieldPolynomial Gcd(PrimeFieldPolynomial a, PrimeFieldPolynomial b)
    {
        CommonPrime(a, b);
        while (!b.IsZero)
            (a, b) = (b, a % b);
        return a.IsZero ? a : a.Monic();
    }

    /// <summary>
    /// The extended Euclidean algorithm: the monic gcd g together with s and t such that
    /// a·s + b·t = g. For Gcd(0, 0) all three are zero.
    /// </summary>
    public static (PrimeFieldPolynomial Gcd, PrimeFieldPolynomial S, PrimeFieldPolynomial T) ExtendedGcd(PrimeFieldPolynomial a, PrimeFieldPolynomial b)
    {
        long p = CommonPrime(a, b);
        var zero = new PrimeFieldPolynomial(new long[0], p);
        var one = new PrimeFieldPolynomial(new long[] { 1 }, p);

        // Invariants: a·s0 + b·t0 = r0 and a·s1 + b·t1 = r1.
        PrimeFieldPolynomial r0 = a, r1 = b, s0 = one, s1 = zero, t0 = zero, t1 = one;
        while (!r1.IsZero)
        {
            var (q, r) = DivRem(r0, r1);
            (r0, r1) = (r1, r);
            (s0, s1) = (s1, s0 - q * s1);
            (t0, t1) = (t1, t0 - q * t1);
        }

        if (r0.IsZero) return (zero, zero, zero);
        long scale = ModularArithmetic.Inverse(r0.LeadingCoefficient, p);
        return (r0.Scale(scale), s0.Scale(scale), t0.Scale(scale));
    }

    /// <summary>This polynomial divided by its leading coefficient.</summary>
    /// <exception cref="InvalidOperationException">The polynomial is zero.</exception>
    public PrimeFieldPolynomial Monic()
    {
        if (IsZero) throw new InvalidOperationException("The zero polynomial has no monic multiple.");
        return IsMonic ? this : Scale(ModularArithmetic.Inverse(LeadingCoefficient, Prime));
    }

    /// <summary>The formal derivative c₁ + 2c₂x + … + n·cₙxⁿ⁻¹.</summary>
    public PrimeFieldPolynomial Derivative()
    {
        if (_coefficients.Length <= 1) return new PrimeFieldPolynomial(new long[0], Prime);
        var c = new long[_coefficients.Length - 1];
        for (int i = 1; i < _coefficients.Length; i++)
            c[i - 1] = ModularArithmetic.MulMod(i % Prime, _coefficients[i], Prime);
        return new PrimeFieldPolynomial(Trim(c), Prime);
    }

    /// <summary>The value at <paramref name="x"/> (reduced modulo p), by Horner's rule.</summary>
    public long Evaluate(long x)
    {
        x = ModularArithmetic.Mod(x, Prime);
        long value = 0;
        for (int i = _coefficients.Length - 1; i >= 0; i--)
            value = ModularArithmetic.AddMod(ModularArithmetic.MulMod(value, x, Prime), _coefficients[i], Prime);
        return value;
    }

    /// <summary>The value at <paramref name="x"/>, an element of GF(p).</summary>
    /// <exception cref="ArgumentException">The modulus of <paramref name="x"/> is not p.</exception>
    public ModInt Evaluate(ModInt x)
    {
        if (x.Modulus != Prime) throw new ArgumentException($"{x} is not an element of GF({Prime}).", nameof(x));
        return new ModInt(Evaluate(x.Value), Prime);
    }

    /// <summary>This polynomial raised to <paramref name="exponent"/>, reduced modulo <paramref name="modulus"/>, by square-and-multiply.</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="exponent"/> is negative.</exception>
    /// <exception cref="DivideByZeroException"><paramref name="modulus"/> is the zero polynomial.</exception>
    public PrimeFieldPolynomial PowMod(long exponent, PrimeFieldPolynomial modulus)
    {
        if (exponent < 0) throw new ArgumentOutOfRangeException(nameof(exponent), exponent, "The exponent must be non-negative.");
        long p = CommonPrime(this, modulus);
        if (modulus.IsZero) throw new DivideByZeroException("Division by the zero polynomial.");
        return new PrimeFieldPolynomial(Trim(PowModRaw(_coefficients, exponent, modulus._coefficients, p)), p);
    }

    #endregion

    #region Irreducibility

    /// <summary>
    /// True when the polynomial has degree ≥ 1 and is not a product of two polynomials of lower
    /// degree. Uses Rabin's test: a polynomial f of degree n is irreducible exactly when
    /// x^(pⁿ) ≡ x (mod f) and gcd(x^(p^(n/q)) − x, f) = 1 for every prime q dividing n.
    /// </summary>
    public bool IsIrreducible()
    {
        int n = Degree;
        if (n < 1) return false;
        if (n == 1) return true;

        long[] f = Monic()._coefficients;
        var nFactors = IntegerFactorization.Factor(n);
        long[] x = { 0, 1 };

        // Frobenius powers x^(p^k) mod f for k = 1 … n.
        long[] power = x;
        for (int k = 1; k <= n; k++)
        {
            power = PowModRaw(power, Prime, f, Prime);
            if (k == n) break;

            foreach (var (q, _) in nFactors)
            {
                if (k != n / q) continue;
                var difference = new PrimeFieldPolynomial(Trim(Subtract(power, x, Prime)), Prime);
                if (Gcd(new PrimeFieldPolynomial(f, Prime), difference).Degree != 0) return false;
            }
        }
        return Trim(Subtract(power, x, Prime)).Length == 0;
    }

    /// <summary>
    /// True when the polynomial is irreducible of degree n and x generates the multiplicative group
    /// of GF(p)[x]/(f), so its order is pⁿ − 1. The powers of x then run through every nonzero
    /// element of GF(pⁿ).
    /// </summary>
    /// <exception cref="OverflowException">pⁿ − 1 does not fit in a <see cref="long"/>.</exception>
    public bool IsPrimitive()
    {
        if (Degree < 1) return false;
        long order = GroupOrder(Prime, Degree);
        return IsPrimitive(Monic(), order, IntegerFactorization.Factor(order));
    }

    /// <summary>
    /// The first monic irreducible polynomial of degree <paramref name="degree"/> over
    /// GF(<paramref name="prime"/>), with the lower coefficients enumerated as a base-p counter.
    /// The result is deterministic: over GF(2) of degree 8 it is x⁸ + x⁴ + x³ + x + 1 (0x11B).
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="degree"/> &lt; 1.</exception>
    /// <exception cref="ArgumentException"><paramref name="prime"/> is not prime.</exception>
    public static PrimeFieldPolynomial FindIrreducible(long prime, int degree) =>
        Find(prime, degree, f => f.IsIrreducible());

    /// <summary>
    /// The first monic primitive polynomial of degree <paramref name="degree"/> over
    /// GF(<paramref name="prime"/>), in the same order as <see cref="FindIrreducible"/>. Over GF(2)
    /// of degree 8 it is x⁸ + x⁴ + x³ + x² + 1 (0x11D).
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="degree"/> &lt; 1.</exception>
    /// <exception cref="ArgumentException"><paramref name="prime"/> is not prime.</exception>
    /// <exception cref="OverflowException">pⁿ − 1 does not fit in a <see cref="long"/>.</exception>
    public static PrimeFieldPolynomial FindPrimitive(long prime, int degree)
    {
        if (degree < 1) throw new ArgumentOutOfRangeException(nameof(degree), degree, "The degree must be at least 1.");
        ValidatePrime(prime);
        long order = GroupOrder(prime, degree);
        var orderFactors = IntegerFactorization.Factor(order);
        return Find(prime, degree, f => IsPrimitive(f, order, orderFactors));
    }

    /// <summary>pⁿ − 1, the order of the multiplicative group of GF(pⁿ).</summary>
    /// <exception cref="OverflowException">pⁿ − 1 does not fit in a <see cref="long"/>.</exception>
    internal static long GroupOrder(long prime, int degree)
    {
        try
        {
            // pⁿ − 1 = (p − 1)(pⁿ⁻¹ + … + p + 1) avoids overflow when pⁿ itself is 2⁶³.
            long sum = 0;
            for (int i = 0; i < degree; i++)
                sum = checked(sum * prime + 1);
            return checked((prime - 1) * sum);
        }
        catch (OverflowException)
        {
            throw new OverflowException($"{prime}^{degree} − 1 does not fit in a long.");
        }
    }

    private static bool IsPrimitive(PrimeFieldPolynomial monic, long order, IReadOnlyList<(long Prime, int Exponent)> orderFactors)
    {
        // x must be invertible modulo f, which rules out f = x.
        if (!monic.IsIrreducible() || monic._coefficients[0] == 0) return false;

        long[] x = Reduce(new long[] { 0, 1 }, monic._coefficients, monic.Prime);
        foreach (var (q, _) in orderFactors)
        {
            long[] power = Trim(PowModRaw(x, order / q, monic._coefficients, monic.Prime));
            if (power.Length == 1 && power[0] == 1) return false;
        }
        return true;
    }

    private static PrimeFieldPolynomial Find(long prime, int degree, Func<PrimeFieldPolynomial, bool> accept)
    {
        if (degree < 1) throw new ArgumentOutOfRangeException(nameof(degree), degree, "The degree must be at least 1.");
        ValidatePrime(prime);

        // Lower coefficients c₀ … cₙ₋₁ as a base-p counter, c₀ the least significant digit.
        var c = new long[degree + 1];
        c[degree] = 1;
        while (true)
        {
            // For degree ≥ 2 a zero constant term means x divides f.
            if (degree == 1 || c[0] != 0)
            {
                var candidate = new PrimeFieldPolynomial((long[])c.Clone(), prime);
                if (accept(candidate)) return candidate;
            }

            int i = 0;
            while (i < degree && ++c[i] == prime) c[i++] = 0;
            if (i == degree)
                throw new InvalidOperationException($"No polynomial of degree {degree} over GF({prime}) was accepted.");
        }
    }

    #endregion

    #region Equality and formatting

    public static bool operator ==(PrimeFieldPolynomial a, PrimeFieldPolynomial b) =>
        ReferenceEquals(a, b) || (a is not null && a.Equals(b));

    public static bool operator !=(PrimeFieldPolynomial a, PrimeFieldPolynomial b) => !(a == b);

    public bool Equals(PrimeFieldPolynomial other)
    {
        if (other is null || Prime != other.Prime || _coefficients.Length != other._coefficients.Length) return false;
        for (int i = 0; i < _coefficients.Length; i++)
            if (_coefficients[i] != other._coefficients[i]) return false;
        return true;
    }

    public override bool Equals(object obj) => obj is PrimeFieldPolynomial other && Equals(other);

    public override int GetHashCode()
    {
        var hash = new HashCode();
        hash.Add(Prime);
        foreach (long c in _coefficients) hash.Add(c);
        return hash.ToHashCode();
    }

    /// <summary>Formats as "x^2 + 2x + 1 (mod 3)".</summary>
    public override string ToString() => $"{FormatTerms(_coefficients, "x")} (mod {Prime})";

    /// <summary>Formats coefficients (constant term first) as "x^2 + 2x + 1", or "0".</summary>
    internal static string FormatTerms(long[] coefficients, string variable)
    {
        var text = new StringBuilder();
        for (int i = coefficients.Length - 1; i >= 0; i--)
        {
            long c = coefficients[i];
            if (c == 0) continue;
            if (text.Length > 0) text.Append(" + ");
            if (c != 1 || i == 0) text.Append(c);
            if (i >= 1) text.Append(variable);
            if (i >= 2) text.Append('^').Append(i);
        }
        return text.Length == 0 ? "0" : text.ToString();
    }

    #endregion

    #region Raw coefficient arrays

    // These helpers work on reduced coefficient arrays, constant term first, and are shared with
    // ExtensionField. Results may have trailing zeros unless trimmed.

    internal static PrimeFieldPolynomial FromTrusted(long[] reduced, long prime) => new PrimeFieldPolynomial(Trim(reduced), prime);

    internal long[] RawCoefficients => _coefficients;

    internal static long[] Trim(long[] c)
    {
        int length = c.Length;
        while (length > 0 && c[length - 1] == 0) length--;
        if (length == c.Length) return c;
        var trimmed = new long[length];
        Array.Copy(c, trimmed, length);
        return trimmed;
    }

    internal static long[] Multiply(long[] a, long[] b, long p)
    {
        if (a.Length == 0 || b.Length == 0) return new long[0];
        var c = new long[a.Length + b.Length - 1];
        for (int i = 0; i < a.Length; i++)
        {
            if (a[i] == 0) continue;
            for (int j = 0; j < b.Length; j++)
                c[i + j] = ModularArithmetic.AddMod(c[i + j], ModularArithmetic.MulMod(a[i], b[j], p), p);
        }
        return c;
    }

    internal static long[] Subtract(long[] a, long[] b, long p)
    {
        var c = new long[Math.Max(a.Length, b.Length)];
        for (int i = 0; i < c.Length; i++)
            c[i] = ModularArithmetic.SubtractMod(i < a.Length ? a[i] : 0, i < b.Length ? b[i] : 0, p);
        return c;
    }

    /// <summary>a mod m as an array of length deg m (zero-padded). m must be trimmed and nonzero.</summary>
    internal static long[] Reduce(long[] a, long[] m, long p)
    {
        var r = (long[])a.Clone();
        DivRemInPlace(r, m, p, null);
        var result = new long[m.Length - 1];
        Array.Copy(r, result, Math.Min(r.Length, result.Length));
        return result;
    }

    internal static long[] PowModRaw(long[] b, long exponent, long[] m, long p)
    {
        long[] result = Reduce(new long[] { 1 }, m, p);
        b = Reduce(b, m, p);
        for (; exponent > 0; exponent >>= 1)
        {
            if ((exponent & 1) == 1) result = Reduce(Multiply(result, b, p), m, p);
            if (exponent > 1) b = Reduce(Multiply(b, b, p), m, p);
        }
        return result;
    }

    /// <summary>Replaces <paramref name="r"/> by r mod d (in its low entries) and fills <paramref name="quotient"/> when given.</summary>
    private static void DivRemInPlace(long[] r, long[] d, long p, long[] quotient)
    {
        int dd = d.Length - 1;
        long leadInverse = d[dd] == 1 ? 1 : ModularArithmetic.Inverse(d[dd], p);
        for (int i = r.Length - 1; i >= dd; i--)
        {
            long c = r[i];
            if (c == 0) continue;
            if (leadInverse != 1) c = ModularArithmetic.MulMod(c, leadInverse, p);
            if (quotient != null) quotient[i - dd] = c;
            for (int j = 0; j <= dd; j++)
                r[i - dd + j] = ModularArithmetic.SubtractMod(r[i - dd + j], ModularArithmetic.MulMod(c, d[j], p), p);
        }
    }

    private long Coefficient(int i) => i < _coefficients.Length ? _coefficients[i] : 0;

    private PrimeFieldPolynomial Scale(long scalar)
    {
        var c = new long[_coefficients.Length];
        for (int i = 0; i < c.Length; i++)
            c[i] = ModularArithmetic.MulMod(_coefficients[i], scalar, Prime);
        return new PrimeFieldPolynomial(Trim(c), Prime);
    }

    private static long CommonPrime(PrimeFieldPolynomial a, PrimeFieldPolynomial b)
    {
        if (a is null) throw new ArgumentNullException(nameof(a));
        if (b is null) throw new ArgumentNullException(nameof(b));
        if (a.Prime != b.Prime)
            throw new ArgumentException($"Cannot combine polynomials over GF({a.Prime}) and GF({b.Prime}).");
        return a.Prime;
    }

    private static void ValidatePrime(long prime)
    {
        if (!PrimalityTesting.IsPrime(prime))
            throw new ArgumentException($"{prime} is not prime.", nameof(prime));
    }

    #endregion
}
