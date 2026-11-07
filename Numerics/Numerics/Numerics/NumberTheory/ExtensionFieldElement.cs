using System;
using System.Collections.Generic;

namespace CSharpNumerics.Numerics.NumberTheory;

/// <summary>
/// An element of a finite field GF(pⁿ): a polynomial of degree &lt; n over GF(p), taken modulo the
/// defining polynomial of its <see cref="ExtensionField"/>.
/// </summary>
/// <remarks>
/// Elements are created by the field (<see cref="ExtensionField.Element(long)"/>,
/// <see cref="ExtensionField.FromCoefficients"/>, …). Operators that combine elements of different
/// fields throw <see cref="ArgumentException"/>. <c>default(ExtensionFieldElement)</c> belongs to no
/// field and is not a valid value.
/// </remarks>
public readonly struct ExtensionFieldElement : IEquatable<ExtensionFieldElement>
{
    private readonly long[] _coefficients; // length Field.Degree, constant term first, never mutated

    internal ExtensionFieldElement(ExtensionField field, long[] coefficients)
    {
        Field = field;
        _coefficients = coefficients;
    }

    /// <summary>The field the element belongs to.</summary>
    public ExtensionField Field { get; }

    /// <summary>The n coefficients, constant term first.</summary>
    public IReadOnlyList<long> Coefficients => Array.AsReadOnly(_coefficients ?? new long[0]);

    /// <summary>True for the zero element.</summary>
    public bool IsZero
    {
        get
        {
            foreach (long c in Raw)
                if (c != 0) return false;
            return true;
        }
    }

    internal long[] Raw => _coefficients ?? new long[0];

    /// <summary>The index whose base-p digits are the coefficients: the inverse of <see cref="ExtensionField.Element(long)"/>.</summary>
    /// <exception cref="OverflowException">pⁿ does not fit in a <see cref="long"/> and the index is too large.</exception>
    public long ToInt64() => FieldOrThrow.Index(_coefficients);

    /// <summary>The element as a polynomial of degree &lt; n over GF(p).</summary>
    public PrimeFieldPolynomial ToPolynomial() =>
        PrimeFieldPolynomial.FromTrusted((long[])_coefficients.Clone(), FieldOrThrow.Characteristic);

    /// <summary>The element raised to <paramref name="exponent"/>. A negative exponent uses the inverse.</summary>
    /// <exception cref="DivideByZeroException">The element is zero and <paramref name="exponent"/> is negative.</exception>
    public ExtensionFieldElement Pow(long exponent) => FieldOrThrow.Pow(this, exponent);

    /// <summary>The multiplicative inverse.</summary>
    /// <exception cref="DivideByZeroException">The element is zero.</exception>
    public ExtensionFieldElement Inverse() => FieldOrThrow.Inverse(this);

    public static ExtensionFieldElement operator +(ExtensionFieldElement a, ExtensionFieldElement b) => a.FieldOrThrow.Add(a, b);

    public static ExtensionFieldElement operator -(ExtensionFieldElement a, ExtensionFieldElement b) => a.FieldOrThrow.Subtract(a, b);

    public static ExtensionFieldElement operator -(ExtensionFieldElement a) => a.FieldOrThrow.Negate(a);

    public static ExtensionFieldElement operator *(ExtensionFieldElement a, ExtensionFieldElement b) => a.FieldOrThrow.Multiply(a, b);

    /// <exception cref="DivideByZeroException"><paramref name="b"/> is zero.</exception>
    public static ExtensionFieldElement operator /(ExtensionFieldElement a, ExtensionFieldElement b) => a.FieldOrThrow.Divide(a, b);

    public static bool operator ==(ExtensionFieldElement a, ExtensionFieldElement b) => a.Equals(b);

    public static bool operator !=(ExtensionFieldElement a, ExtensionFieldElement b) => !a.Equals(b);

    public bool Equals(ExtensionFieldElement other)
    {
        if (Field is null || other.Field is null) return Field is null && other.Field is null;
        if (!Field.Equals(other.Field)) return false;
        for (int i = 0; i < _coefficients.Length; i++)
            if (_coefficients[i] != other._coefficients[i]) return false;
        return true;
    }

    public override bool Equals(object obj) => obj is ExtensionFieldElement other && Equals(other);

    public override int GetHashCode()
    {
        var hash = new HashCode();
        foreach (long c in Raw) hash.Add(c);
        return hash.ToHashCode();
    }

    /// <summary>Formats as "x^6 + x^4 + x^2 + x + 1 in GF(2^8)".</summary>
    public override string ToString() => Field is null ? "default(ExtensionFieldElement)" : Field.Format(_coefficients);

    private ExtensionField FieldOrThrow =>
        Field ?? throw new ArgumentException("default(ExtensionFieldElement) is not a valid element.");
}
