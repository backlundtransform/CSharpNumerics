using CSharpNumerics.Numerics.NumberTheory.Interfaces;
using System;
using System.Collections.Generic;
using System.Numerics;

namespace CSharpNumerics.Numerics.NumberTheory;

/// <summary>
/// The prime field GF(p) = ℤ/pℤ, with <see cref="ModInt"/> elements of modulus p.
/// </summary>
/// <remarks>
/// <para>
/// <see cref="ModInt"/> already has the arithmetic operators. This class adds the field view:
/// it checks that elements belong to GF(p), reports division by zero as
/// <see cref="DivideByZeroException"/>, and implements <see cref="IFiniteField{TElement}"/>.
/// </para>
/// <para>
/// Not constant-time: do not use for cryptographic purposes.
/// </para>
/// </remarks>
/// <example>
/// <code>
/// var gf7 = new PrimeField(7);
/// ModInt a = gf7.Element(3);
/// ModInt b = gf7.Inverse(a);     // 5 (mod 7), since 3 · 5 = 15 ≡ 1
/// ModInt g = gf7.PrimitiveElement; // 3 (mod 7)
/// </code>
/// </example>
public sealed class PrimeField : IFiniteField<ModInt>, IEquatable<PrimeField>
{
    private ModInt? _primitiveElement;

    /// <summary>Creates GF(<paramref name="prime"/>).</summary>
    /// <exception cref="ArgumentException"><paramref name="prime"/> is not prime.</exception>
    public PrimeField(long prime)
    {
        if (!PrimalityTesting.IsPrime(prime))
            throw new ArgumentException($"{prime} is not prime.", nameof(prime));
        Prime = prime;
    }

    /// <summary>The prime p.</summary>
    public long Prime { get; }

    public long Characteristic => Prime;

    public int Degree => 1;

    public BigInteger Order => Prime;

    public ModInt Zero => new ModInt(0, Prime);

    public ModInt One => new ModInt(1, Prime);

    /// <summary>The smallest primitive root modulo p.</summary>
    public ModInt PrimitiveElement =>
        _primitiveElement ??= new ModInt(ModularArithmetic.PrimitiveRoot(Prime), Prime);

    /// <summary>The residue of <paramref name="value"/> modulo p. Negative values are reduced into [0, p).</summary>
    public ModInt Element(long value) => new ModInt(value, Prime);

    /// <summary>The polynomial with the given coefficients (constant term first) over this field.</summary>
    public PrimeFieldPolynomial Polynomial(params long[] coefficients) => new PrimeFieldPolynomial(coefficients, Prime);

    public ModInt Add(ModInt a, ModInt b) => Check(a) + Check(b);

    public ModInt Subtract(ModInt a, ModInt b) => Check(a) - Check(b);

    public ModInt Negate(ModInt a) => -Check(a);

    public ModInt Multiply(ModInt a, ModInt b) => Check(a) * Check(b);

    public ModInt Divide(ModInt a, ModInt b) => Check(a) * Inverse(b);

    public ModInt Inverse(ModInt a)
    {
        if (Check(a).Value == 0) throw new DivideByZeroException("Zero has no inverse.");
        return a.Inverse();
    }

    public ModInt Pow(ModInt a, long exponent)
    {
        if (Check(a).Value == 0 && exponent < 0) throw new DivideByZeroException("Zero has no inverse.");
        return a.Pow(exponent);
    }

    /// <summary>The order of <paramref name="a"/> in the multiplicative group: the smallest k ≥ 1 with aᵏ = 1. It divides p − 1.</summary>
    /// <exception cref="ArgumentException"><paramref name="a"/> is zero.</exception>
    public long MultiplicativeOrder(ModInt a)
    {
        if (Check(a).Value == 0) throw new ArgumentException("Zero has no multiplicative order.", nameof(a));
        return ModularArithmetic.MultiplicativeOrder(a.Value, Prime);
    }

    public IEnumerable<ModInt> Elements()
    {
        for (long v = 0; v < Prime; v++)
            yield return new ModInt(v, Prime);
    }

    public bool Equals(PrimeField other) => other is not null && Prime == other.Prime;

    public override bool Equals(object obj) => obj is PrimeField other && Equals(other);

    public override int GetHashCode() => Prime.GetHashCode();

    /// <summary>Formats as "GF(7)".</summary>
    public override string ToString() => $"GF({Prime})";

    private ModInt Check(ModInt a)
    {
        if (a.Modulus != Prime)
            throw new ArgumentException($"{a} is not an element of GF({Prime}).");
        return a;
    }
}
