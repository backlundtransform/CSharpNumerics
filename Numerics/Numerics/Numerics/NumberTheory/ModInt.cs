using System;

namespace CSharpNumerics.Numerics.NumberTheory;

/// <summary>
/// An integer residue modulo <see cref="Modulus"/>: an element of ℤ/mℤ.
/// With a prime modulus it is an element of the finite field GF(p).
/// </summary>
/// <remarks>
/// <para>
/// <see cref="Value"/> is always in [0, Modulus). There are no implicit conversions, so the
/// modulus is always explicit. Operations that mix two different moduli throw
/// <see cref="ArgumentException"/>.
/// </para>
/// <para>
/// Division multiplies by the inverse, which exists only when the divisor is coprime to the
/// modulus. <c>default(ModInt)</c> has modulus 0 and is not a valid value.
/// </para>
/// </remarks>
public readonly struct ModInt : IEquatable<ModInt>
{
    /// <summary>The residue, in [0, <see cref="Modulus"/>).</summary>
    public long Value { get; }

    /// <summary>The modulus m &gt; 0.</summary>
    public long Modulus { get; }

    /// <summary>
    /// Creates the residue of <paramref name="value"/> modulo <paramref name="modulus"/>.
    /// Negative values are reduced into [0, modulus).
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="modulus"/> ≤ 0.</exception>
    public ModInt(long value, long modulus)
    {
        Value = ModularArithmetic.Mod(value, modulus);
        Modulus = modulus;
    }

    // Trusted constructor for values already in [0, modulus).
    private ModInt(long value, long modulus, bool _)
    {
        Value = value;
        Modulus = modulus;
    }

    /// <summary>The value raised to <paramref name="exponent"/>. A negative exponent uses the inverse.</summary>
    /// <exception cref="ArgumentException"><paramref name="exponent"/> is negative and the value has no inverse.</exception>
    public ModInt Pow(long exponent) => new ModInt(ModularArithmetic.Pow(Value, exponent, Modulus), Modulus, true);

    /// <summary>The multiplicative inverse.</summary>
    /// <exception cref="ArgumentException">gcd(Value, Modulus) ≠ 1, so no inverse exists.</exception>
    public ModInt Inverse() => new ModInt(ModularArithmetic.Inverse(Value, Modulus), Modulus, true);

    /// <summary>True when the value has a multiplicative inverse, that is gcd(Value, Modulus) = 1.</summary>
    public bool IsUnit => EuclideanAlgorithm.Gcd(Value, Modulus) == 1;

    public static ModInt operator +(ModInt a, ModInt b)
    {
        long m = CommonModulus(a, b);
        return new ModInt(ModularArithmetic.AddMod(a.Value, b.Value, m), m, true);
    }

    public static ModInt operator -(ModInt a, ModInt b)
    {
        long m = CommonModulus(a, b);
        return new ModInt(ModularArithmetic.SubtractMod(a.Value, b.Value, m), m, true);
    }

    public static ModInt operator -(ModInt a) =>
        new ModInt(a.Value == 0 ? 0 : a.Modulus - a.Value, a.Modulus, true);

    public static ModInt operator *(ModInt a, ModInt b)
    {
        long m = CommonModulus(a, b);
        return new ModInt(ModularArithmetic.MulMod(a.Value, b.Value, m), m, true);
    }

    /// <exception cref="ArgumentException">The moduli differ, or <paramref name="b"/> has no inverse.</exception>
    public static ModInt operator /(ModInt a, ModInt b)
    {
        long m = CommonModulus(a, b);
        return new ModInt(ModularArithmetic.MulMod(a.Value, ModularArithmetic.Inverse(b.Value, m), m), m, true);
    }

    public static bool operator ==(ModInt a, ModInt b) => a.Equals(b);

    public static bool operator !=(ModInt a, ModInt b) => !a.Equals(b);

    public bool Equals(ModInt other) => Value == other.Value && Modulus == other.Modulus;

    public override bool Equals(object obj) => obj is ModInt other && Equals(other);

    public override int GetHashCode() => HashCode.Combine(Value, Modulus);

    /// <summary>Formats as "3 (mod 7)".</summary>
    public override string ToString() => $"{Value} (mod {Modulus})";

    private static long CommonModulus(ModInt a, ModInt b)
    {
        if (a.Modulus != b.Modulus)
            throw new ArgumentException($"Cannot combine residues with different moduli ({a.Modulus} and {b.Modulus}).");
        if (a.Modulus <= 0)
            throw new ArgumentException("default(ModInt) is not a valid residue.");
        return a.Modulus;
    }
}
