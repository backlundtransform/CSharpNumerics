using System;
using System.Collections.Generic;
using System.Numerics;

namespace CSharpNumerics.Numerics.NumberTheory;

/// <summary>
/// Solves systems of simultaneous congruences x ≡ rᵢ (mod mᵢ) with the Chinese remainder theorem.
/// </summary>
/// <remarks>
/// <para>
/// The solution is returned as (Remainder, Modulus): every solution is x ≡ Remainder (mod Modulus),
/// with Remainder in [0, Modulus). Modulus is the least common multiple of the moduli, which is
/// their product when they are pairwise coprime. An empty system gives (0, 1).
/// </para>
/// <para>
/// <see cref="Solve(IReadOnlyList{long}, IReadOnlyList{long})"/> requires pairwise coprime moduli.
/// <see cref="TrySolve(IReadOnlyList{long}, IReadOnlyList{long}, out long, out long)"/> handles
/// any moduli and reports when the system has no solution.
/// </para>
/// </remarks>
public static class ChineseRemainderTheorem
{
    /// <summary>
    /// Solves x ≡ remainders[i] (mod moduli[i]) for pairwise coprime moduli.
    /// </summary>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">The lists differ in length, or the moduli are not pairwise coprime.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A modulus is ≤ 0.</exception>
    /// <exception cref="OverflowException">The product of the moduli does not fit in a <see cref="long"/>.</exception>
    public static (long Remainder, long Modulus) Solve(IReadOnlyList<long> remainders, IReadOnlyList<long> moduli)
    {
        ValidateInput(remainders, moduli);

        long r = 0, m = 1;
        for (int i = 0; i < moduli.Count; i++)
        {
            ValidateModulus(moduli[i], i);
            if (EuclideanAlgorithm.Gcd(m, moduli[i]) != 1)
                throw new ArgumentException(
                    $"The moduli are not pairwise coprime (modulus {moduli[i]} at index {i}). Use TrySolve for general moduli.",
                    nameof(moduli));
            Combine(ref r, ref m, remainders[i], moduli[i]);
        }
        return (r, m);
    }

    /// <summary>
    /// Solves x ≡ remainders[i] (mod moduli[i]) for pairwise coprime moduli.
    /// </summary>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">The lists differ in length, or the moduli are not pairwise coprime.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A modulus is ≤ 0.</exception>
    public static (BigInteger Remainder, BigInteger Modulus) Solve(IReadOnlyList<BigInteger> remainders, IReadOnlyList<BigInteger> moduli)
    {
        ValidateInput(remainders, moduli);

        BigInteger r = BigInteger.Zero, m = BigInteger.One;
        for (int i = 0; i < moduli.Count; i++)
        {
            ValidateModulus(moduli[i], i);
            if (!BigInteger.GreatestCommonDivisor(m, moduli[i]).IsOne)
                throw new ArgumentException(
                    $"The moduli are not pairwise coprime (modulus {moduli[i]} at index {i}). Use TrySolve for general moduli.",
                    nameof(moduli));
            Combine(ref r, ref m, remainders[i], moduli[i]);
        }
        return (r, m);
    }

    /// <summary>
    /// Solves x ≡ remainders[i] (mod moduli[i]) for any positive moduli. Returns false when the
    /// congruences are inconsistent, for example x ≡ 1 (mod 4) and x ≡ 2 (mod 6).
    /// </summary>
    /// <param name="remainders">The remainders rᵢ.</param>
    /// <param name="moduli">The moduli mᵢ.</param>
    /// <param name="remainder">The solution, in [0, <paramref name="modulus"/>).</param>
    /// <param name="modulus">The least common multiple of the moduli.</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">The lists differ in length.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A modulus is ≤ 0.</exception>
    /// <exception cref="OverflowException">The least common multiple of the moduli does not fit in a <see cref="long"/>.</exception>
    public static bool TrySolve(IReadOnlyList<long> remainders, IReadOnlyList<long> moduli, out long remainder, out long modulus)
    {
        ValidateInput(remainders, moduli);

        long r = 0, m = 1;
        for (int i = 0; i < moduli.Count; i++)
        {
            ValidateModulus(moduli[i], i);
            if (!Combine(ref r, ref m, remainders[i], moduli[i]))
            {
                remainder = 0;
                modulus = 0;
                return false;
            }
        }
        remainder = r;
        modulus = m;
        return true;
    }

    /// <summary>
    /// Solves x ≡ remainders[i] (mod moduli[i]) for any positive moduli. Returns false when the
    /// congruences are inconsistent.
    /// </summary>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">The lists differ in length.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A modulus is ≤ 0.</exception>
    public static bool TrySolve(IReadOnlyList<BigInteger> remainders, IReadOnlyList<BigInteger> moduli, out BigInteger remainder, out BigInteger modulus)
    {
        ValidateInput(remainders, moduli);

        BigInteger r = BigInteger.Zero, m = BigInteger.One;
        for (int i = 0; i < moduli.Count; i++)
        {
            ValidateModulus(moduli[i], i);
            if (!Combine(ref r, ref m, remainders[i], moduli[i]))
            {
                remainder = BigInteger.Zero;
                modulus = BigInteger.Zero;
                return false;
            }
        }
        remainder = r;
        modulus = m;
        return true;
    }

    /// <summary>
    /// Merges x ≡ r (mod m) with x ≡ r2 (mod m2) into x ≡ r (mod lcm(m, m2)).
    /// Returns false if the two congruences are inconsistent.
    /// </summary>
    private static bool Combine(ref long r, ref long m, long r2, long m2)
    {
        r2 = ModularArithmetic.Mod(r2, m2);
        long g = EuclideanAlgorithm.Gcd(m, m2);

        // r and r2 are both non-negative and below their moduli, so the difference cannot overflow.
        long diff = r2 - r;
        if (diff % g != 0) return false;

        long lcm = checked(m / g * m2);
        long m2g = m2 / g;

        // Solve m·k ≡ diff (mod m2): k ≡ (diff / g)·(m / g)⁻¹ (mod m2 / g).
        long k = m2g == 1
            ? 0
            : ModularArithmetic.MulMod(diff / g, ModularArithmetic.Inverse(m / g, m2g), m2g);

        // k < m2 / g, so m·k + r < lcm and nothing overflows.
        r += m * k;
        m = lcm;
        return true;
    }

    private static bool Combine(ref BigInteger r, ref BigInteger m, BigInteger r2, BigInteger m2)
    {
        r2 = ModularArithmetic.Mod(r2, m2);
        BigInteger g = BigInteger.GreatestCommonDivisor(m, m2);

        BigInteger diff = r2 - r;
        if (!(diff % g).IsZero) return false;

        BigInteger m2g = m2 / g;
        BigInteger k = m2g.IsOne
            ? BigInteger.Zero
            : ModularArithmetic.Mod(diff / g * ModularArithmetic.Inverse(m / g, m2g), m2g);

        r += m * k;
        m = m / g * m2;
        return true;
    }

    private static void ValidateModulus(long modulus, int index)
    {
        if (modulus <= 0)
            throw new ArgumentOutOfRangeException("moduli", modulus, $"Modulus at index {index} must be positive.");
    }

    private static void ValidateModulus(BigInteger modulus, int index)
    {
        if (modulus.Sign <= 0)
            throw new ArgumentOutOfRangeException("moduli", modulus, $"Modulus at index {index} must be positive.");
    }

    private static void ValidateInput<T>(IReadOnlyList<T> remainders, IReadOnlyList<T> moduli)
    {
        if (remainders == null) throw new ArgumentNullException(nameof(remainders));
        if (moduli == null) throw new ArgumentNullException(nameof(moduli));
        if (remainders.Count != moduli.Count)
            throw new ArgumentException("remainders and moduli must have the same length.", nameof(moduli));
    }
}
