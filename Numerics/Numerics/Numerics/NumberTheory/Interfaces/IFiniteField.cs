using System;
using System.Collections.Generic;
using System.Numerics;

namespace CSharpNumerics.Numerics.NumberTheory.Interfaces;

/// <summary>
/// A finite field GF(q) with q = pⁿ elements: the field operations on elements of type
/// <typeparamref name="TElement"/>, so algorithms can be written once for any finite field.
/// </summary>
/// <remarks>
/// Implementations: <see cref="PrimeField"/> (GF(p), elements <see cref="ModInt"/>) and
/// <see cref="ExtensionField"/> (GF(pⁿ), elements <see cref="ExtensionFieldElement"/>).
/// Every method throws <see cref="ArgumentException"/> for an element of another field.
/// </remarks>
/// <typeparam name="TElement">The element type.</typeparam>
public interface IFiniteField<TElement>
{
    /// <summary>The characteristic p, a prime: 1 + 1 + … + 1 (p times) = 0.</summary>
    long Characteristic { get; }

    /// <summary>The degree n over the prime field GF(p).</summary>
    int Degree { get; }

    /// <summary>The number of elements q = pⁿ.</summary>
    BigInteger Order { get; }

    /// <summary>The additive identity.</summary>
    TElement Zero { get; }

    /// <summary>The multiplicative identity.</summary>
    TElement One { get; }

    /// <summary>A generator of the multiplicative group, which is cyclic of order q − 1.</summary>
    /// <exception cref="OverflowException">q − 1 does not fit in a <see cref="long"/>, so it cannot be factored.</exception>
    TElement PrimitiveElement { get; }

    TElement Add(TElement a, TElement b);

    TElement Subtract(TElement a, TElement b);

    TElement Negate(TElement a);

    TElement Multiply(TElement a, TElement b);

    /// <exception cref="DivideByZeroException"><paramref name="b"/> is zero.</exception>
    TElement Divide(TElement a, TElement b);

    /// <exception cref="DivideByZeroException"><paramref name="a"/> is zero.</exception>
    TElement Inverse(TElement a);

    /// <summary><paramref name="a"/> raised to <paramref name="exponent"/>. A negative exponent uses the inverse, and a⁰ = 1, also for a = 0.</summary>
    /// <exception cref="DivideByZeroException"><paramref name="a"/> is zero and <paramref name="exponent"/> is negative.</exception>
    TElement Pow(TElement a, long exponent);

    /// <summary>All q elements, lazily, starting with zero.</summary>
    IEnumerable<TElement> Elements();
}
