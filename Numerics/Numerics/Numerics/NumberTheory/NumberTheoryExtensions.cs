namespace CSharpNumerics.Numerics.NumberTheory;

/// <summary>
/// Extension methods for the common number-theory operations on <see cref="long"/>.
/// The full API is on <see cref="EuclideanAlgorithm"/> and <see cref="ModularArithmetic"/>.
/// </summary>
/// <remarks>
/// C# does not apply numeric conversions to the receiver of an extension method, so an
/// <see cref="int"/> receiver needs a cast or an L suffix: <c>84L.Gcd(36)</c>.
/// </remarks>
public static class NumberTheoryExtensions
{
    /// <inheritdoc cref="EuclideanAlgorithm.Gcd(long, long)"/>
    public static long Gcd(this long a, long b) => EuclideanAlgorithm.Gcd(a, b);

    /// <inheritdoc cref="EuclideanAlgorithm.Lcm(long, long)"/>
    public static long Lcm(this long a, long b) => EuclideanAlgorithm.Lcm(a, b);

    /// <inheritdoc cref="ModularArithmetic.Pow(long, long, long)"/>
    public static long ModPow(this long b, long exponent, long m) => ModularArithmetic.Pow(b, exponent, m);

    /// <inheritdoc cref="ModularArithmetic.Inverse(long, long)"/>
    public static long ModInverse(this long a, long m) => ModularArithmetic.Inverse(a, m);

    /// <inheritdoc cref="PrimalityTesting.IsPrime(long)"/>
    public static bool IsPrime(this long n) => PrimalityTesting.IsPrime(n);
}
