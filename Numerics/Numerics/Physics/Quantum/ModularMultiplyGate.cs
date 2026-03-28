using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Modular multiplication gate: performs the map |y⟩ → |ay mod N⟩ for y &lt; N,
/// and leaves |y⟩ unchanged for y ≥ N.
///
/// This is the unitary oracle used in Shor's algorithm. The gate acts on n qubits
/// where n = ⌈log₂(N)⌉ (or a specified qubit count ≥ that).
///
/// Requires gcd(a, N) = 1 to ensure the permutation is reversible.
/// The gate matrix is a permutation matrix constructed directly.
/// </summary>
public class ModularMultiplyGate : QuantumGate
{
    private readonly int _qubitCount;

    /// <summary>The multiplier a.</summary>
    public long Multiplier { get; }

    /// <summary>The modulus N.</summary>
    public long Modulus { get; }

    /// <summary>
    /// Creates a modular multiplication gate |y⟩ → |ay mod N⟩.
    /// </summary>
    /// <param name="a">Multiplier (must satisfy gcd(a, N) = 1).</param>
    /// <param name="N">Modulus.</param>
    /// <param name="qubitCount">Number of qubits (must be ≥ ⌈log₂(N)⌉).</param>
    public ModularMultiplyGate(long a, long N, int qubitCount)
    {
        if (N <= 1)
            throw new ArgumentException("Modulus must be > 1.", nameof(N));
        if (a <= 0 || a >= N)
            throw new ArgumentException($"Multiplier must be in range (0, {N}).", nameof(a));
        if (GCD(a, N) != 1)
            throw new ArgumentException($"gcd({a}, {N}) ≠ 1 — multiplier must be coprime to modulus.", nameof(a));

        int minQubits = (int)Math.Ceiling(Math.Log(N, 2));
        if (qubitCount < minQubits)
            throw new ArgumentException(
                $"Need at least {minQubits} qubits to represent values mod {N}.", nameof(qubitCount));

        Multiplier = a;
        Modulus = N;
        _qubitCount = qubitCount;
    }

    public override int QubitCount => _qubitCount;

    public override ComplexMatrix GetMatrix()
    {
        int size = 1 << _qubitCount;
        var zero = new ComplexNumber(0, 0);
        var one = new ComplexNumber(1, 0);

        var m = new ComplexNumber[size, size];
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                m[i, j] = zero;

        for (int y = 0; y < size; y++)
        {
            if (y < Modulus)
            {
                // |y⟩ → |a*y mod N⟩
                int target = (int)((Multiplier * y) % Modulus);
                m[target, y] = one;
            }
            else
            {
                // y ≥ N: identity
                m[y, y] = one;
            }
        }

        return new ComplexMatrix(m);
    }

    private static long GCD(long a, long b)
    {
        while (b != 0)
        {
            long t = b;
            b = a % b;
            a = t;
        }
        return a;
    }
}
