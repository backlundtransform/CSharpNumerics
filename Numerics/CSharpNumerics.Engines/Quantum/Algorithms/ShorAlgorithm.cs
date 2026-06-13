using CSharpNumerics.Physics.Quantum;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Quantum.Algorithms;

/// <summary>
/// Result of a Shor factorization attempt.
/// </summary>
public class ShorResult
{
    /// <summary>The number that was factored.</summary>
    public long N { get; }

    /// <summary>First non-trivial factor, or -1 if factoring failed.</summary>
    public long Factor1 { get; }

    /// <summary>Second non-trivial factor (N / Factor1), or -1 if factoring failed.</summary>
    public long Factor2 { get; }

    /// <summary>Whether a non-trivial factoring was found.</summary>
    public bool Success => Factor1 > 1 && Factor1 < N;

    /// <summary>The random base a that was used in the successful attempt.</summary>
    public long Base { get; }

    /// <summary>The order r that was found (a^r ≡ 1 mod N).</summary>
    public long Order { get; }

    /// <summary>Number of attempts made before success (or exhaustion).</summary>
    public int Attempts { get; }

    public ShorResult(long n, long factor1, long factor2, long baseA, long order, int attempts)
    {
        N = n;
        Factor1 = factor1;
        Factor2 = factor2;
        Base = baseA;
        Order = order;
        Attempts = attempts;
    }
}

/// <summary>
/// Shor's factoring algorithm for composite integers.
///
/// Combines classical pre/post-processing with quantum order-finding via QPE:
/// 1. Classical: pick random a, check gcd(a, N)
/// 2. Quantum: use QPE with modular multiplication U_a|y⟩ = |ay mod N⟩ to find the order r
/// 3. Classical: extract r via continued fractions, check if r is even, compute gcd(a^(r/2) ± 1, N)
///
/// The quantum circuit uses t counting qubits (precision) and n target qubits
/// where n = ⌈log₂(N)⌉.
/// </summary>
public static class ShorAlgorithm
{
    /// <summary>
    /// Factors the integer N using Shor's algorithm.
    /// </summary>
    /// <param name="N">The integer to factor (must be composite, odd, and not a prime power).</param>
    /// <param name="random">Random number generator.</param>
    /// <param name="maxAttempts">Maximum number of random base attempts (default: 20).</param>
    /// <param name="precisionQubits">Number of counting qubits for QPE (default: 2n where n = ⌈log₂(N)⌉).</param>
    /// <returns>A <see cref="ShorResult"/> with the factors if successful.</returns>
    public static ShorResult Factor(long N, Random random, int maxAttempts = 20, int? precisionQubits = null)
    {
        if (N <= 1)
            throw new ArgumentException("N must be > 1.", nameof(N));
        if (N % 2 == 0)
            return new ShorResult(N, 2, N / 2, 2, 0, 1);
        if (random == null)
            throw new ArgumentNullException(nameof(random));

        // Check for prime power
        for (int k = 2; k <= (int)Math.Log(N, 2) + 1; k++)
        {
            double root = Math.Pow(N, 1.0 / k);
            long r = (long)Math.Round(root);
            if (Power(r, k) == N)
                return new ShorResult(N, r, N / r, r, 0, 1);
        }

        int n = (int)Math.Ceiling(Math.Log(N, 2));
        int t = precisionQubits ?? 2 * n;

        for (int attempt = 1; attempt <= maxAttempts; attempt++)
        {
            long a = 2 + (long)(random.NextDouble() * (N - 3)); // random in [2, N-2]
            long g = GCD(a, N);
            if (g > 1 && g < N)
                return new ShorResult(N, g, N / g, a, 0, attempt);

            // Quantum order-finding
            long order = FindOrder(a, N, t, n);
            if (order <= 0 || order % 2 != 0)
                continue;

            long halfPow = ModPow(a, order / 2, N);
            if (halfPow == N - 1) // a^(r/2) ≡ -1 mod N
                continue;

            long f1 = GCD(halfPow - 1, N);
            long f2 = GCD(halfPow + 1, N);

            if (f1 > 1 && f1 < N)
                return new ShorResult(N, f1, N / f1, a, order, attempt);
            if (f2 > 1 && f2 < N)
                return new ShorResult(N, f2, N / f2, a, order, attempt);
        }

        return new ShorResult(N, -1, -1, 0, 0, maxAttempts);
    }

    /// <summary>
    /// Creates the quantum circuit for order-finding of a modulo N.
    /// This is the quantum core of Shor's algorithm: QPE applied to
    /// the modular multiplication gate U_a|y⟩ = |ay mod N⟩.
    /// </summary>
    /// <param name="a">The base for modular multiplication (coprime to N).</param>
    /// <param name="N">The modulus.</param>
    /// <param name="countingQubits">Number of counting (precision) qubits.</param>
    /// <param name="targetQubits">Number of target qubits (≥ ⌈log₂(N)⌉).</param>
    /// <returns>A quantum circuit. Counting qubits are indices 0..t-1, target qubits are t..t+n-1.
    /// The target register is initialised to |1⟩ (eigenstate of U_a).</returns>
    public static QuantumCircuit CreateOrderFindingCircuit(long a, long N,
        int countingQubits, int targetQubits)
    {
        int totalQubits = countingQubits + targetQubits;
        var circuit = new QuantumCircuit(totalQubits);

        // Initialise target register to |1⟩ (= |00...01⟩): flip qubit t (LSB of target)
        circuit.AddInstruction(new QuantumInstruction(
            new PauliXGate(), new List<int> { countingQubits }));

        // Set up qubit index arrays
        var counting = new int[countingQubits];
        for (int i = 0; i < countingQubits; i++)
            counting[i] = i;

        var target = new int[targetQubits];
        for (int i = 0; i < targetQubits; i++)
            target[i] = countingQubits + i;

        // Build modular multiplication gate and run QPE
        var mulGate = new ModularMultiplyGate(a, N, targetQubits);
        var qpeCircuit = QPE.CreateCircuit(totalQubits, counting, target, mulGate);

        // Copy QPE instructions (skip the initial X gate already added)
        foreach (var instruction in qpeCircuit.Instructions)
            circuit.AddInstruction(instruction);

        return circuit;
    }

    /// <summary>
    /// Extracts the order from QPE measurement results using continued fraction expansion.
    /// </summary>
    /// <param name="measuredValue">The measured value from the counting register.</param>
    /// <param name="precision">Total number of counting qubits (denominator = 2^precision).</param>
    /// <param name="N">The modulus.</param>
    /// <param name="a">The base for modular multiplication.</param>
    /// <returns>The estimated order, or -1 if extraction fails.</returns>
    public static long ExtractOrder(long measuredValue, int precision, long N, long a)
    {
        if (measuredValue == 0)
            return -1;

        long denominator = 1L << precision;
        // Continued fraction expansion of measuredValue / denominator
        var convergents = ContinuedFractionConvergents(measuredValue, denominator);

        foreach (var (_, den) in convergents)
        {
            if (den <= 0 || den >= N)
                continue;

            // Check if den is a valid order: a^den ≡ 1 mod N
            if (ModPow(a, den, N) == 1)
                return den;

            // Try small multiples
            for (int k = 2; k * den < N; k++)
            {
                if (ModPow(a, k * den, N) == 1)
                    return k * den;
            }
        }

        return -1;
    }

    private static long FindOrder(long a, long N, int t, int n)
    {
        var circuit = CreateOrderFindingCircuit(a, N, t, n);
        var sim = new QuantumSimulator();
        var state = sim.Run(circuit);

        // Sample the counting register multiple times to find the order
        int totalQubits = t + n;
        var rng = new Random(42);

        for (int sample = 0; sample < 10; sample++)
        {
            // Measure via sampling
            var counts = state.Sample(1, rng);
            foreach (var kv in counts)
            {
                int fullResult = kv.Key;
                // Extract counting register bits (qubits 0..t-1)
                int countingResult = fullResult & ((1 << t) - 1);

                long order = ExtractOrder(countingResult, t, N, a);
                if (order > 0 && order % 2 == 0)
                    return order;
                if (order > 0)
                    return order;
            }
        }

        return -1;
    }

    /// <summary>
    /// Computes the convergents of the continued fraction expansion of num/den.
    /// Returns a list of (numerator, denominator) pairs.
    /// </summary>
    public static List<(long num, long den)> ContinuedFractionConvergents(long num, long den)
    {
        var convergents = new List<(long, long)>();
        long h0 = 0, h1 = 1;
        long k0 = 1, k1 = 0;

        long n = num, d = den;
        while (d != 0)
        {
            long a = n / d;
            long h2 = a * h1 + h0;
            long k2 = a * k1 + k0;

            convergents.Add((h2, k2));

            h0 = h1; h1 = h2;
            k0 = k1; k1 = k2;

            long temp = n - a * d;
            n = d;
            d = temp;
        }

        return convergents;
    }

    public static long GCD(long a, long b)
    {
        a = Math.Abs(a);
        b = Math.Abs(b);
        while (b != 0)
        {
            long t = b;
            b = a % b;
            a = t;
        }
        return a;
    }

    public static long ModPow(long baseVal, long exp, long mod)
    {
        long result = 1;
        baseVal %= mod;
        while (exp > 0)
        {
            if ((exp & 1) == 1)
                result = result * baseVal % mod;
            exp >>= 1;
            baseVal = baseVal * baseVal % mod;
        }
        return result;
    }

    private static long Power(long baseVal, int exp)
    {
        long result = 1;
        for (int i = 0; i < exp; i++)
        {
            result *= baseVal;
            if (result > 1_000_000_000)
                return -1; // overflow guard
        }
        return result;
    }
}
