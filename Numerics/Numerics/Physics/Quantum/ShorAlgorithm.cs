using CSharpNumerics.Numerics.NumberTheory;
using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Shor's factoring algorithm on a simulated state vector.
/// </summary>
/// <remarks>
/// <para>
/// To factor N, pick a base a coprime to N and find the period r of aˣ mod N, the smallest r with
/// aʳ ≡ 1 (mod N). If r is even and a^(r/2) ≢ −1, then gcd(a^(r/2) ± 1, N) are nontrivial factors.
/// </para>
/// <para>
/// The period is found by quantum phase estimation. A counting register of t qubits with
/// 2ᵗ ≥ N² is put in uniform superposition. Each counting qubit j controls the modular
/// multiplication |y⟩ → |a^(2ʲ)·y mod N⟩ on a work register that starts in |1⟩. An inverse QFT then
/// concentrates the counting register near multiples of 2ᵗ/r. A measured value m gives m/2ᵗ ≈ s/r,
/// and the continued-fraction convergents of m/2ᵗ recover r (<see cref="TryRecoverPeriod"/>).
/// </para>
/// <para>
/// The state vector has 2^(t + n) amplitudes for n = ⌈log₂ N⌉ work qubits, so the simulation is
/// limited to <see cref="MaxSimulatedQubits"/> qubits in total, which covers N ≤ 64. N = 15 needs
/// 12 qubits and N = 21 needs 14.
/// </para>
/// </remarks>
public static class ShorAlgorithm
{
    /// <summary>Largest total qubit count (counting + work) that <see cref="Factor"/> will simulate.</summary>
    public const int MaxSimulatedQubits = 18;

    /// <summary>Work-register size n: the smallest n with 2ⁿ ≥ N, as required by <see cref="ModularMultiplyGate"/>.</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> &lt; 2.</exception>
    public static int WorkQubitCount(long n)
    {
        if (n < 2) throw new ArgumentOutOfRangeException(nameof(n), n, "N must be at least 2.");
        int bits = 0;
        while (bits < 62 && (1L << bits) < n) bits++;
        return bits;
    }

    /// <summary>Counting-register size t: the smallest t with 2ᵗ ≥ N², which guarantees the period can be recovered.</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> &lt; 2 or N² overflows.</exception>
    public static int CountingQubitCount(long n)
    {
        if (n < 2 || n > 3037000499) throw new ArgumentOutOfRangeException(nameof(n), n, "N must be in [2, 3037000499].");
        long square = n * n;
        int bits = 0;
        while (bits < 62 && (1L << bits) < square) bits++;
        return bits;
    }

    /// <summary>
    /// Builds the period-finding state for base <paramref name="a"/> modulo <paramref name="n"/>:
    /// the state just before measurement. Qubits 0 … t − 1 are the counting register (least
    /// significant first) and qubits t … t + n − 1 the work register.
    /// </summary>
    /// <exception cref="ArgumentException">a is not in [2, N), or gcd(a, N) ≠ 1.</exception>
    /// <exception cref="ArgumentOutOfRangeException">The circuit needs more than <see cref="MaxSimulatedQubits"/> qubits.</exception>
    public static QuantumState PeriodFindingState(long a, long n)
    {
        int t = CountingQubitCount(n);
        int w = WorkQubitCount(n);
        if (t + w > MaxSimulatedQubits)
            throw new ArgumentOutOfRangeException(nameof(n), n,
                $"Factoring {n} needs {t + w} qubits; the simulation is limited to {MaxSimulatedQubits}.");
        if (a < 2 || a >= n)
            throw new ArgumentException($"The base must be in [2, {n}).", nameof(a));
        if (EuclideanAlgorithm.Gcd(a, n) != 1)
            throw new ArgumentException($"gcd({a}, {n}) ≠ 1, so {a} has no period modulo {n}.", nameof(a));

        int total = t + w;
        var workQubits = new int[w];
        for (int i = 0; i < w; i++) workQubits[i] = t + i;

        // |0…0⟩ on the counting register, |1⟩ on the work register.
        var state = new ComplexVectorN(1 << total);
        state[1 << t] = new ComplexNumber(1, 0);

        var hadamard = new HadamardGate();
        for (int j = 0; j < t; j++)
            state = hadamard.Apply(state, new[] { j }, total);

        // Counting qubit j controls multiplication by a^(2ʲ) mod N.
        long power = a;
        var targets = new int[w + 1];
        Array.Copy(workQubits, 0, targets, 1, w);
        for (int j = 0; j < t; j++)
        {
            if (power != 1)
            {
                targets[0] = j;
                state = new ControlledGate(new ModularMultiplyGate(power, n, w)).Apply(state, targets, total);
            }
            power = ModularArithmetic.MulMod(power, power, n);
        }

        var countingQubits = new int[t];
        for (int j = 0; j < t; j++) countingQubits[j] = j;
        state = new QuantumFourierTransform(t, inverse: true).ApplyAsCircuit(state, countingQubits, total);

        return new QuantumState(state);
    }

    /// <summary>
    /// Classical post-processing of a measurement: recovers the period r of aˣ mod N from a
    /// counting-register value m, using the convergents of m/2ᵗ with denominator below N.
    /// Returns false when m carries no usable information, for example m = 0.
    /// </summary>
    /// <remarks>
    /// A convergent with denominator k gives r/gcd(s, r), so the multiples k, 2k, … up to
    /// ⌈log₂ N⌉·k are tested with aʳ ≡ 1 (mod N). The first hit is then reduced to the exact order.
    /// </remarks>
    /// <param name="measurement">The measured counting-register value m, 0 ≤ m &lt; 2ᵗ.</param>
    /// <param name="countingQubits">The counting-register size t.</param>
    /// <param name="a">The base.</param>
    /// <param name="n">The modulus N.</param>
    /// <param name="period">The exact multiplicative order of a modulo N, when found.</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="countingQubits"/> is not in [1, 62], or <paramref name="measurement"/> is outside [0, 2ᵗ).</exception>
    public static bool TryRecoverPeriod(long measurement, int countingQubits, long a, long n, out long period)
    {
        if (countingQubits < 1 || countingQubits > 62)
            throw new ArgumentOutOfRangeException(nameof(countingQubits), countingQubits, "countingQubits must be in [1, 62].");
        long q = 1L << countingQubits;
        if (measurement < 0 || measurement >= q)
            throw new ArgumentOutOfRangeException(nameof(measurement), measurement, $"measurement must be in [0, {q}).");

        period = 0;
        if (measurement == 0 || EuclideanAlgorithm.Gcd(a, n) != 1) return false;

        // Only a few multiples of each denominator are tried, so the search stays far from simply
        // trying every r < N, which would find the period classically and make the measurement irrelevant.
        int maxMultiple = WorkQubitCount(n);

        foreach (var (_, denominator) in ContinuedFractions.Convergents(ContinuedFractions.Expand(measurement, q)))
        {
            if (denominator >= n) break;
            if (denominator == 1) continue; // s/r ≈ 0 or 1 says nothing about r

            // The convergent gives r / gcd(s, r), so the period is a small multiple of the denominator.
            for (long r = denominator; r < n && r <= maxMultiple * denominator; r += denominator)
            {
                if (ModularArithmetic.Pow(a, r, n) == 1)
                {
                    // r is a multiple of the order. Strip prime factors while aʳ stays 1 to reach the order itself.
                    foreach (var (p, _) in IntegerFactorization.Factor(r))
                        while (r % p == 0 && ModularArithmetic.Pow(a, r / p, n) == 1)
                            r /= p;
                    period = r;
                    return true;
                }
            }
        }
        return false;
    }

    /// <summary>
    /// Classical post-processing of a period: if r is even and a^(r/2) ≢ −1 (mod N), returns the
    /// nontrivial factors gcd(a^(r/2) − 1, N) and N divided by it.
    /// </summary>
    /// <returns>False when the period does not lead to a factor, so another base must be tried.</returns>
    public static bool TryFactorFromPeriod(long a, long period, long n, out long factor1, out long factor2)
    {
        factor1 = factor2 = 0;
        if (period <= 0 || (period & 1) == 1 || n < 4) return false;

        long x = ModularArithmetic.Pow(a, period / 2, n);
        if (x == n - 1 || x == 1) return false;

        long f = EuclideanAlgorithm.Gcd(x - 1, n);
        if (f == 1 || f == n) f = EuclideanAlgorithm.Gcd(x + 1, n);
        if (f == 1 || f == n) return false;

        factor1 = Math.Min(f, n / f);
        factor2 = Math.Max(f, n / f);
        return true;
    }

    /// <summary>
    /// Factors a composite N with Shor's algorithm: classical shortcuts first (N even, N a perfect
    /// power, gcd(a, N) &gt; 1), otherwise simulated quantum period finding with random bases.
    /// </summary>
    /// <param name="n">An odd or even composite number.</param>
    /// <param name="random">Source of the random bases and the measurement outcomes.</param>
    /// <param name="maxAttempts">How many bases to try before giving up.</param>
    /// <exception cref="ArgumentNullException"><paramref name="random"/> is null.</exception>
    /// <exception cref="ArgumentException">N is less than 4 or prime.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Period finding for N needs more than <see cref="MaxSimulatedQubits"/> qubits.</exception>
    /// <exception cref="InvalidOperationException">No factor was found within <paramref name="maxAttempts"/> attempts.</exception>
    public static ShorResult Factor(long n, Random random, int maxAttempts = 20)
    {
        if (random == null) throw new ArgumentNullException(nameof(random));
        if (n < 4 || PrimalityTesting.IsPrime(n))
            throw new ArgumentException($"{n} is not composite.", nameof(n));

        if ((n & 1) == 0) return new ShorResult(n, 2, n / 2, false, 0, 0, 0, 0);

        // Shor's method cannot split a prime power pᵏ, but a perfect power is easy to detect classically.
        if (TryPerfectPower(n, out long root)) return new ShorResult(n, root, n / root, false, 0, 0, 0, 0);

        int t = CountingQubitCount(n);
        int w = WorkQubitCount(n);
        if (t + w > MaxSimulatedQubits)
            throw new ArgumentOutOfRangeException(nameof(n), n,
                $"Factoring {n} needs {t + w} qubits; the simulation is limited to {MaxSimulatedQubits}.");

        long mask = (1L << t) - 1;
        for (int attempt = 1; attempt <= maxAttempts; attempt++)
        {
            long a = random.Next(2, (int)n);
            long g = EuclideanAlgorithm.Gcd(a, n);
            if (g > 1) return new ShorResult(n, g, n / g, false, a, 0, 0, attempt);

            var state = PeriodFindingState(a, n);
            long measurement = state.Measure(random) & mask;

            if (TryRecoverPeriod(measurement, t, a, n, out long period) &&
                TryFactorFromPeriod(a, period, n, out long f1, out long f2))
                return new ShorResult(n, f1, f2, true, a, period, measurement, attempt);
        }

        throw new InvalidOperationException($"No factor of {n} found in {maxAttempts} attempts.");
    }

    /// <summary>Finds a nontrivial root b with bᵏ = n for some k ≥ 2, if n is a perfect power.</summary>
    private static bool TryPerfectPower(long n, out long root)
    {
        for (int k = 2; (1L << k) <= n; k++)
        {
            long b = (long)Math.Round(Math.Pow(n, 1.0 / k));
            for (long candidate = Math.Max(2, b - 1); candidate <= b + 1; candidate++)
            {
                long value = 1;
                int e = 0;
                while (e < k && value <= n / candidate)
                {
                    value *= candidate;
                    e++;
                }
                if (e == k && value == n)
                {
                    root = candidate;
                    return true;
                }
            }
        }
        root = 0;
        return false;
    }
}
