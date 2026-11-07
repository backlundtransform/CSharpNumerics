using System;
using System.Linq;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Quantum;

namespace NumericsTests;

[TestClass]
public class ShorAlgorithmTests
{
    #region QuantumFourierTransform

    [TestMethod]
    public void Qft_OnBasisState_GivesUniformPhases()
    {
        // QFT|x⟩ = 2^(−n/2) Σ e^(2πi·xy/2ⁿ)|y⟩
        int n = 3, size = 1 << n, x = 5;
        var state = BasisState(n, x);
        var result = new QuantumFourierTransform(n).Apply(state, Enumerable.Range(0, n).ToArray(), n);

        for (int y = 0; y < size; y++)
        {
            var expected = ComplexNumber.FromPolarCoordinates(1 / Math.Sqrt(size), 2 * Math.PI * x * y / size);
            AssertClose(expected, result[y], $"y = {y}");
        }
    }

    [TestMethod]
    public void Qft_CircuitMatchesDenseMatrix()
    {
        var rng = new Random(11);
        foreach (bool inverse in new[] { false, true })
        {
            for (int n = 1; n <= 5; n++)
            {
                // Embed the register in a larger system on scrambled qubits to check the index handling.
                int total = n + 2;
                var qubits = Enumerable.Range(0, total).OrderBy(_ => rng.Next()).Take(n).ToArray();
                var state = RandomState(total, rng);

                var qft = new QuantumFourierTransform(n, inverse);
                var dense = qft.Apply(state, qubits, total);
                var circuit = qft.ApplyAsCircuit(state, qubits, total);

                for (int i = 0; i < state.Length; i++)
                    AssertClose(dense[i], circuit[i], $"n = {n}, inverse = {inverse}, i = {i}");
            }
        }
    }

    [TestMethod]
    public void Qft_InverseUndoesForward()
    {
        var rng = new Random(12);
        int n = 6;
        var qubits = Enumerable.Range(0, n).ToArray();
        var state = RandomState(n, rng);

        var forward = new QuantumFourierTransform(n).ApplyAsCircuit(state, qubits, n);
        var back = new QuantumFourierTransform(n, inverse: true).ApplyAsCircuit(forward, qubits, n);

        for (int i = 0; i < state.Length; i++)
            AssertClose(state[i], back[i], $"i = {i}");
    }

    [TestMethod]
    public void Qft_InvalidArguments()
    {
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => new QuantumFourierTransform(0));
        var qft = new QuantumFourierTransform(2);
        Assert.ThrowsException<ArgumentException>(() => qft.ApplyAsCircuit(BasisState(3, 0), new[] { 0 }, 3));
    }

    #endregion

    #region Classical post-processing

    [TestMethod]
    public void TryRecoverPeriod_FromIdealMeasurements()
    {
        // a = 7 mod 15 has period 4; with t = 8 the peaks are exactly at multiples of 256/4 = 64.
        Assert.IsTrue(ShorAlgorithm.TryRecoverPeriod(64, 8, 7, 15, out long r));
        Assert.AreEqual(4, r);
        Assert.IsTrue(ShorAlgorithm.TryRecoverPeriod(192, 8, 7, 15, out r));
        Assert.AreEqual(4, r);

        // s/r = 2/4 reduces to 1/2: the convergent gives 2, and the multiple search restores 4.
        Assert.IsTrue(ShorAlgorithm.TryRecoverPeriod(128, 8, 7, 15, out r));
        Assert.AreEqual(4, r);

        // m = 0 says nothing about r, and neither does a value far from every peak.
        Assert.IsFalse(ShorAlgorithm.TryRecoverPeriod(0, 8, 7, 15, out _));
        Assert.IsFalse(ShorAlgorithm.TryRecoverPeriod(1, 8, 7, 15, out _));

        // a = 2 mod 21 has period 6, which does not divide 2⁹: the peaks are only near 512·s/6.
        Assert.IsTrue(ShorAlgorithm.TryRecoverPeriod(85, 9, 2, 21, out r));   // 512/6 ≈ 85.3
        Assert.AreEqual(6, r);
        Assert.IsTrue(ShorAlgorithm.TryRecoverPeriod(427, 9, 2, 21, out r));  // 5·512/6 ≈ 426.7
        Assert.AreEqual(6, r);

        Assert.ThrowsException<ArgumentOutOfRangeException>(() => ShorAlgorithm.TryRecoverPeriod(256, 8, 7, 15, out _));
    }

    [TestMethod]
    public void TryFactorFromPeriod_KnownCases()
    {
        // 7⁴ ≡ 1 (mod 15), 7² = 49 ≡ 4: gcd(3, 15) = 3 and gcd(5, 15) = 5.
        Assert.IsTrue(ShorAlgorithm.TryFactorFromPeriod(7, 4, 15, out long f1, out long f2));
        Assert.AreEqual((3L, 5L), (f1, f2));

        // 2⁶ ≡ 1 (mod 21), 2³ = 8: gcd(7, 21) = 7 and gcd(9, 21) = 3.
        Assert.IsTrue(ShorAlgorithm.TryFactorFromPeriod(2, 6, 21, out f1, out f2));
        Assert.AreEqual((3L, 7L), (f1, f2));

        // Odd period: no luck with this base.
        Assert.IsFalse(ShorAlgorithm.TryFactorFromPeriod(4, 3, 21, out _, out _));
        // a^(r/2) ≡ −1: 14 ≡ −1 (mod 15) has period 2, and 14¹ = N − 1 gives only trivial factors.
        Assert.IsFalse(ShorAlgorithm.TryFactorFromPeriod(14, 2, 15, out _, out _));
    }

    #endregion

    #region End-to-end

    [TestMethod]
    public void PeriodFindingState_PeaksAtMultiplesOfQOverR()
    {
        // a = 7 mod 15, period 4, t = 8: all probability on m ∈ {0, 64, 128, 192}, 1/4 each.
        var state = ShorAlgorithm.PeriodFindingState(7, 15);
        int t = ShorAlgorithm.CountingQubitCount(15);
        Assert.AreEqual(8, t);
        Assert.AreEqual(12, state.QubitCount);

        var counting = CountingDistribution(state, t);
        foreach (int m in new[] { 0, 64, 128, 192 })
            Assert.AreEqual(0.25, counting[m], 1e-9, $"m = {m}");
        Assert.AreEqual(1.0, counting.Sum(), 1e-9);
    }

    [TestMethod]
    public void PeriodFindingState_N21_RecoversPeriodWithHighProbability()
    {
        // a = 2 mod 21 has period 6. Sum the probability of every m from which the period is recovered.
        var state = ShorAlgorithm.PeriodFindingState(2, 21);
        int t = ShorAlgorithm.CountingQubitCount(21);
        Assert.AreEqual(14, state.QubitCount);

        var counting = CountingDistribution(state, t);
        double success = 0;
        for (int m = 0; m < counting.Length; m++)
            if (ShorAlgorithm.TryRecoverPeriod(m, t, 2, 21, out long r) && r == 6)
                success += counting[m];

        // s coprime to 6 (s = 1, 5) always works and s = 2, 3, 4 are rescued by the multiple search.
        Assert.IsTrue(success > 0.6, $"success probability {success}");
    }

    [TestMethod]
    public void Factor_15And21()
    {
        // N = 21 simulates 14 qubits, so each seed costs most of a second; N = 15 (12 qubits) is cheap.
        foreach (long n in new long[] { 15, 21 })
        {
            // Some seed must go through the quantum step rather than a lucky gcd.
            var quantum = Enumerable.Range(0, 20)
                .Select(seed => ShorAlgorithm.Factor(n, new Random(seed)))
                .First(r => r.UsedQuantumPeriodFinding);

            Assert.AreEqual((3L, n / 3), (quantum.Factor1, quantum.Factor2), quantum.ToString()); // 15 = 3·5, 21 = 3·7
            Assert.AreEqual(0, quantum.Period % 2);
            Assert.AreEqual(1, CSharpNumerics.Numerics.NumberTheory.ModularArithmetic.Pow(quantum.Base, quantum.Period, n));
        }

        for (int seed = 0; seed < 5; seed++)
        {
            var result = ShorAlgorithm.Factor(15, new Random(seed));
            Assert.AreEqual((3L, 5L), (result.Factor1, result.Factor2), result.ToString());
        }
    }

    [TestMethod]
    public void Factor_ClassicalShortcutsAndValidation()
    {
        var even = ShorAlgorithm.Factor(22, new Random(0));
        Assert.AreEqual((2L, 11L), (even.Factor1, even.Factor2));
        Assert.IsFalse(even.UsedQuantumPeriodFinding);

        var power = ShorAlgorithm.Factor(27, new Random(0));
        Assert.AreEqual((3L, 9L), (power.Factor1, power.Factor2));
        Assert.IsFalse(power.UsedQuantumPeriodFinding);

        Assert.ThrowsException<ArgumentException>(() => ShorAlgorithm.Factor(13, new Random(0)));
        Assert.ThrowsException<ArgumentException>(() => ShorAlgorithm.Factor(3, new Random(0)));
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => ShorAlgorithm.Factor(1003, new Random(0))); // 17·59
        Assert.ThrowsException<ArgumentNullException>(() => ShorAlgorithm.Factor(15, null!));
        Assert.ThrowsException<ArgumentException>(() => ShorAlgorithm.PeriodFindingState(5, 15)); // gcd 5
    }

    #endregion

    #region Helpers

    private static double[] CountingDistribution(QuantumState state, int countingQubits)
    {
        int mask = (1 << countingQubits) - 1;
        var counting = new double[1 << countingQubits];
        for (int i = 0; i < state.Amplitudes.Length; i++)
            counting[i & mask] += state.GetProbability(i);
        return counting;
    }

    private static ComplexVectorN BasisState(int qubits, int index)
    {
        var state = new ComplexVectorN(1 << qubits);
        state[index] = new ComplexNumber(1, 0);
        return state;
    }

    private static ComplexVectorN RandomState(int qubits, Random rng)
    {
        var state = new ComplexVectorN(1 << qubits);
        double norm = 0;
        for (int i = 0; i < state.Length; i++)
        {
            state[i] = new ComplexNumber(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);
            norm += state[i].realPart * state[i].realPart + state[i].imaginaryPart * state[i].imaginaryPart;
        }
        for (int i = 0; i < state.Length; i++)
            state[i] = state[i] * (1 / Math.Sqrt(norm));
        return state;
    }

    private static void AssertClose(ComplexNumber expected, ComplexNumber actual, string message)
    {
        Assert.AreEqual(expected.realPart, actual.realPart, 1e-9, message);
        Assert.AreEqual(expected.imaginaryPart, actual.imaginaryPart, 1e-9, message);
    }

    #endregion
}
