using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Quantum;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Quantum;

/// <summary>
/// Represents the quantum state of a multi-qubit system as a vector of complex amplitudes.
/// For n qubits the state has 2^n amplitudes in the computational basis.
/// </summary>
public class QuantumState
{
    /// <summary>Complex amplitudes in the computational basis.</summary>
    public ComplexVectorN Amplitudes { get; }

    /// <summary>Number of qubits in the system.</summary>
    public int QubitCount { get; }

    public QuantumState(ComplexVectorN amplitudes)
    {
        int n = amplitudes.Length;
        if ((n & (n - 1)) != 0)
            throw new ArgumentException("Amplitude count must be a power of 2.", nameof(amplitudes));

        Amplitudes = amplitudes;
        QubitCount = BitLength(n) - 1;
    }

    /// <summary>
    /// Returns the measurement probability of a specific computational basis state.
    /// </summary>
    public double GetProbability(int basisStateIndex)
    {
        if (basisStateIndex < 0 || basisStateIndex >= Amplitudes.Length)
            throw new ArgumentOutOfRangeException(nameof(basisStateIndex));

        var a = Amplitudes[basisStateIndex];
        return a.realPart * a.realPart + a.imaginaryPart * a.imaginaryPart;
    }

    /// <summary>
    /// Returns measurement probabilities for all computational basis states as a VectorN.
    /// </summary>
    public VectorN GetProbabilities()
    {
        var probs = new double[Amplitudes.Length];
        for (int i = 0; i < Amplitudes.Length; i++)
        {
            var a = Amplitudes[i];
            probs[i] = a.realPart * a.realPart + a.imaginaryPart * a.imaginaryPart;
        }
        return new VectorN(probs);
    }

    /// <summary>
    /// Returns the Bloch vector for a single-qubit state.
    /// Throws if the state has more than one qubit.
    /// </summary>
    public BlochVector GetBlochVector()
    {
        if (QubitCount != 1)
            throw new InvalidOperationException(
                $"Bloch vector is only defined for single-qubit states (this state has {QubitCount} qubits).");

        return BlochVector.FromAmplitudes(Amplitudes[0], Amplitudes[1]);
    }

    /// <summary>
    /// Performs a full measurement in the computational basis, collapsing the state
    /// to a single basis state. Returns the measured basis-state index.
    /// </summary>
    /// <param name="random">Random number generator for sampling.</param>
    /// <returns>Index of the measured computational basis state.</returns>
    public int Measure(Random random)
    {
        if (random == null) throw new ArgumentNullException(nameof(random));

        double r = random.NextDouble();
        double cumulative = 0.0;

        for (int i = 0; i < Amplitudes.Length; i++)
        {
            var a = Amplitudes[i];
            cumulative += a.realPart * a.realPart + a.imaginaryPart * a.imaginaryPart;
            if (r < cumulative)
            {
                Collapse(i);
                return i;
            }
        }

        // Numerical edge case: return last state
        int last = Amplitudes.Length - 1;
        Collapse(last);
        return last;
    }

    /// <summary>
    /// Measures a single qubit, collapsing it to |0⟩ or |1⟩ and renormalising
    /// the remaining state. Returns 0 or 1.
    /// </summary>
    /// <param name="qubitIndex">Index of the qubit to measure (0-based).</param>
    /// <param name="random">Random number generator for sampling.</param>
    /// <returns>0 or 1 — the measurement outcome for the specified qubit.</returns>
    public int MeasureQubit(int qubitIndex, Random random)
    {
        if (qubitIndex < 0 || qubitIndex >= QubitCount)
            throw new ArgumentOutOfRangeException(nameof(qubitIndex));
        if (random == null) throw new ArgumentNullException(nameof(random));

        // Calculate probability of qubit being |0⟩
        double prob0 = 0.0;
        for (int i = 0; i < Amplitudes.Length; i++)
        {
            if ((i & (1 << qubitIndex)) == 0)
            {
                var a = Amplitudes[i];
                prob0 += a.realPart * a.realPart + a.imaginaryPart * a.imaginaryPart;
            }
        }

        int outcome = random.NextDouble() < prob0 ? 0 : 1;

        // Collapse: zero-out amplitudes inconsistent with outcome, renormalise
        var vals = Amplitudes.Values;
        double norm = 0.0;
        for (int i = 0; i < vals.Length; i++)
        {
            bool qubitIs1 = (i & (1 << qubitIndex)) != 0;
            if ((outcome == 0 && qubitIs1) || (outcome == 1 && !qubitIs1))
            {
                vals[i] = new ComplexNumber(0, 0);
            }
            else
            {
                var a = vals[i];
                norm += a.realPart * a.realPart + a.imaginaryPart * a.imaginaryPart;
            }
        }

        // Renormalise surviving amplitudes
        double scale = 1.0 / Math.Sqrt(norm);
        for (int i = 0; i < vals.Length; i++)
        {
            var a = vals[i];
            if (a.realPart != 0 || a.imaginaryPart != 0)
                vals[i] = new ComplexNumber(a.realPart * scale, a.imaginaryPart * scale);
        }

        return outcome;
    }

    /// <summary>
    /// Samples the measurement distribution <paramref name="shots"/> times
    /// without collapsing the state. Returns a dictionary mapping basis-state index to count.
    /// </summary>
    /// <param name="shots">Number of measurement samples.</param>
    /// <param name="random">Random number generator for sampling.</param>
    /// <returns>Dictionary of basis-state index → measurement count.</returns>
    public Dictionary<int, int> Sample(int shots, Random random)
    {
        if (shots <= 0) throw new ArgumentException("Shots must be positive.", nameof(shots));
        if (random == null) throw new ArgumentNullException(nameof(random));

        // Pre-compute cumulative distribution
        var cdf = new double[Amplitudes.Length];
        double cumulative = 0.0;
        for (int i = 0; i < Amplitudes.Length; i++)
        {
            var a = Amplitudes[i];
            cumulative += a.realPart * a.realPart + a.imaginaryPart * a.imaginaryPart;
            cdf[i] = cumulative;
        }

        var counts = new Dictionary<int, int>();
        for (int s = 0; s < shots; s++)
        {
            double r = random.NextDouble();
            int result = 0;
            for (int i = 0; i < cdf.Length; i++)
            {
                if (r < cdf[i])
                {
                    result = i;
                    break;
                }
                if (i == cdf.Length - 1)
                    result = i;
            }

            if (counts.ContainsKey(result))
                counts[result]++;
            else
                counts[result] = 1;
        }

        return counts;
    }

    private void Collapse(int basisStateIndex)
    {
        var vals = Amplitudes.Values;
        for (int i = 0; i < vals.Length; i++)
            vals[i] = i == basisStateIndex
                ? new ComplexNumber(1, 0)
                : new ComplexNumber(0, 0);
    }

    private static int BitLength(int value)
    {
        int bits = 0;
        while (value > 0) { value >>= 1; bits++; }
        return bits;
    }
}
