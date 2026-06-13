using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Quantum;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Quantum;

/// <summary>
/// Quantum circuit simulator with noise. After each gate application, noise
/// channels are applied to the affected qubits using Monte Carlo trajectory
/// selection (stochastic Kraus operator sampling).
/// </summary>
public class NoisyQuantumSimulator
{
    private readonly List<INoiseChannel> _noiseChannels = new List<INoiseChannel>();
    private readonly Random _random;

    /// <summary>
    /// Creates a noisy simulator with the given random number generator.
    /// </summary>
    public NoisyQuantumSimulator(Random random)
    {
        _random = random ?? throw new ArgumentNullException(nameof(random));
    }

    /// <summary>
    /// Adds a noise channel that will be applied to every qubit after each gate.
    /// </summary>
    public NoisyQuantumSimulator WithNoise(INoiseChannel channel)
    {
        if (channel == null) throw new ArgumentNullException(nameof(channel));
        _noiseChannels.Add(channel);
        return this;
    }

    /// <summary>
    /// Runs the circuit with noise applied after each gate instruction.
    /// For each gate, every noise channel is applied independently to each
    /// qubit the gate acted on.
    /// </summary>
    public QuantumState Run(QuantumCircuit circuit)
    {
        int stateSize = 1 << circuit.QubitCount;

        // Initialise to |0…0⟩
        var amplitudes = new ComplexVectorN(stateSize);
        amplitudes[0] = new ComplexNumber(1, 0);

        // Apply each gate + noise
        foreach (var instruction in circuit.Instructions)
        {
            // Apply the ideal gate
            amplitudes = instruction.Gate.Apply(
                amplitudes,
                instruction.QubitIndices.ToArray(),
                circuit.QubitCount);

            // Apply noise channels to each affected qubit
            foreach (var channel in _noiseChannels)
            {
                if (channel.QubitCount == 1)
                {
                    foreach (var qubit in instruction.QubitIndices)
                    {
                        amplitudes = ApplySingleQubitNoise(amplitudes, channel, qubit, circuit.QubitCount);
                    }
                }
            }
        }

        return new QuantumState(amplitudes);
    }

    /// <summary>
    /// Applies a single-qubit noise channel to a specific qubit using Monte Carlo
    /// trajectory selection. Selects a Kraus operator E_k with probability
    /// p_k = ⟨ψ|E_k† E_k|ψ⟩ and applies E_k|ψ⟩ / √p_k.
    /// </summary>
    private ComplexVectorN ApplySingleQubitNoise(
        ComplexVectorN amplitudes, INoiseChannel channel, int qubitIndex, int totalQubits)
    {
        var kraus = channel.GetKrausOperators();
        int stateSize = amplitudes.Length;

        // Compute the resulting state for each Kraus operator and its probability
        var candidates = new ComplexVectorN[kraus.Length];
        var probabilities = new double[kraus.Length];

        for (int k = 0; k < kraus.Length; k++)
        {
            candidates[k] = ApplyMatrixToQubit(amplitudes, kraus[k], qubitIndex, totalQubits);

            // p_k = ⟨ψ'|ψ'⟩ where ψ' = E_k|ψ⟩ (before normalisation)
            double norm = 0;
            for (int i = 0; i < stateSize; i++)
            {
                var a = candidates[k][i];
                norm += a.realPart * a.realPart + a.imaginaryPart * a.imaginaryPart;
            }
            probabilities[k] = norm;
        }

        // Monte Carlo selection
        double r = _random.NextDouble();
        double cumulative = 0;
        int selected = kraus.Length - 1;

        for (int k = 0; k < kraus.Length; k++)
        {
            cumulative += probabilities[k];
            if (r < cumulative)
            {
                selected = k;
                break;
            }
        }

        // Normalise the selected trajectory
        double scale = 1.0 / Math.Sqrt(probabilities[selected]);
        var result = new ComplexVectorN(stateSize);
        for (int i = 0; i < stateSize; i++)
        {
            var a = candidates[selected][i];
            result[i] = new ComplexNumber(a.realPart * scale, a.imaginaryPart * scale);
        }

        return result;
    }

    /// <summary>
    /// Applies a 2×2 matrix to a single qubit within a multi-qubit state vector.
    /// </summary>
    private static ComplexVectorN ApplyMatrixToQubit(
        ComplexVectorN amplitudes, ComplexMatrix matrix, int qubitIndex, int totalQubits)
    {
        int stateSize = amplitudes.Length;
        var result = new ComplexVectorN(stateSize);

        int mask = 1 << qubitIndex;

        for (int i = 0; i < stateSize; i++)
        {
            if ((i & mask) != 0) continue; // skip — handled by pair below

            int i0 = i;          // qubit = |0⟩
            int i1 = i | mask;   // qubit = |1⟩

            var a0 = amplitudes[i0];
            var a1 = amplitudes[i1];

            // [result_0]   [m00  m01] [a0]
            // [result_1] = [m10  m11] [a1]
            result[i0] = matrix.values[0, 0] * a0 + matrix.values[0, 1] * a1;
            result[i1] = matrix.values[1, 0] * a0 + matrix.values[1, 1] * a1;
        }

        return result;
    }
}
