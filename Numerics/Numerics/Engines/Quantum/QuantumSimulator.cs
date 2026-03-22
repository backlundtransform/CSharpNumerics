using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Engines.Quantum;

/// <summary>
/// Executes a <see cref="QuantumCircuit"/> by initialising the |0…0⟩ state
/// and applying each instruction sequentially, returning the final <see cref="QuantumState"/>.
/// </summary>
public class QuantumSimulator
{
    /// <summary>
    /// Runs the circuit and returns the resulting quantum state.
    /// </summary>
    public QuantumState Run(QuantumCircuit circuit)
    {
        int stateSize = 1 << circuit.QubitCount;

        // Initialise to |0…0⟩
        var amplitudes = new ComplexVectorN(stateSize);
        amplitudes[0] = new ComplexNumber(1, 0);

        // Apply each gate in sequence
        foreach (var instruction in circuit.Instructions)
        {
            amplitudes = instruction.Gate.Apply(
                amplitudes,
                instruction.QubitIndices.ToArray(),
                circuit.QubitCount);
        }

        return new QuantumState(amplitudes);
    }
}
