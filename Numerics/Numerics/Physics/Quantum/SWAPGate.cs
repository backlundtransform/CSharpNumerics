using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// SWAP gate — two-qubit gate that swaps the states of two qubits.
/// |01⟩ ↔ |10⟩
/// </summary>
public class SWAPGate : QuantumGate
{
    public override int QubitCount => 2;

    public override ComplexMatrix GetMatrix()
    {
        var zero = new ComplexNumber(0, 0);
        var one = new ComplexNumber(1, 0);

        // Gate-local: bit 0 = qubitIndices[0], bit 1 = qubitIndices[1]
        // |00⟩→|00⟩ (0→0), |10⟩↔|01⟩ (1↔2), |11⟩→|11⟩ (3→3)
        return new ComplexMatrix(new ComplexNumber[,]
        {
            { one,  zero, zero, zero },
            { zero, zero, one,  zero },
            { zero, one,  zero, zero },
            { zero, zero, zero, one  }
        });
    }
}
