using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Controlled-NOT (CNOT) gate — two-qubit gate that flips the target qubit
/// when the control qubit is |1⟩.
/// First qubit index = control, second = target.
/// </summary>
public class CNOTGate : QuantumGate
{
    public override int QubitCount => 2;

    public override ComplexMatrix GetMatrix()
    {
        var zero = new ComplexNumber(0, 0);
        var one = new ComplexNumber(1, 0);

        // Gate-local index: bit 0 = control (qubitIndices[0]), bit 1 = target (qubitIndices[1])
        // |c=0,t=0⟩ → |c=0,t=0⟩   (0 → 0)
        // |c=1,t=0⟩ → |c=1,t=1⟩   (1 → 3)
        // |c=0,t=1⟩ → |c=0,t=1⟩   (2 → 2)
        // |c=1,t=1⟩ → |c=1,t=0⟩   (3 → 1)
        return new ComplexMatrix(new ComplexNumber[,]
        {
            { one,  zero, zero, zero },
            { zero, zero, zero, one  },
            { zero, zero, one,  zero },
            { zero, one,  zero, zero }
        });
    }
}
