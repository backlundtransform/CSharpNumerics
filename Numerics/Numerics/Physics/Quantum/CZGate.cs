using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Controlled-Z (CZ) gate — two-qubit gate that applies a phase flip
/// to the |11⟩ state. First qubit index = control, second = target.
/// </summary>
public class CZGate : QuantumGate
{
    public override int QubitCount => 2;

    public override ComplexMatrix GetMatrix()
    {
        var zero = new ComplexNumber(0, 0);
        var one = new ComplexNumber(1, 0);
        var mone = new ComplexNumber(-1, 0);

        // |00⟩ → |00⟩, |01⟩ → |01⟩, |10⟩ → |10⟩, |11⟩ → −|11⟩
        // Gate-local: bit 0 = qubitIndices[0], bit 1 = qubitIndices[1]
        // 0→0, 1→1, 2→2, 3→−3
        return new ComplexMatrix(new ComplexNumber[,]
        {
            { one,  zero, zero, zero },
            { zero, one,  zero, zero },
            { zero, zero, one,  zero },
            { zero, zero, zero, mone }
        });
    }
}
