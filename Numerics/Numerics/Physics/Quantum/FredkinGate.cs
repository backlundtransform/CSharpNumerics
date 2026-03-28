using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Fredkin gate (CSWAP) — three-qubit gate that swaps the two target qubits
/// only when the control qubit is |1⟩.
///
/// qubitIndices[0] = control, qubitIndices[1] = target₁, qubitIndices[2] = target₂.
/// Gate-local index: bit 0 = control, bit 1 = target₁, bit 2 = target₂.
/// </summary>
public class FredkinGate : QuantumGate
{
    public override int QubitCount => 3;

    public override ComplexMatrix GetMatrix()
    {
        var zero = new ComplexNumber(0, 0);
        var one = new ComplexNumber(1, 0);

        // 8×8 identity except swap target₁ and target₂ when control=1:
        // Gate-local indices: idx = control + 2*target₁ + 4*target₂
        //   |c=1,t1=1,t2=0⟩ = 3  ↔  |c=1,t1=0,t2=1⟩ = 5
        var m = new ComplexNumber[8, 8];
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 8; j++)
                m[i, j] = zero;

        // Identity for all states except 3 and 5
        for (int i = 0; i < 8; i++)
        {
            if (i != 3 && i != 5)
                m[i, i] = one;
        }

        // Swap |3⟩ ↔ |5⟩ (swap targets when control = 1)
        m[3, 5] = one;
        m[5, 3] = one;

        return new ComplexMatrix(m);
    }
}
