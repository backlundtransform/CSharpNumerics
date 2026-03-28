using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Toffoli gate (CCNOT) — three-qubit gate that flips the target qubit
/// only when both control qubits are |1⟩.
///
/// qubitIndices[0] = control₁, qubitIndices[1] = control₂, qubitIndices[2] = target.
/// Gate-local index: bit 0 = control₁, bit 1 = control₂, bit 2 = target.
///
/// Essential for Grover oracles and reversible arithmetic (Shor).
/// </summary>
public class ToffoliGate : QuantumGate
{
    public override int QubitCount => 3;

    public override ComplexMatrix GetMatrix()
    {
        var zero = new ComplexNumber(0, 0);
        var one = new ComplexNumber(1, 0);

        // 8×8 identity except |c1=1,c2=1,t=0⟩ ↔ |c1=1,c2=1,t=1⟩
        // Gate-local indices:
        //   |c1,c2,t⟩ → idx = c1 + 2*c2 + 4*t
        //   |1,1,0⟩ = 3  ↔  |1,1,1⟩ = 7
        var m = new ComplexNumber[8, 8];
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 8; j++)
                m[i, j] = zero;

        // Identity for all states except 3 and 7
        for (int i = 0; i < 8; i++)
        {
            if (i != 3 && i != 7)
                m[i, i] = one;
        }

        // Swap |3⟩ ↔ |7⟩ (flip target when both controls are 1)
        m[3, 7] = one;
        m[7, 3] = one;

        return new ComplexMatrix(m);
    }
}
