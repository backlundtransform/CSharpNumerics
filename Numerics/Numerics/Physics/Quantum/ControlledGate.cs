using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Wraps any <see cref="QuantumGate"/> with a control qubit, producing an
/// (n+1)-qubit controlled gate. When the control qubit is |1⟩ the inner gate
/// is applied; when |0⟩ the target qubits are left unchanged.
///
/// Qubit ordering: qubitIndices[0] = control, qubitIndices[1..n] = inner gate targets.
/// Gate-local bit 0 = control bit.
/// </summary>
public class ControlledGate : QuantumGate
{
    private readonly QuantumGate _innerGate;

    /// <summary>Creates a controlled version of the given gate.</summary>
    /// <param name="innerGate">The gate to control.</param>
    public ControlledGate(QuantumGate innerGate)
    {
        _innerGate = innerGate ?? throw new ArgumentNullException(nameof(innerGate));
    }

    /// <summary>Inner gate plus one control qubit.</summary>
    public override int QubitCount => _innerGate.QubitCount + 1;

    public override ComplexMatrix GetMatrix()
    {
        int innerSize = 1 << _innerGate.QubitCount;
        int totalSize = innerSize * 2; // 2^(n+1)
        var innerMatrix = _innerGate.GetMatrix();

        var zero = new ComplexNumber(0, 0);
        var one = new ComplexNumber(1, 0);
        var m = new ComplexNumber[totalSize, totalSize];

        for (int i = 0; i < totalSize; i++)
            for (int j = 0; j < totalSize; j++)
                m[i, j] = zero;

        // Gate-local bit 0 = control.
        // For each pair (row, col):
        //   control bit of row = row & 1, inner bits of row = row >> 1
        //   control bit of col = col & 1, inner bits of col = col >> 1
        // If control bits differ → 0
        // If control = 0 → identity on inner bits
        // If control = 1 → inner gate matrix on inner bits

        for (int row = 0; row < totalSize; row++)
        {
            int rowCtrl = row & 1;
            int rowInner = row >> 1;

            for (int col = 0; col < totalSize; col++)
            {
                int colCtrl = col & 1;
                int colInner = col >> 1;

                if (rowCtrl != colCtrl)
                    continue; // control bits must match

                if (rowCtrl == 0)
                {
                    // Control = |0⟩ → identity
                    m[row, col] = (rowInner == colInner) ? one : zero;
                }
                else
                {
                    // Control = |1⟩ → apply inner gate
                    m[row, col] = innerMatrix.values[rowInner, colInner];
                }
            }
        }

        return new ComplexMatrix(m);
    }
}
