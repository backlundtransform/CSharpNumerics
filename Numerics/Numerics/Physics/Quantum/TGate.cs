using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// T gate — applies a π/4 phase shift to |1⟩.
/// T = [[1, 0], [0, e^(iπ/4)]]
/// T² = S
/// </summary>
public class TGate : QuantumGate
{
    public override int QubitCount => 1;

    public override ComplexMatrix GetMatrix()
    {
        double angle = Math.PI / 4.0;
        return new ComplexMatrix(new ComplexNumber[,]
        {
            { new ComplexNumber(1, 0), new ComplexNumber(0, 0) },
            { new ComplexNumber(0, 0), new ComplexNumber(Math.Cos(angle), Math.Sin(angle)) }
        });
    }
}
