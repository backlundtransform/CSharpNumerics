using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Hadamard gate — creates equal superposition from a basis state.
/// H = (1/√2) * [[1, 1], [1, -1]]
/// </summary>
public class HadamardGate : QuantumGate
{
    private static readonly double InvSqrt2 = 1.0 / Math.Sqrt(2.0);

    public override int QubitCount => 1;

    public override ComplexMatrix GetMatrix()
    {
        return new ComplexMatrix(new ComplexNumber[,]
        {
            { new ComplexNumber(InvSqrt2, 0), new ComplexNumber(InvSqrt2, 0) },
            { new ComplexNumber(InvSqrt2, 0), new ComplexNumber(-InvSqrt2, 0) }
        });
    }
}
