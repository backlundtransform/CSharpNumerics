using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// S gate (Phase gate) — applies a π/2 phase shift to |1⟩.
/// S = [[1, 0], [0, i]]
/// S² = Z
/// </summary>
public class SGate : QuantumGate
{
    public override int QubitCount => 1;

    public override ComplexMatrix GetMatrix()
    {
        return new ComplexMatrix(new ComplexNumber[,]
        {
            { new ComplexNumber(1, 0), new ComplexNumber(0, 0) },
            { new ComplexNumber(0, 0), new ComplexNumber(0, 1) }
        });
    }
}
