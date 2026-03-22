using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Pauli-X gate — quantum NOT gate, flips |0⟩ ↔ |1⟩.
/// X = [[0, 1], [1, 0]]
/// </summary>
public class PauliXGate : QuantumGate
{
    public override int QubitCount => 1;

    public override ComplexMatrix GetMatrix()
    {
        return new ComplexMatrix(new ComplexNumber[,]
        {
            { new ComplexNumber(0, 0), new ComplexNumber(1, 0) },
            { new ComplexNumber(1, 0), new ComplexNumber(0, 0) }
        });
    }
}
