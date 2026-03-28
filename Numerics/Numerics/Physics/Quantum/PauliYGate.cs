using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Pauli-Y gate.
/// Y = [[0, -i], [i, 0]]
///
/// Y|0⟩ = i|1⟩, Y|1⟩ = -i|0⟩.
/// Completes the Pauli group {I, X, Y, Z}.
/// </summary>
public class PauliYGate : QuantumGate
{
    public override int QubitCount => 1;

    public override ComplexMatrix GetMatrix()
    {
        return new ComplexMatrix(new ComplexNumber[,]
        {
            { new ComplexNumber(0, 0), new ComplexNumber(0, -1) },
            { new ComplexNumber(0, 1), new ComplexNumber(0, 0) }
        });
    }
}
