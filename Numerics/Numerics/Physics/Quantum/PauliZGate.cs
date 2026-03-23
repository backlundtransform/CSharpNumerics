using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Pauli-Z gate — phase flip gate, leaves |0⟩ unchanged and maps |1⟩ → −|1⟩.
/// Z = [[1, 0], [0, -1]]
/// </summary>
public class PauliZGate : QuantumGate
{
    public override int QubitCount => 1;

    public override ComplexMatrix GetMatrix()
    {
        return new ComplexMatrix(new ComplexNumber[,]
        {
            { new ComplexNumber(1, 0), new ComplexNumber(0, 0) },
            { new ComplexNumber(0, 0), new ComplexNumber(-1, 0) }
        });
    }
}
