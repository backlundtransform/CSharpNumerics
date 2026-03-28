using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// General single-qubit phase gate P(θ).
/// Applies a phase of e^(iθ) to the |1⟩ component.
/// P(θ) = [[1, 0], [0, e^(iθ)]]
///
/// Special cases: P(π/2) = S, P(π/4) = T, P(π) = Z.
/// </summary>
public class PhaseGate : QuantumGate
{
    /// <summary>Phase angle θ in radians.</summary>
    public double Theta { get; }

    public PhaseGate(double theta)
    {
        Theta = theta;
    }

    public override int QubitCount => 1;

    public override ComplexMatrix GetMatrix()
    {
        return new ComplexMatrix(new ComplexNumber[,]
        {
            { new ComplexNumber(1, 0), new ComplexNumber(0, 0) },
            { new ComplexNumber(0, 0), ComplexNumber.FromPolarCoordinates(1, Theta) }
        });
    }
}
