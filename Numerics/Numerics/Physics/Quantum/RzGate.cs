using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Rotation about the Z-axis by angle θ.
/// Rz(θ) = [[e^(-iθ/2), 0], [0, e^(iθ/2)]]
/// </summary>
public class RzGate : QuantumGate
{
    public double Theta { get; }

    public RzGate(double theta)
    {
        Theta = theta;
    }

    public override int QubitCount => 1;

    public override ComplexMatrix GetMatrix()
    {
        double half = Theta / 2.0;
        return new ComplexMatrix(new ComplexNumber[,]
        {
            { new ComplexNumber(Math.Cos(-half), Math.Sin(-half)), new ComplexNumber(0, 0) },
            { new ComplexNumber(0, 0), new ComplexNumber(Math.Cos(half), Math.Sin(half)) }
        });
    }
}
