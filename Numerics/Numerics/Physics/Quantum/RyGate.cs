using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Rotation about the Y-axis by angle θ.
/// Ry(θ) = [[cos(θ/2), -sin(θ/2)], [sin(θ/2), cos(θ/2)]]
/// </summary>
public class RyGate : QuantumGate
{
    public double Theta { get; }

    public RyGate(double theta)
    {
        Theta = theta;
    }

    public override int QubitCount => 1;

    public override ComplexMatrix GetMatrix()
    {
        double c = Math.Cos(Theta / 2.0);
        double s = Math.Sin(Theta / 2.0);
        return new ComplexMatrix(new ComplexNumber[,]
        {
            { new ComplexNumber(c, 0),  new ComplexNumber(-s, 0) },
            { new ComplexNumber(s, 0),  new ComplexNumber(c, 0)  }
        });
    }
}
