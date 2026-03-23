using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Rotation about the X-axis by angle θ.
/// Rx(θ) = [[cos(θ/2), -i·sin(θ/2)], [-i·sin(θ/2), cos(θ/2)]]
/// </summary>
public class RxGate : QuantumGate
{
    public double Theta { get; }

    public RxGate(double theta)
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
            { new ComplexNumber(c, 0),  new ComplexNumber(0, -s) },
            { new ComplexNumber(0, -s), new ComplexNumber(c, 0)  }
        });
    }
}
