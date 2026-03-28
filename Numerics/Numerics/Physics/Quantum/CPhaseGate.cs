using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Controlled-Phase gate CP(θ). Applies a phase of e^(iθ) to the |11⟩ component.
/// CP(θ) = diag(1, 1, 1, e^(iθ))
///
/// First qubit index = control, second = target.
/// Essential for constructing the Quantum Fourier Transform.
/// Special case: CP(π) = CZ.
/// </summary>
public class CPhaseGate : QuantumGate
{
    /// <summary>Phase angle θ in radians.</summary>
    public double Theta { get; }

    public CPhaseGate(double theta)
    {
        Theta = theta;
    }

    public override int QubitCount => 2;

    public override ComplexMatrix GetMatrix()
    {
        var zero = new ComplexNumber(0, 0);
        var one = new ComplexNumber(1, 0);
        var phase = ComplexNumber.FromPolarCoordinates(1, Theta);

        // Gate-local: bit 0 = control (qubitIndices[0]), bit 1 = target (qubitIndices[1])
        // |00⟩→|00⟩, |01⟩→|01⟩, |10⟩→|10⟩, |11⟩→e^(iθ)|11⟩
        return new ComplexMatrix(new ComplexNumber[,]
        {
            { one,  zero, zero, zero  },
            { zero, one,  zero, zero  },
            { zero, zero, one,  zero  },
            { zero, zero, zero, phase }
        });
    }
}
