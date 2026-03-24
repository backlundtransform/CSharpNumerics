using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Represents a single-qubit quantum state on the Bloch sphere.
/// Given amplitudes α (|0⟩) and β (|1⟩):
///   x = 2·Re(α*·β),  y = 2·Im(α*·β),  z = |α|² − |β|²
/// The Bloch vector (x, y, z) lies on or inside the unit sphere.
/// </summary>
public class BlochVector
{
    public double X { get; }
    public double Y { get; }
    public double Z { get; }

    public BlochVector(double x, double y, double z)
    {
        X = x;
        Y = y;
        Z = z;
    }

    /// <summary>
    /// Constructs a Bloch vector from the two amplitudes of a single-qubit state |ψ⟩ = α|0⟩ + β|1⟩.
    /// </summary>
    public static BlochVector FromAmplitudes(ComplexNumber alpha, ComplexNumber beta)
    {
        // α* · β
        var alphaConj = alpha.GetConjugate();
        var product = alphaConj * beta;

        double x = 2.0 * product.realPart;
        double y = 2.0 * product.imaginaryPart;
        double z = alpha.GetMagnitude() * alpha.GetMagnitude()
                 - beta.GetMagnitude() * beta.GetMagnitude();

        return new BlochVector(x, y, z);
    }

    /// <summary>
    /// Polar angle θ ∈ [0, π] from the +Z axis.
    /// |0⟩ → θ=0 (north pole), |1⟩ → θ=π (south pole).
    /// </summary>
    public double Theta => Math.Acos(Math.Max(-1.0, Math.Min(1.0, Z / Radius)));

    /// <summary>
    /// Azimuthal angle φ ∈ (−π, π] in the XY plane.
    /// </summary>
    public double Phi => Math.Atan2(Y, X);

    /// <summary>
    /// Length of the Bloch vector. Equals 1 for pure states, &lt; 1 for mixed states.
    /// </summary>
    public double Radius => Math.Sqrt(X * X + Y * Y + Z * Z);

    /// <summary>
    /// Returns the Bloch vector as a 3D <see cref="Vector"/> for visualization.
    /// </summary>
    public Vector ToVector() => new Vector(X, Y, Z);

    public override string ToString() => $"Bloch({X:F4}, {Y:F4}, {Z:F4})";
}
