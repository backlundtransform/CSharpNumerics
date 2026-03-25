using CSharpNumerics.Engines.Quantum;
using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Static methods for computing quantum state fidelity — the overlap
/// between two quantum states. Fidelity ranges from 0 (orthogonal) to 1 (identical).
/// </summary>
public static class QuantumFidelity
{
    /// <summary>
    /// State fidelity between two pure states: F = |⟨ψ|φ⟩|².
    /// </summary>
    /// <param name="state1">First quantum state.</param>
    /// <param name="state2">Second quantum state.</param>
    /// <returns>Fidelity in [0, 1].</returns>
    public static double Fidelity(QuantumState state1, QuantumState state2)
    {
        if (state1 == null) throw new ArgumentNullException(nameof(state1));
        if (state2 == null) throw new ArgumentNullException(nameof(state2));
        if (state1.Amplitudes.Length != state2.Amplitudes.Length)
            throw new ArgumentException("States must have the same number of qubits.");

        var inner = state1.Amplitudes.HermitianDot(state2.Amplitudes);
        double mag = inner.GetMagnitude();
        return mag * mag;
    }

    /// <summary>
    /// State fidelity between two pure states given as amplitude vectors: F = |⟨ψ|φ⟩|².
    /// </summary>
    public static double Fidelity(ComplexVectorN psi, ComplexVectorN phi)
    {
        if (psi.Length != phi.Length)
            throw new ArgumentException("Vectors must have the same length.");

        var inner = psi.HermitianDot(phi);
        double mag = inner.GetMagnitude();
        return mag * mag;
    }

    /// <summary>
    /// Bloch sphere fidelity for single-qubit states:
    /// F = ½(1 + n̂₁·n̂₂), where n̂ is the Bloch vector.
    /// For pure states this equals |⟨ψ|φ⟩|².
    /// </summary>
    public static double BlochFidelity(BlochVector b1, BlochVector b2)
    {
        if (b1 == null) throw new ArgumentNullException(nameof(b1));
        if (b2 == null) throw new ArgumentNullException(nameof(b2));

        double dot = b1.X * b2.X + b1.Y * b2.Y + b1.Z * b2.Z;
        return 0.5 * (1.0 + dot);
    }
}
