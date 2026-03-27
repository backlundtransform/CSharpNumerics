using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Represents a quantum noise channel described by a set of Kraus operators {E_k}
/// satisfying the completeness relation Σ E_k† E_k = I.
/// The channel transforms a state |ψ⟩ via Monte Carlo trajectory selection:
/// with probability p_k = ⟨ψ|E_k† E_k|ψ⟩, the state becomes E_k|ψ⟩ / √p_k.
/// </summary>
public interface INoiseChannel
{
    /// <summary>Number of qubits this channel acts on.</summary>
    int QubitCount { get; }

    /// <summary>
    /// Returns the Kraus operators for this channel.
    /// The set {E_k} must satisfy Σ E_k† E_k = I.
    /// </summary>
    ComplexMatrix[] GetKrausOperators();
}
