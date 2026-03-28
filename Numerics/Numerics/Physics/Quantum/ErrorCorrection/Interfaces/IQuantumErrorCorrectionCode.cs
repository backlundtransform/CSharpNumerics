using System.Collections.Generic;

namespace CSharpNumerics.Physics.Quantum.ErrorCorrection;

/// <summary>
/// Defines a quantum error-correcting code [[n, k, d]] that encodes k logical qubits
/// into n physical qubits with code distance d.
///
/// Implementations provide the encoding/decoding circuits, syndrome extraction circuits,
/// and a correction map that maps syndrome bit patterns to the corrective gate(s).
/// </summary>
public interface IQuantumErrorCorrectionCode
{
    /// <summary>Number of physical qubits (n).</summary>
    int PhysicalQubits { get; }

    /// <summary>Number of logical qubits (k).</summary>
    int LogicalQubits { get; }

    /// <summary>Code distance (d) — minimum weight of a non-trivial logical operator.</summary>
    int Distance { get; }

    /// <summary>Number of syndrome (ancilla) qubits needed for one error-correction round.</summary>
    int SyndromeQubits { get; }

    /// <summary>
    /// Returns the stabilizer generators as lists of (qubit index, Pauli type) pairs.
    /// Pauli types: 'X', 'Y', 'Z'.
    /// Each generator is a list indicating which Pauli operator acts on which data qubit.
    /// </summary>
    List<List<(int qubit, char pauli)>> GetStabilizers();

    /// <summary>
    /// Returns the correction map: syndrome bit pattern → list of corrective operations.
    /// Each corrective operation is (qubit index, Pauli type) — e.g., (0, 'X') means apply X to qubit 0.
    /// Syndrome 0 (no error) maps to an empty list.
    /// </summary>
    Dictionary<int, List<(int qubit, char pauli)>> GetCorrectionMap();

    /// <summary>
    /// Returns the logical X operator for logical qubit <paramref name="logicalIndex"/>
    /// as a list of (qubit, Pauli) pairs.
    /// </summary>
    List<(int qubit, char pauli)> GetLogicalX(int logicalIndex = 0);

    /// <summary>
    /// Returns the logical Z operator for logical qubit <paramref name="logicalIndex"/>
    /// as a list of (qubit, Pauli) pairs.
    /// </summary>
    List<(int qubit, char pauli)> GetLogicalZ(int logicalIndex = 0);
}
