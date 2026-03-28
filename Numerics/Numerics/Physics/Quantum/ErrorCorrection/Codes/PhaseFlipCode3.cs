using System.Collections.Generic;

namespace CSharpNumerics.Physics.Quantum.ErrorCorrection;

/// <summary>
/// Three-qubit phase-flip code [[3,1,1]].
/// Encodes |ψ⟩ = α|0⟩ + β|1⟩  into  α|+++⟩ + β|---⟩.
/// Corrects any single-qubit Z (phase-flip) error.
///
/// The encoding applies H to each qubit after the bit-flip encoding,
/// converting bit-flip protection into phase-flip protection.
///
/// Stabilizers: X₀X₁ and X₁X₂.
/// Syndrome 00 → no error, 10 → Z₀, 11 → Z₁, 01 → Z₂.
/// </summary>
public class PhaseFlipCode3 : IQuantumErrorCorrectionCode
{
    public int PhysicalQubits => 3;
    public int LogicalQubits => 1;
    public int Distance => 1;
    public int SyndromeQubits => 2;

    public List<List<(int qubit, char pauli)>> GetStabilizers()
    {
        return new List<List<(int qubit, char pauli)>>
        {
            new List<(int, char)> { (0, 'X'), (1, 'X') },
            new List<(int, char)> { (1, 'X'), (2, 'X') }
        };
    }

    public Dictionary<int, List<(int qubit, char pauli)>> GetCorrectionMap()
    {
        // Syndrome is a 2-bit integer: bit 0 = stabilizer 0 (X₀X₁), bit 1 = stabilizer 1 (X₁X₂)
        // Z₀ ⇒ X₀X₁=-1, X₁X₂=+1 ⇒ s0=1, s1=0 ⇒ 0b01
        // Z₁ ⇒ X₀X₁=-1, X₁X₂=-1 ⇒ s0=1, s1=1 ⇒ 0b11
        // Z₂ ⇒ X₀X₁=+1, X₁X₂=-1 ⇒ s0=0, s1=1 ⇒ 0b10
        return new Dictionary<int, List<(int qubit, char pauli)>>
        {
            { 0b00, new List<(int, char)>() },                     // no error
            { 0b01, new List<(int, char)> { (0, 'Z') } },          // X₀X₁ triggered → phase error on q0
            { 0b10, new List<(int, char)> { (2, 'Z') } },          // X₁X₂ triggered → phase error on q2
            { 0b11, new List<(int, char)> { (1, 'Z') } }           // both triggered → phase error on q1
        };
    }

    public List<(int qubit, char pauli)> GetLogicalX(int logicalIndex = 0)
    {
        // Logical X = X₀ (or X₁ or X₂ — equivalent up to stabilizers)
        return new List<(int, char)> { (0, 'X') };
    }

    public List<(int qubit, char pauli)> GetLogicalZ(int logicalIndex = 0)
    {
        // Logical Z = Z₀Z₁Z₂ (phase-flip all three)
        return new List<(int, char)> { (0, 'Z'), (1, 'Z'), (2, 'Z') };
    }
}
