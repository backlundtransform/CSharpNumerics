using System.Collections.Generic;

namespace CSharpNumerics.Physics.Quantum.ErrorCorrection;

/// <summary>
/// Three-qubit bit-flip code [[3,1,1]].
/// Encodes |ψ⟩ = α|0⟩ + β|1⟩  into  α|000⟩ + β|111⟩.
/// Corrects any single-qubit X (bit-flip) error.
///
/// Stabilizers: Z₀Z₁ and Z₁Z₂.
/// Syndrome 00 → no error, 10 → X₀, 11 → X₁, 01 → X₂.
/// </summary>
public class BitFlipCode3 : IQuantumErrorCorrectionCode
{
    public int PhysicalQubits => 3;
    public int LogicalQubits => 1;
    public int Distance => 1;
    public int SyndromeQubits => 2;

    public List<List<(int qubit, char pauli)>> GetStabilizers()
    {
        return new List<List<(int qubit, char pauli)>>
        {
            new List<(int, char)> { (0, 'Z'), (1, 'Z') },
            new List<(int, char)> { (1, 'Z'), (2, 'Z') }
        };
    }

    public Dictionary<int, List<(int qubit, char pauli)>> GetCorrectionMap()
    {
        // Syndrome is a 2-bit integer: bit 0 = stabilizer 0 (Z₀Z₁), bit 1 = stabilizer 1 (Z₁Z₂)
        // X₀ ⇒ Z₀Z₁=-1, Z₁Z₂=+1 ⇒ s0=1, s1=0 ⇒ 0b01
        // X₁ ⇒ Z₀Z₁=-1, Z₁Z₂=-1 ⇒ s0=1, s1=1 ⇒ 0b11
        // X₂ ⇒ Z₀Z₁=+1, Z₁Z₂=-1 ⇒ s0=0, s1=1 ⇒ 0b10
        return new Dictionary<int, List<(int qubit, char pauli)>>
        {
            { 0b00, new List<(int, char)>() },                     // no error
            { 0b01, new List<(int, char)> { (0, 'X') } },          // Z₀Z₁ triggered → error on q0
            { 0b10, new List<(int, char)> { (2, 'X') } },          // Z₁Z₂ triggered → error on q2
            { 0b11, new List<(int, char)> { (1, 'X') } }           // both triggered → error on q1
        };
    }

    public List<(int qubit, char pauli)> GetLogicalX(int logicalIndex = 0)
    {
        // Logical X = X₀X₁X₂ (bit-flip all three)
        return new List<(int, char)> { (0, 'X'), (1, 'X'), (2, 'X') };
    }

    public List<(int qubit, char pauli)> GetLogicalZ(int logicalIndex = 0)
    {
        // Logical Z = Z₀ (or Z₁ or Z₂ — all equivalent up to stabilizers)
        return new List<(int, char)> { (0, 'Z') };
    }
}
