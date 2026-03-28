using System.Collections.Generic;

namespace CSharpNumerics.Physics.Quantum.ErrorCorrection;

/// <summary>
/// The Steane [[7,1,3]] code — the smallest CSS (Calderbank-Shor-Steane) code 
/// capable of correcting any single-qubit error (X, Z, or Y).
/// Based on the classical [7,4,3] Hamming code and its dual [7,3,4] code.
///
/// Physical qubits: 7 (indices 0–6)
/// Syndrome qubits: 6 (3 for Z-stabilizers detecting X errors,
///                      3 for X-stabilizers detecting Z errors)
///
/// Stabilizer generators (ordered to match syndrome bit assignment):
///   bit 0: Z₃Z₄Z₅Z₆    (Z-type, from Hamming parity check row 0)
///   bit 1: Z₁Z₂Z₅Z₆    (Z-type, from Hamming parity check row 1)
///   bit 2: Z₀Z₂Z₄Z₆    (Z-type, from Hamming parity check row 2)
///   bit 3: X₃X₄X₅X₆    (X-type, same support as bit 0)
///   bit 4: X₁X₂X₅X₆    (X-type, same support as bit 1)
///   bit 5: X₀X₂X₄X₆    (X-type, same support as bit 2)
///
/// The Z-stabilizer syndrome (bits 0–2) identifies X errors:
///   syndrome column for qubit j = binary(j+1), giving the Hamming decoding map
///   {1→q3, 2→q1, 3→q5, 4→q0, 5→q4, 6→q2, 7→q6}.
///
/// The X-stabilizer syndrome (bits 3–5) identifies Z errors with the same map.
///
/// Logical operators:
///   X̄ = X₀X₁X₂X₃X₄X₅X₆   (X on all 7 qubits)
///   Z̄ = Z₀Z₁Z₂Z₃Z₄Z₅Z₆   (Z on all 7 qubits)
/// </summary>
public class SteaneCode7 : IQuantumErrorCorrectionCode
{
    public int PhysicalQubits => 7;
    public int LogicalQubits => 1;
    public int Distance => 3;
    public int SyndromeQubits => 6;

    public List<List<(int qubit, char pauli)>> GetStabilizers()
    {
        return new List<List<(int qubit, char pauli)>>
        {
            // Z-type stabilizers (detect X errors)
            new List<(int, char)> { (3, 'Z'), (4, 'Z'), (5, 'Z'), (6, 'Z') },  // bit 0: Z₃Z₄Z₅Z₆
            new List<(int, char)> { (1, 'Z'), (2, 'Z'), (5, 'Z'), (6, 'Z') },  // bit 1: Z₁Z₂Z₅Z₆
            new List<(int, char)> { (0, 'Z'), (2, 'Z'), (4, 'Z'), (6, 'Z') },  // bit 2: Z₀Z₂Z₄Z₆
            // X-type stabilizers (detect Z errors)
            new List<(int, char)> { (3, 'X'), (4, 'X'), (5, 'X'), (6, 'X') },  // bit 3: X₃X₄X₅X₆
            new List<(int, char)> { (1, 'X'), (2, 'X'), (5, 'X'), (6, 'X') },  // bit 4: X₁X₂X₅X₆
            new List<(int, char)> { (0, 'X'), (2, 'X'), (4, 'X'), (6, 'X') },  // bit 5: X₀X₂X₄X₆
        };
    }

    public Dictionary<int, List<(int qubit, char pauli)>> GetCorrectionMap()
    {
        // The 6-bit syndrome splits into:
        //   Z-part = bits 0–2: identifies which qubit has an X error
        //   X-part = bits 3–5: identifies which qubit has a Z error
        //
        // Both parts use the same Hamming decoding map (self-dual CSS structure):
        //   syndrome value → qubit index
        //   0→none, 1→3, 2→1, 3→5, 4→0, 5→4, 6→2, 7→6
        int[] syndToQubit = { -1, 3, 1, 5, 0, 4, 2, 6 };

        var map = new Dictionary<int, List<(int qubit, char pauli)>>();

        for (int s = 0; s < 64; s++)
        {
            int zPart = s & 7;          // bits 0-2: X-error identification
            int xPart = (s >> 3) & 7;   // bits 3-5: Z-error identification

            var corrections = new List<(int qubit, char pauli)>();

            if (zPart != 0)
                corrections.Add((syndToQubit[zPart], 'X'));

            if (xPart != 0)
                corrections.Add((syndToQubit[xPart], 'Z'));

            map[s] = corrections;
        }

        return map;
    }

    public List<(int qubit, char pauli)> GetLogicalX(int logicalIndex = 0)
    {
        // X̄ = X on all 7 qubits
        var ops = new List<(int, char)>();
        for (int i = 0; i < 7; i++)
            ops.Add((i, 'X'));
        return ops;
    }

    public List<(int qubit, char pauli)> GetLogicalZ(int logicalIndex = 0)
    {
        // Z̄ = Z on all 7 qubits
        var ops = new List<(int, char)>();
        for (int i = 0; i < 7; i++)
            ops.Add((i, 'Z'));
        return ops;
    }
}
