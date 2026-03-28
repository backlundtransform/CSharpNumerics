using System.Collections.Generic;

namespace CSharpNumerics.Physics.Quantum.ErrorCorrection;

/// <summary>
/// Shor's 9-qubit code [[9,1,3]].
/// The first quantum error-correcting code — corrects ANY single-qubit error
/// (bit-flip X, phase-flip Z, or combined Y = iXZ).
///
/// Encoding:
///   |0⟩_L = (|000⟩+|111⟩)(|000⟩+|111⟩)(|000⟩+|111⟩) / 2√2
///   |1⟩_L = (|000⟩-|111⟩)(|000⟩-|111⟩)(|000⟩-|111⟩) / 2√2
///
/// Qubits are arranged in 3 blocks of 3:
///   Block 0: qubits 0,1,2
///   Block 1: qubits 3,4,5
///   Block 2: qubits 6,7,8
///
/// 8 stabilizer generators:
///   Bit-flip stabilizers (within blocks):
///     g₁ = Z₀Z₁,  g₂ = Z₁Z₂
///     g₃ = Z₃Z₄,  g₄ = Z₄Z₅
///     g₅ = Z₆Z₇,  g₆ = Z₇Z₈
///   Phase-flip stabilizers (between blocks):
///     g₇ = X₀X₁X₂X₃X₄X₅
///     g₈ = X₃X₄X₅X₆X₇X₈
/// </summary>
public class ShorCode9 : IQuantumErrorCorrectionCode
{
    public int PhysicalQubits => 9;
    public int LogicalQubits => 1;
    public int Distance => 3;
    public int SyndromeQubits => 8;

    public List<List<(int qubit, char pauli)>> GetStabilizers()
    {
        return new List<List<(int qubit, char pauli)>>
        {
            // Bit-flip stabilizers (within-block ZZ pairs)
            new List<(int, char)> { (0, 'Z'), (1, 'Z') },     // g₁: Z₀Z₁
            new List<(int, char)> { (1, 'Z'), (2, 'Z') },     // g₂: Z₁Z₂
            new List<(int, char)> { (3, 'Z'), (4, 'Z') },     // g₃: Z₃Z₄
            new List<(int, char)> { (4, 'Z'), (5, 'Z') },     // g₄: Z₄Z₅
            new List<(int, char)> { (6, 'Z'), (7, 'Z') },     // g₅: Z₆Z₇
            new List<(int, char)> { (7, 'Z'), (8, 'Z') },     // g₆: Z₇Z₈
            // Phase-flip stabilizers (between blocks, 6-qubit X strings)
            new List<(int, char)> { (0, 'X'), (1, 'X'), (2, 'X'), (3, 'X'), (4, 'X'), (5, 'X') },  // g₇
            new List<(int, char)> { (3, 'X'), (4, 'X'), (5, 'X'), (6, 'X'), (7, 'X'), (8, 'X') }   // g₈
        };
    }

    public Dictionary<int, List<(int qubit, char pauli)>> GetCorrectionMap()
    {
        // Syndrome bits 0-5: bit-flip stabilizers g₁–g₆
        // Syndrome bits 6-7: phase-flip stabilizers g₇–g₈
        //
        // Bit-flip detection (within each block, same as 3-qubit bit-flip code):
        //   Block 0 (g₁,g₂):  s0s1 = 01→X₀, 11→X₁, 10→X₂
        //   Block 1 (g₃,g₄):  s2s3 = 01→X₃, 11→X₄, 10→X₅
        //   Block 2 (g₅,g₆):  s4s5 = 01→X₆, 11→X₇, 10→X₈
        //
        // Phase-flip detection (between blocks, uses g₇,g₈):
        //   s6s7 = 01→Z on block 0, 11→Z on block 1, 10→Z on block 2
        //   For the Z correction, applying Z to any qubit in the block suffices (they're equivalent mod stabilizers).

        var map = new Dictionary<int, List<(int qubit, char pauli)>>();

        // Build all 256 syndrome combinations (8 bits)
        for (int syndrome = 0; syndrome < 256; syndrome++)
        {
            var corrections = new List<(int qubit, char pauli)>();

            // Bit-flip corrections for each block
            int bf0 = syndrome & 0b11;         // bits 0,1 → block 0
            int bf1 = (syndrome >> 2) & 0b11;  // bits 2,3 → block 1
            int bf2 = (syndrome >> 4) & 0b11;  // bits 4,5 → block 2

            AddBitFlipCorrection(corrections, bf0, 0);  // block 0: qubits 0,1,2
            AddBitFlipCorrection(corrections, bf1, 3);  // block 1: qubits 3,4,5
            AddBitFlipCorrection(corrections, bf2, 6);  // block 2: qubits 6,7,8

            // Phase-flip corrections
            int pf = (syndrome >> 6) & 0b11;  // bits 6,7
            AddPhaseFlipCorrection(corrections, pf);

            map[syndrome] = corrections;
        }

        return map;
    }

    public List<(int qubit, char pauli)> GetLogicalX(int logicalIndex = 0)
    {
        // Logical X̄ = X₀X₁X₂X₃X₄X₅X₆X₇X₈ (X on all 9 qubits)
        var ops = new List<(int, char)>();
        for (int i = 0; i < 9; i++)
            ops.Add((i, 'X'));
        return ops;
    }

    public List<(int qubit, char pauli)> GetLogicalZ(int logicalIndex = 0)
    {
        // Logical Z̄ = Z₀Z₃Z₆ (one Z per block)
        return new List<(int, char)> { (0, 'Z'), (3, 'Z'), (6, 'Z') };
    }

    /// <summary>
    /// Adds the X correction for a single 3-qubit block based on the 2-bit syndrome.
    /// Same lookup as the 3-qubit bit-flip code:
    ///   00→none, 01→X on first qubit, 10→X on third qubit, 11→X on second qubit.
    /// </summary>
    private static void AddBitFlipCorrection(List<(int qubit, char pauli)> corrections, int syndrome2bit, int blockStart)
    {
        switch (syndrome2bit)
        {
            case 0b01: corrections.Add((blockStart + 0, 'X')); break;
            case 0b10: corrections.Add((blockStart + 2, 'X')); break;
            case 0b11: corrections.Add((blockStart + 1, 'X')); break;
            // 0b00: no bit-flip error in this block
        }
    }

    /// <summary>
    /// Adds the Z correction for the phase-flip between blocks.
    ///   00→none, 01→Z on block 0 (qubit 0), 10→Z on block 2 (qubit 6), 11→Z on block 1 (qubit 3).
    /// </summary>
    private static void AddPhaseFlipCorrection(List<(int qubit, char pauli)> corrections, int syndrome2bit)
    {
        switch (syndrome2bit)
        {
            case 0b01: corrections.Add((0, 'Z')); break;  // block 0
            case 0b10: corrections.Add((6, 'Z')); break;  // block 2
            case 0b11: corrections.Add((3, 'Z')); break;  // block 1
            // 0b00: no phase-flip error
        }
    }
}
