using CSharpNumerics.Physics.Quantum.ErrorCorrection;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Quantum.ErrorCorrection;

/// <summary>
/// Classical syndrome decoder: given a syndrome bit pattern, looks up the
/// corrective operation(s) from a <see cref="IQuantumErrorCorrectionCode"/>'s
/// correction map. Returns the list of (qubit, Pauli) operations to apply.
/// </summary>
public class SyndromeDecoder
{
    private readonly Dictionary<int, List<(int qubit, char pauli)>> _correctionMap;

    /// <summary>
    /// Creates a decoder from a QEC code's correction map.
    /// </summary>
    public SyndromeDecoder(IQuantumErrorCorrectionCode code)
    {
        _correctionMap = code.GetCorrectionMap();
    }

    /// <summary>
    /// Creates a decoder from an explicit correction map.
    /// </summary>
    public SyndromeDecoder(Dictionary<int, List<(int qubit, char pauli)>> correctionMap)
    {
        _correctionMap = correctionMap;
    }

    /// <summary>
    /// Decodes a syndrome value (integer formed from syndrome bits) and returns
    /// the corrective operations. Returns an empty list if the syndrome is not
    /// recognised (i.e. the error is beyond the code's correction capability).
    /// </summary>
    /// <param name="syndrome">Syndrome value as an integer (bit 0 = first stabilizer, etc.).</param>
    public List<(int qubit, char pauli)> Decode(int syndrome)
    {
        if (_correctionMap.TryGetValue(syndrome, out var corrections))
            return corrections;

        return new List<(int qubit, char pauli)>();
    }

    /// <summary>
    /// Decodes a syndrome given as an array of individual measurement bits.
    /// bit[0] is the least significant bit (stabilizer 0).
    /// </summary>
    /// <param name="syndromeBits">Array of 0/1 measurement outcomes for each ancilla qubit.</param>
    public List<(int qubit, char pauli)> Decode(int[] syndromeBits)
    {
        int syndrome = 0;
        for (int i = 0; i < syndromeBits.Length; i++)
        {
            if (syndromeBits[i] == 1)
                syndrome |= (1 << i);
        }
        return Decode(syndrome);
    }
}
