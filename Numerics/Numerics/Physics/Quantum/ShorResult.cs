namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// The outcome of <see cref="ShorAlgorithm.Factor"/>: a nontrivial factorisation N = Factor1 · Factor2
/// and how it was found.
/// </summary>
public class ShorResult
{
    /// <summary>The number that was factored.</summary>
    public long N { get; }

    /// <summary>The smaller factor, 1 &lt; Factor1 ≤ Factor2.</summary>
    public long Factor1 { get; }

    /// <summary>The larger factor, N / Factor1.</summary>
    public long Factor2 { get; }

    /// <summary>
    /// True when the factors came from quantum period finding. False when a classical shortcut
    /// found them first: N even, N a perfect power, or a random base sharing a factor with N.
    /// </summary>
    public bool UsedQuantumPeriodFinding { get; }

    /// <summary>The base a of the successful attempt, or 0 for the even and perfect-power shortcuts.</summary>
    public long Base { get; }

    /// <summary>The period r of aˣ mod N found by the quantum step, or 0 if it was not used.</summary>
    public long Period { get; }

    /// <summary>The counting-register value measured in the successful attempt, or 0 if the quantum step was not used.</summary>
    public long Measurement { get; }

    /// <summary>Number of random bases tried, including the successful one. 0 for the even and perfect-power shortcuts.</summary>
    public int Attempts { get; }

    public ShorResult(long n, long factor1, long factor2, bool usedQuantumPeriodFinding,
        long @base, long period, long measurement, int attempts)
    {
        N = n;
        Factor1 = factor1 <= factor2 ? factor1 : factor2;
        Factor2 = factor1 <= factor2 ? factor2 : factor1;
        UsedQuantumPeriodFinding = usedQuantumPeriodFinding;
        Base = @base;
        Period = period;
        Measurement = measurement;
        Attempts = attempts;
    }

    public override string ToString() => $"{N} = {Factor1} × {Factor2}";
}
