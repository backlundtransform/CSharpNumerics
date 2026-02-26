namespace CSharpNumerics.Physics.Oscillations
{
    /// <summary>
    /// Describes the damping regime of a damped oscillator.
    /// </summary>
    public enum DampingRegime
    {
        /// <summary>γ &lt; ω₀ — oscillates with exponentially decaying amplitude.</summary>
        Underdamped,

        /// <summary>γ = ω₀ — returns to equilibrium as fast as possible without oscillating.</summary>
        CriticallyDamped,

        /// <summary>γ &gt; ω₀ — returns to equilibrium without oscillating, slower than critical.</summary>
        Overdamped
    }
}
