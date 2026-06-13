using System;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Extension methods for quantum tunnelling through a one-dimensional rectangular potential
/// barrier of height V₀ and width a. Closed-form transmission/reflection coefficients for a
/// particle of energy E, covering the tunnelling regime (E &lt; V₀), the over-barrier regime
/// (E &gt; V₀, with resonances), and the E = V₀ limit.
/// </summary>
public static class QuantumTunnelingExtensions
{
    /// <summary>
    /// Transmission coefficient T ∈ (0, 1] for a rectangular barrier.
    /// <list type="bullet">
    /// <item>E &lt; V₀:  T = [1 + V₀²sinh²(κa) / (4E(V₀−E))]⁻¹,  κ = √(2m(V₀−E))/ħ</item>
    /// <item>E &gt; V₀:  T = [1 + V₀²sin²(ka) / (4E(E−V₀))]⁻¹,   k = √(2m(E−V₀))/ħ</item>
    /// <item>E = V₀:  T = [1 + mV₀a²/(2ħ²)]⁻¹</item>
    /// </list>
    /// </summary>
    /// <param name="energy">Particle energy E (&gt; 0).</param>
    /// <param name="barrierHeight">Barrier height V₀.</param>
    /// <param name="barrierWidth">Barrier width a (&gt; 0).</param>
    /// <param name="mass">Particle mass (default 1).</param>
    /// <param name="hbar">Reduced Planck constant (default 1).</param>
    public static double RectangularBarrierTransmission(
        this double energy, double barrierHeight, double barrierWidth,
        double mass = 1.0, double hbar = 1.0)
    {
        if (energy <= 0) throw new ArgumentOutOfRangeException(nameof(energy), "Energy must be positive.");
        if (barrierWidth <= 0) throw new ArgumentOutOfRangeException(nameof(barrierWidth));
        if (mass <= 0) throw new ArgumentOutOfRangeException(nameof(mass));
        if (hbar <= 0) throw new ArgumentOutOfRangeException(nameof(hbar));

        // No barrier → full transmission.
        if (barrierHeight <= 0) return 1.0;

        double a = barrierWidth;

        if (Math.Abs(energy - barrierHeight) < 1e-15)
        {
            // E = V₀ limit.
            double term = (mass * barrierHeight * a * a) / (2.0 * hbar * hbar);
            return 1.0 / (1.0 + term);
        }

        if (energy < barrierHeight)
        {
            double kappa = Math.Sqrt(2.0 * mass * (barrierHeight - energy)) / hbar;
            double sinh = Math.Sinh(kappa * a);
            double denom = 1.0 + (barrierHeight * barrierHeight * sinh * sinh)
                / (4.0 * energy * (barrierHeight - energy));
            return 1.0 / denom;
        }
        else
        {
            double k = Math.Sqrt(2.0 * mass * (energy - barrierHeight)) / hbar;
            double sin = Math.Sin(k * a);
            double denom = 1.0 + (barrierHeight * barrierHeight * sin * sin)
                / (4.0 * energy * (energy - barrierHeight));
            return 1.0 / denom;
        }
    }

    /// <summary>
    /// Reflection coefficient R = 1 − T for a rectangular barrier.
    /// </summary>
    public static double RectangularBarrierReflection(
        this double energy, double barrierHeight, double barrierWidth,
        double mass = 1.0, double hbar = 1.0)
        => 1.0 - energy.RectangularBarrierTransmission(barrierHeight, barrierWidth, mass, hbar);
}
