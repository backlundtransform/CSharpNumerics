using System;

namespace CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

/// <summary>
/// Compressible drag model for rocket-shaped bodies.
/// Returns drag coefficient Cd as a function of Mach number using a piecewise model:
///   - Subsonic (M &lt; 0.8): constant Cd₀
///   - Transonic (0.8–1.2): sharp rise to peak (wave drag onset)
///   - Supersonic (M &gt; 1.2): Prandtl-Glauert decay Cd ~ 1/√(M²-1)
/// 
/// This models the characteristic "drag divergence" near Mach 1 observed in all
/// blunt and slender bodies in compressible flow.
/// </summary>
public class CompressibleDragModel
{
    /// <summary>Subsonic drag coefficient (M &lt; 0.8).</summary>
    public double Cd0 { get; }

    /// <summary>Peak drag coefficient at transonic (M ≈ 1.0–1.2). Typically 2–4× Cd0.</summary>
    public double CdPeak { get; }

    /// <summary>Mach number at which peak Cd occurs (default 1.1).</summary>
    public double MachPeak { get; }

    /// <summary>
    /// Creates a compressible drag model.
    /// </summary>
    /// <param name="cd0">Subsonic drag coefficient (default 0.3 for typical rocket).</param>
    /// <param name="cdPeak">Peak transonic drag coefficient (default 0.8).</param>
    /// <param name="machPeak">Mach at peak drag (default 1.1).</param>
    public CompressibleDragModel(double cd0 = 0.3, double cdPeak = 0.8, double machPeak = 1.1)
    {
        if (cd0 <= 0) throw new ArgumentOutOfRangeException(nameof(cd0));
        if (cdPeak <= 0) throw new ArgumentOutOfRangeException(nameof(cdPeak));
        if (machPeak <= 0) throw new ArgumentOutOfRangeException(nameof(machPeak));

        Cd0 = cd0;
        CdPeak = cdPeak;
        MachPeak = machPeak;
    }

    /// <summary>
    /// Evaluates drag coefficient at the given Mach number.
    /// </summary>
    /// <param name="mach">Mach number (≥ 0).</param>
    /// <returns>Drag coefficient Cd.</returns>
    public double Evaluate(double mach)
    {
        if (mach < 0) mach = 0;

        // Subsonic regime: constant Cd0
        if (mach < 0.8)
            return Cd0;

        // Transonic rise: smooth ramp from Cd0 at M=0.8 to CdPeak at MachPeak
        if (mach <= MachPeak)
        {
            double t = (mach - 0.8) / (MachPeak - 0.8);
            // Smooth cubic interpolation (Hermite)
            double s = t * t * (3.0 - 2.0 * t);
            return Cd0 + (CdPeak - Cd0) * s;
        }

        // Supersonic decay: Prandtl-Glauert approximation
        // Cd ≈ CdPeak / √(M² - 1) normalized so that Cd(MachPeak) = CdPeak
        double betaPeak = Math.Sqrt(MachPeak * MachPeak - 1.0);
        double beta = Math.Sqrt(mach * mach - 1.0);
        return CdPeak * betaPeak / beta;
    }

    /// <summary>
    /// Creates a drag model typical for a slender rocket body (fineness ratio ~10-15).
    /// </summary>
    public static CompressibleDragModel SlenderRocket()
        => new CompressibleDragModel(cd0: 0.25, cdPeak: 0.7, machPeak: 1.05);

    /// <summary>
    /// Creates a drag model for a blunt capsule (Apollo/Soyuz-style reentry shape).
    /// </summary>
    public static CompressibleDragModel BluntCapsule()
        => new CompressibleDragModel(cd0: 0.4, cdPeak: 1.0, machPeak: 1.15);
}
