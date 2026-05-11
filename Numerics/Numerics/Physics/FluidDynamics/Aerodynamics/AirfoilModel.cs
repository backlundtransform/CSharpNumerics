using System;

namespace CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

/// <summary>
/// Provides lift and drag coefficient lookup as a function of angle of attack (AoA).
/// Supports built-in profiles (flat plate, NACA symmetric) and custom lookup tables
/// with linear interpolation.
/// 
/// Coordinate conventions:
///   α > 0  ⇒  nose up, positive lift
///   Cl positive = upward lift, Cd positive = drag opposing motion
/// </summary>
public class AirfoilModel
{
    private readonly double[] _alphas;  // sorted AoA breakpoints in radians
    private readonly double[] _cls;     // Cl at each breakpoint
    private readonly double[] _cds;     // Cd at each breakpoint

    /// <summary>Display name of this airfoil profile.</summary>
    public string Name { get; }

    /// <summary>
    /// Creates a custom airfoil model from lookup tables.
    /// Arrays must be the same length and <paramref name="alphas"/> must be sorted ascending.
    /// </summary>
    /// <param name="name">Profile name.</param>
    /// <param name="alphas">Angle of attack breakpoints in radians (sorted ascending).</param>
    /// <param name="cls">Lift coefficients at each breakpoint.</param>
    /// <param name="cds">Drag coefficients at each breakpoint.</param>
    public AirfoilModel(string name, double[] alphas, double[] cls, double[] cds)
    {
        if (alphas.Length != cls.Length || alphas.Length != cds.Length)
            throw new ArgumentException("Alpha, Cl, and Cd arrays must have the same length.");
        if (alphas.Length < 2)
            throw new ArgumentException("At least two data points are required.");

        Name = name;
        _alphas = alphas;
        _cls = cls;
        _cds = cds;
    }

    /// <summary>
    /// Interpolated lift coefficient at the given angle of attack.
    /// Clamps to the table endpoints for AoA outside the defined range.
    /// </summary>
    /// <param name="alpha">Angle of attack in radians.</param>
    public double Cl(double alpha) => Interpolate(_alphas, _cls, alpha);

    /// <summary>
    /// Interpolated drag coefficient at the given angle of attack.
    /// Clamps to the table endpoints for AoA outside the defined range.
    /// </summary>
    /// <param name="alpha">Angle of attack in radians.</param>
    public double Cd(double alpha) => Interpolate(_alphas, _cds, alpha);

    /// <summary>
    /// Lift-to-drag ratio at the given angle of attack.
    /// </summary>
    /// <param name="alpha">Angle of attack in radians.</param>
    public double LiftToDrag(double alpha)
    {
        double cd = Cd(alpha);
        return cd > 1e-12 ? Cl(alpha) / cd : 0.0;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Built-in profiles
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Flat plate aerodynamic model.
    /// Cl ≈ 2π·sin(α)·cos(α), Cd ≈ 2·sin²(α) + Cd0.
    /// Simple analytical model useful for testing.
    /// </summary>
    /// <param name="cd0">Parasitic drag coefficient at zero AoA (default 0.02).</param>
    /// <param name="points">Number of table points to generate (default 91).</param>
    public static AirfoilModel FlatPlate(double cd0 = 0.02, int points = 91)
    {
        var alphas = new double[points];
        var cls = new double[points];
        var cds = new double[points];

        double aMin = -Math.PI / 2;
        double aMax = Math.PI / 2;
        double step = (aMax - aMin) / (points - 1);

        for (int i = 0; i < points; i++)
        {
            double a = aMin + i * step;
            alphas[i] = a;
            cls[i] = 2.0 * Math.PI * Math.Sin(a) * Math.Cos(a);
            cds[i] = 2.0 * Math.Sin(a) * Math.Sin(a) + cd0;
        }

        return new AirfoilModel("Flat Plate", alphas, cls, cds);
    }

    /// <summary>
    /// Simplified NACA 0012 symmetric airfoil model.
    /// Linear Cl region up to stall (≈15°), then Cl drops.
    /// Cd modelled as parabolic drag polar: Cd = Cd0 + k·Cl².
    /// </summary>
    /// <param name="clAlpha">Lift curve slope dCl/dα in 1/rad (default 2π ≈ 6.283).</param>
    /// <param name="alphaStall">Stall angle in radians (default 15° ≈ 0.2618).</param>
    /// <param name="clMax">Maximum Cl at stall (default 1.5).</param>
    /// <param name="cd0">Zero-lift drag coefficient (default 0.008).</param>
    /// <param name="k">Induced drag factor (default 0.04).</param>
    /// <param name="points">Number of table points (default 181).</param>
    public static AirfoilModel NACASymmetric(
        double clAlpha = 2.0 * Math.PI,
        double alphaStall = 0.2618,
        double clMax = 1.5,
        double cd0 = 0.008,
        double k = 0.04,
        int points = 181)
    {
        var alphas = new double[points];
        var cls = new double[points];
        var cds = new double[points];

        double aMin = -Math.PI / 2;
        double aMax = Math.PI / 2;
        double step = (aMax - aMin) / (points - 1);

        for (int i = 0; i < points; i++)
        {
            double a = aMin + i * step;
            alphas[i] = a;

            double absA = Math.Abs(a);
            double sign = Math.Sign(a);

            if (absA <= alphaStall)
            {
                // Linear region
                cls[i] = clAlpha * a;
            }
            else
            {
                // Post-stall: Cl drops from clMax following a sinusoidal decay
                // toward the flat-plate value at 90°. Models the separated flow regime.
                double range = Math.PI / 2 - alphaStall;
                double postStallFrac = (absA - alphaStall) / range;
                // Flat plate contribution at this AoA
                double clFlat = 2.0 * Math.Sin(absA) * Math.Cos(absA);
                // Decay from clMax toward flat-plate value using cosine blend
                double decay = 0.5 * (1.0 + Math.Cos(Math.PI * postStallFrac));
                cls[i] = sign * (clMax * decay + clFlat * (1.0 - decay));
            }

            cds[i] = cd0 + k * cls[i] * cls[i];
        }

        return new AirfoilModel("NACA Symmetric", alphas, cls, cds);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Helpers
    // ═══════════════════════════════════════════════════════════════

    private static double Interpolate(double[] xs, double[] ys, double x)
    {
        if (x <= xs[0]) return ys[0];
        if (x >= xs[xs.Length - 1]) return ys[ys.Length - 1];

        // Binary search for interval
        int lo = 0, hi = xs.Length - 1;
        while (hi - lo > 1)
        {
            int mid = (lo + hi) / 2;
            if (xs[mid] <= x) lo = mid;
            else hi = mid;
        }

        double t = (x - xs[lo]) / (xs[hi] - xs[lo]);
        return ys[lo] + t * (ys[hi] - ys[lo]);
    }
}
