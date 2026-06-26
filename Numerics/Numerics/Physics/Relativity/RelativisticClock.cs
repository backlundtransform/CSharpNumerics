using System;
using CSharpNumerics.Physics.Constants;

namespace CSharpNumerics.Physics.Relativity;

/// <summary>
/// Combined special- and general-relativistic clock rates in the weak field —
/// the effect that must be corrected for in satellite navigation. A moving clock
/// runs slow (special relativity, −v²/2c²); a clock higher in a gravitational
/// potential runs fast (general relativity, +ΔΦ/c²). For a GPS satellite the net
/// effect is about +38 microseconds per day.
/// </summary>
public static class RelativisticClock
{
    private const double G = PhysicsConstants.GravitationalConstant;
    private const double C = PhysicsConstants.SpeedOfLight;
    private const double C2 = C * C;

    /// <summary>
    /// Weak-field gravitational fractional rate by which a static clock at
    /// <paramref name="upperRadius"/> runs faster than one at
    /// <paramref name="lowerRadius"/>: ΔΦ/c² ≈ (GM/c²)·(1/r_lower − 1/r_upper).
    /// A positive value means the upper clock ticks faster.
    /// </summary>
    /// <param name="mass">Central mass (kg).</param>
    /// <param name="lowerRadius">Radius of the lower (deeper) clock (m).</param>
    /// <param name="upperRadius">Radius of the higher clock (m).</param>
    public static double GravitationalFractionalRate(double mass, double lowerRadius, double upperRadius)
    {
        if (lowerRadius <= 0 || upperRadius <= 0)
            throw new ArgumentOutOfRangeException(nameof(lowerRadius), "Radii must be positive.");

        return G * mass / C2 * (1.0 / lowerRadius - 1.0 / upperRadius);
    }

    /// <summary>
    /// Velocity (special-relativistic) fractional rate by which a moving clock
    /// runs slow relative to a static one: −v²/(2c²) to leading order.
    /// </summary>
    public static double VelocityFractionalRate(double speed)
    {
        return -(speed * speed) / (2.0 * C2);
    }

    /// <summary>
    /// Net fractional rate of an orbiting clock relative to one on the ground,
    /// combining the gravitational blueshift (it is higher up) and the
    /// special-relativistic time dilation (it is moving):
    /// (GM/c²)(1/r_ground − 1/r_orbit) − v²/(2c²).
    /// Positive means the orbiting clock runs fast. Multiply by elapsed time to
    /// get the accumulated offset (≈ +38 µs/day for GPS).
    /// </summary>
    /// <param name="mass">Mass of the central body (kg).</param>
    /// <param name="groundRadius">Radius of the ground clock from the centre (m).</param>
    /// <param name="orbitRadius">Orbital radius of the satellite from the centre (m).</param>
    /// <param name="orbitSpeed">Orbital speed of the satellite (m/s).</param>
    public static double OrbitingClockFractionalRate(double mass, double groundRadius, double orbitRadius, double orbitSpeed)
    {
        return GravitationalFractionalRate(mass, groundRadius, orbitRadius) + VelocityFractionalRate(orbitSpeed);
    }
}
