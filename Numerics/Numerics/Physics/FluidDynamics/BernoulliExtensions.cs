using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics;
using System;

namespace CSharpNumerics.Physics.FluidDynamics;

/// <summary>
/// Extension methods for Bernoulli's principle and the continuity equation:
/// Bernoulli constant, dynamic/stagnation pressure, mass flux,
/// volume flow rate, and continuity speed.
/// </summary>
public static class BernoulliExtensions
{
    // ═══════════════════════════════════════════════════════════════
    //  Bernoulli's principle
    //
    //  p + ½ρv² + ρgh = const   (along a streamline, steady, inviscid)
    // ═══════════════════════════════════════════════════════════════

    #region Bernoulli

    /// <summary>
    /// Computes Bernoulli's constant along a streamline:
    /// B = p + ½ρv² + ρgh.
    /// </summary>
    /// <param name="pressure">Static pressure at the point in Pa.</param>
    /// <param name="density">Fluid density ρ in kg/m³.</param>
    /// <param name="velocity">Flow speed at the point in m/s.</param>
    /// <param name="height">Height above reference datum in m.</param>
    public static double BernoulliConstant(
        this double pressure, double density, double velocity, double height = 0)
    {
        return pressure
            + 0.5 * density * velocity * velocity
            + density * PhysicsConstants.GravitationalAcceleration * height;
    }

    /// <summary>
    /// Computes the pressure at a second point along a streamline using
    /// Bernoulli's equation, given conditions at a first point.
    /// p₂ = p₁ + ½ρ(v₁² − v₂²) + ρg(h₁ − h₂).
    /// </summary>
    /// <param name="p1">Pressure at point 1 in Pa.</param>
    /// <param name="density">Fluid density ρ in kg/m³.</param>
    /// <param name="v1">Speed at point 1 in m/s.</param>
    /// <param name="v2">Speed at point 2 in m/s.</param>
    /// <param name="h1">Height at point 1 in m.</param>
    /// <param name="h2">Height at point 2 in m.</param>
    public static double BernoulliPressure(
        this double p1, double density, double v1, double v2,
        double h1 = 0, double h2 = 0)
    {
        return p1
            + 0.5 * density * (v1 * v1 - v2 * v2)
            + density * PhysicsConstants.GravitationalAcceleration * (h1 - h2);
    }

    /// <summary>
    /// Computes the dynamic pressure: q = ½ρv².
    /// </summary>
    /// <param name="density">Fluid density ρ in kg/m³.</param>
    /// <param name="speed">Flow speed in m/s.</param>
    public static double DynamicPressure(this double density, double speed)
    {
        return 0.5 * density * speed * speed;
    }

    /// <summary>
    /// Computes the stagnation (total) pressure: p₀ = p + ½ρv².
    /// </summary>
    /// <param name="staticPressure">Static pressure in Pa.</param>
    /// <param name="density">Fluid density ρ in kg/m³.</param>
    /// <param name="speed">Flow speed in m/s.</param>
    public static double StagnationPressure(this double staticPressure, double density, double speed)
    {
        return staticPressure + 0.5 * density * speed * speed;
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Continuity equation
    //
    //  ∂ρ/∂t + ∇·(ρv) = 0
    //  For incompressible flow: ∇·v = 0
    // ═══════════════════════════════════════════════════════════════

    #region Continuity

    /// <summary>
    /// Computes the mass flux ρv as a vector field from a constant-density
    /// velocity field.
    /// </summary>
    /// <param name="velocity">The velocity vector field v(r).</param>
    /// <param name="density">Fluid density ρ in kg/m³.</param>
    public static VectorField MassFlux(this VectorField velocity, double density)
    {
        return new VectorField(
            r => density * velocity.fx(r),
            r => density * velocity.fy(r),
            r => density * velocity.fz(r));
    }

    /// <summary>
    /// Computes the mass flux ρv as a vector field with spatially varying density.
    /// </summary>
    /// <param name="velocity">The velocity vector field v(r).</param>
    /// <param name="density">The density scalar field ρ(r).</param>
    public static VectorField MassFlux(this VectorField velocity, ScalarField density)
    {
        return new VectorField(
            r => density.f(r) * velocity.fx(r),
            r => density.f(r) * velocity.fy(r),
            r => density.f(r) * velocity.fz(r));
    }

    /// <summary>
    /// Evaluates the continuity equation residual ∇·(ρv) for constant density.
    /// For incompressible flow this equals ρ∇·v, which should be zero.
    /// </summary>
    /// <param name="velocity">The velocity vector field v(r).</param>
    /// <param name="density">Fluid density ρ in kg/m³.</param>
    /// <param name="point">The point at which to evaluate.</param>
    public static double ContinuityResidual(
        this VectorField velocity, double density, (double, double, double) point)
    {
        return velocity.MassFlux(density).Divergence(point);
    }

    /// <summary>
    /// Computes the volume flow rate Q = A · v through a cross-section
    /// of known area, given a uniform speed perpendicular to the section.
    /// </summary>
    /// <param name="area">Cross-sectional area in m².</param>
    /// <param name="speed">Average flow speed normal to the section in m/s.</param>
    public static double VolumeFlowRate(this double area, double speed)
    {
        return area * speed;
    }

    /// <summary>
    /// Computes the speed at a second cross-section from the continuity
    /// equation for incompressible flow: A₁v₁ = A₂v₂ → v₂ = A₁v₁ / A₂.
    /// </summary>
    /// <param name="area1">Cross-sectional area at section 1 in m².</param>
    /// <param name="speed1">Speed at section 1 in m/s.</param>
    /// <param name="area2">Cross-sectional area at section 2 in m².</param>
    public static double ContinuitySpeed(this double area1, double speed1, double area2)
    {
        if (area2 <= 0) throw new ArgumentException("Cross-sectional area must be greater than zero.");
        return area1 * speed1 / area2;
    }

    #endregion
}
