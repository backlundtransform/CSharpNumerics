using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.FluidDynamics;

/// <summary>
/// Extension methods for viscous/pipe flow, dimensionless numbers,
/// hydrostatics, and fluid energy/momentum density.
/// </summary>
public static class ViscousFlowExtensions
{
    // ═══════════════════════════════════════════════════════════════
    //  Dimensionless numbers
    // ═══════════════════════════════════════════════════════════════

    #region Dimensionless Numbers

    /// <summary>
    /// Computes the Reynolds number: Re = ρvL / μ = vL / ν.
    /// Characterises the ratio of inertial to viscous forces.
    /// </summary>
    /// <param name="density">Fluid density ρ in kg/m³.</param>
    /// <param name="speed">Characteristic flow speed in m/s.</param>
    /// <param name="characteristicLength">Characteristic length L in m.</param>
    /// <param name="dynamicViscosity">Dynamic viscosity μ in Pa·s.</param>
    public static double ReynoldsNumber(
        this double density, double speed, double characteristicLength, double dynamicViscosity)
    {
        if (dynamicViscosity <= 0) throw new ArgumentException("Dynamic viscosity must be positive.");
        return density * speed * characteristicLength / dynamicViscosity;
    }

    /// <summary>
    /// Computes the Mach number: Ma = v / c.
    /// </summary>
    /// <param name="speed">Flow speed in m/s.</param>
    /// <param name="speedOfSound">Speed of sound in the medium in m/s
    /// (defaults to <see cref="PhysicsConstants.SpeedOfSoundAir"/>).</param>
    public static double MachNumber(this double speed, double speedOfSound = 0)
    {
        if (speedOfSound <= 0) speedOfSound = PhysicsConstants.SpeedOfSoundAir;
        return speed / speedOfSound;
    }

    /// <summary>
    /// Computes the Froude number: Fr = v / √(gL).
    /// Characterises the ratio of flow inertia to gravity.
    /// </summary>
    /// <param name="speed">Flow speed in m/s.</param>
    /// <param name="characteristicLength">Characteristic length L in m (e.g. water depth).</param>
    public static double FroudeNumber(this double speed, double characteristicLength)
    {
        if (characteristicLength <= 0) throw new ArgumentException("Characteristic length must be positive.");
        return speed / Math.Sqrt(PhysicsConstants.GravitationalAcceleration * characteristicLength);
    }

    /// <summary>
    /// Computes the Strouhal number: St = fL / v.
    /// Characterises oscillating flow mechanisms (e.g. vortex shedding).
    /// </summary>
    /// <param name="frequency">Oscillation frequency f in Hz.</param>
    /// <param name="characteristicLength">Characteristic length L in m.</param>
    /// <param name="speed">Flow speed v in m/s.</param>
    public static double StrouhalNumber(this double frequency, double characteristicLength, double speed)
    {
        if (speed <= 0) throw new ArgumentException("Speed must be positive.");
        return frequency * characteristicLength / speed;
    }

    /// <summary>
    /// Computes the Weber number: We = ρv²L / σ.
    /// Characterises the ratio of inertia to surface tension.
    /// </summary>
    /// <param name="density">Fluid density ρ in kg/m³.</param>
    /// <param name="speed">Flow speed in m/s.</param>
    /// <param name="characteristicLength">Characteristic length L in m.</param>
    /// <param name="surfaceTension">Surface tension σ in N/m.</param>
    public static double WeberNumber(
        this double density, double speed, double characteristicLength, double surfaceTension)
    {
        if (surfaceTension <= 0) throw new ArgumentException("Surface tension must be positive.");
        return density * speed * speed * characteristicLength / surfaceTension;
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Viscous / pipe flow
    // ═══════════════════════════════════════════════════════════════

    #region Viscous Flow

    /// <summary>
    /// Computes the Hagen–Poiseuille volumetric flow rate through a
    /// circular pipe of radius R and length L under pressure drop ΔP:
    /// Q = πR⁴ΔP / (8μL).
    /// </summary>
    /// <param name="radius">Pipe inner radius R in m.</param>
    /// <param name="pressureDrop">Pressure drop ΔP across the pipe in Pa.</param>
    /// <param name="dynamicViscosity">Dynamic viscosity μ in Pa·s.</param>
    /// <param name="length">Pipe length L in m.</param>
    public static double PoiseuilleFlowRate(
        this double radius, double pressureDrop, double dynamicViscosity, double length)
    {
        if (radius <= 0 || dynamicViscosity <= 0 || length <= 0)
            throw new ArgumentException("Radius, viscosity, and length must be positive.");
        return Math.PI * Math.Pow(radius, 4) * pressureDrop / (8.0 * dynamicViscosity * length);
    }

    /// <summary>
    /// Computes the Poiseuille velocity profile in a circular pipe:
    /// v(r) = (ΔP / 4μL)(R² − r²), where r is the radial distance from centre.
    /// </summary>
    /// <param name="radius">Pipe inner radius R in m.</param>
    /// <param name="pressureDrop">Pressure drop ΔP in Pa.</param>
    /// <param name="dynamicViscosity">Dynamic viscosity μ in Pa·s.</param>
    /// <param name="length">Pipe length L in m.</param>
    /// <param name="radialDistance">Radial distance r from the pipe centre in m.</param>
    public static double PoiseuilleVelocity(
        this double radius, double pressureDrop, double dynamicViscosity,
        double length, double radialDistance)
    {
        if (radialDistance < 0 || radialDistance > radius)
            throw new ArgumentException("Radial distance must be between 0 and the pipe radius.");
        return pressureDrop / (4.0 * dynamicViscosity * length) * (radius * radius - radialDistance * radialDistance);
    }

    /// <summary>
    /// Computes the Stokes drag force on a sphere in creeping flow:
    /// F = 6πμRv.
    /// </summary>
    /// <param name="dynamicViscosity">Dynamic viscosity μ in Pa·s.</param>
    /// <param name="radius">Sphere radius R in m.</param>
    /// <param name="speed">Relative speed v in m/s.</param>
    public static double StokesDrag(this double dynamicViscosity, double radius, double speed)
    {
        return 6.0 * Math.PI * dynamicViscosity * radius * speed;
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Hydrostatics
    // ═══════════════════════════════════════════════════════════════

    #region Hydrostatics

    /// <summary>
    /// Computes the hydrostatic pressure at a depth below the surface:
    /// p = p₀ + ρgh.
    /// </summary>
    /// <param name="depth">Depth h below the surface in m.</param>
    /// <param name="density">Fluid density ρ in kg/m³.</param>
    /// <param name="surfacePressure">Pressure at the surface p₀ in Pa
    /// (defaults to <see cref="PhysicsConstants.StandardAtmosphericPressure"/>).</param>
    public static double HydrostaticPressure(
        this double depth, double density,
        double surfacePressure = 0)
    {
        if (surfacePressure <= 0) surfacePressure = PhysicsConstants.StandardAtmosphericPressure;
        return surfacePressure + density * PhysicsConstants.GravitationalAcceleration * depth;
    }

    /// <summary>
    /// Computes the buoyant force (Archimedes' principle):
    /// F_b = ρ_fluid · V_displaced · g.
    /// </summary>
    /// <param name="fluidDensity">Density of the surrounding fluid in kg/m³.</param>
    /// <param name="displacedVolume">Volume of fluid displaced in m³.</param>
    public static double BuoyantForce(this double fluidDensity, double displacedVolume)
    {
        return fluidDensity * displacedVolume * PhysicsConstants.GravitationalAcceleration;
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Kinetic energy & momentum
    // ═══════════════════════════════════════════════════════════════

    #region Fluid Energy

    /// <summary>
    /// Computes the kinetic energy density of a flow: e_k = ½ρ|v|².
    /// </summary>
    /// <param name="velocity">Velocity vector at the point.</param>
    /// <param name="density">Fluid density ρ in kg/m³.</param>
    public static double KineticEnergyDensity(this Vector velocity, double density)
    {
        return 0.5 * density * velocity.Dot(velocity);
    }

    /// <summary>
    /// Computes the momentum density of a flow: ρv.
    /// </summary>
    /// <param name="velocity">Velocity vector at the point.</param>
    /// <param name="density">Fluid density ρ in kg/m³.</param>
    public static Vector MomentumDensity(this Vector velocity, double density)
    {
        return density * velocity;
    }

    #endregion
}
