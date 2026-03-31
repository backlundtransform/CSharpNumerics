using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics;
using System;

namespace CSharpNumerics.Physics.Electromagnetism;

/// <summary>
/// Extension methods for Maxwell's equations (differential form),
/// the electromagnetic wave equation, and field-from-potential calculations.
/// </summary>
public static class MaxwellExtensions
{
    // ═══════════════════════════════════════════════════════════════
    //  Maxwell's equations — differential form
    //
    //  (1) ∇·E = ρ/ε₀          Gauss (electric)
    //  (2) ∇·B = 0              Gauss (magnetic)
    //  (3) ∇×E = −∂B/∂t        Faraday
    //  (4) ∇×B = μ₀J + μ₀ε₀∂E/∂t   Ampère–Maxwell
    // ═══════════════════════════════════════════════════════════════

    #region Maxwell's Equations

    /// <summary>
    /// Gauss's law for electricity: computes the charge density at a point
    /// from the electric field divergence.
    /// ρ = ε₀ ∇·E.
    /// </summary>
    /// <param name="electricField">The electric vector field.</param>
    /// <param name="point">The point at which to evaluate.</param>
    public static double ChargeDensity(this VectorField electricField, (double, double, double) point)
    {
        return PhysicsConstants.VacuumPermittivity * electricField.Divergence(point);
    }

    /// <summary>
    /// Gauss's law for magnetism: computes ∇·B which should be zero
    /// for any physical magnetic field. Returns the numerical divergence
    /// for verification.
    /// </summary>
    /// <param name="magneticField">The magnetic vector field.</param>
    /// <param name="point">The point at which to evaluate.</param>
    public static double GaussMagnetic(this VectorField magneticField, (double, double, double) point)
    {
        return magneticField.Divergence(point);
    }

    /// <summary>
    /// Faraday's law of induction: computes −∇×E, which equals ∂B/∂t.
    /// For electrostatic fields, the result should be approximately zero.
    /// </summary>
    /// <param name="electricField">The electric vector field.</param>
    /// <param name="point">The point at which to evaluate.</param>
    /// <returns>The time rate of change of the magnetic field: ∂B/∂t = −∇×E.</returns>
    public static Vector FaradayLaw(this VectorField electricField, (double, double, double) point)
    {
        var curl = electricField.Curl(point);
        return new Vector(-curl.x, -curl.y, -curl.z);
    }

    /// <summary>
    /// Ampère–Maxwell law: computes (1/μ₀)∇×B, which equals J + ε₀ ∂E/∂t.
    /// For magnetostatic fields with no changing E, the result gives
    /// the current density J.
    /// </summary>
    /// <param name="magneticField">The magnetic vector field.</param>
    /// <param name="point">The point at which to evaluate.</param>
    /// <returns>J + ε₀ ∂E/∂t = (1/μ₀)∇×B.</returns>
    public static Vector AmpereLaw(this VectorField magneticField, (double, double, double) point)
    {
        var curl = magneticField.Curl(point);
        var invMu = 1.0 / PhysicsConstants.VacuumPermeability;
        return invMu * curl;
    }

    /// <summary>
    /// Verifies the wave equation for a scalar component of an electromagnetic field:
    /// returns ∇²f − μ₀ε₀ ∂²f/∂t², which should be zero in free space.
    /// </summary>
    /// <param name="field">A scalar field component f(r, t) parameterised as f(t)(r).</param>
    /// <param name="t">Time at which to evaluate.</param>
    /// <param name="point">Spatial point at which to evaluate.</param>
    /// <param name="dt">Time step for numerical second derivative (default 1 × 10⁻⁶).</param>
    public static double WaveEquationResidual(
        this Func<double, Func<Vector, double>> field,
        double t, (double, double, double) point, double dt = 1e-6)
    {
        var laplacian = field(t).Laplacian(point);

        var fPlus = field(t + dt)(new Vector(point));
        var fMid = field(t)(new Vector(point));
        var fMinus = field(t - dt)(new Vector(point));
        var d2fdt2 = (fPlus - 2 * fMid + fMinus) / (dt * dt);

        var mu0eps0 = PhysicsConstants.VacuumPermeability * PhysicsConstants.VacuumPermittivity;
        return laplacian - mu0eps0 * d2fdt2;
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Potentials
    // ═══════════════════════════════════════════════════════════════

    #region Potentials

    /// <summary>
    /// Creates an electric vector field from a scalar potential: E = −∇V.
    /// The returned <see cref="VectorField"/> can be used with ∇· and ∇× operators.
    /// </summary>
    /// <param name="potential">The electric potential V(r).</param>
    public static VectorField ElectricFieldFromPotential(this Func<Vector, double> potential)
    {
        return (-new ScalarField(potential)).GradientField();
    }

    /// <summary>
    /// Creates an electric vector field from a <see cref="ScalarField"/> potential: E = −∇V.
    /// </summary>
    /// <param name="potential">The electric potential V(r) as a ScalarField.</param>
    public static VectorField ElectricFieldFromPotential(this ScalarField potential)
    {
        return (-potential).GradientField();
    }

    /// <summary>
    /// Computes the magnetic field from a vector potential at a point: B = ∇ × A.
    /// </summary>
    /// <param name="vectorPotential">The magnetic vector potential A(r).</param>
    /// <param name="point">The point at which to evaluate.</param>
    public static Vector MagneticFieldFromVectorPotential(
        this VectorField vectorPotential, (double, double, double) point)
    {
        return vectorPotential.Curl(point);
    }

    #endregion
}
