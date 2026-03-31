using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics;
using System;

namespace CSharpNumerics.Physics.Environmental;

/// <summary>
/// Extension methods for advective and combined advection–diffusion transport:
/// advection rate, the advection–diffusion equation, and the Péclet number.
/// </summary>
public static class TransportExtensions
{
    // ═══════════════════════════════════════════════════════════════
    //  Advection
    // ═══════════════════════════════════════════════════════════════

    #region Advection

    /// <summary>
    /// Computes the advective rate of change of concentration: ∂C/∂t = −v · ∇C.
    /// Returns a <see cref="ScalarField"/> representing the local time derivative
    /// due to advective transport by the velocity field.
    /// </summary>
    /// <param name="concentration">Concentration field C(r).</param>
    /// <param name="velocity">Velocity (or wind) field v(r).</param>
    public static ScalarField AdvectionRate(
        this ScalarField concentration,
        VectorField velocity)
    {
        var func = concentration.f;
        return new ScalarField(r =>
        {
            var grad = func.Gradient((r.x, r.y, r.z));
            return -(velocity.fx(r) * grad.x + velocity.fy(r) * grad.y + velocity.fz(r) * grad.z);
        });
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Advection–diffusion
    // ═══════════════════════════════════════════════════════════════

    #region Advection–Diffusion

    /// <summary>
    /// Computes the combined advection–diffusion rate of change:
    /// ∂C/∂t = D ∇²C − v · ∇C (+ S if a source term is provided).
    /// <para>
    /// Bridges <see cref="ScalarField.Laplacian"/> (diffusion) and
    /// <see cref="VectorFieldExtensions.Gradient(Func{Vector, double}, ValueTuple{double, double, double})"/> (advection)
    /// into a single transport operator.
    /// </para>
    /// </summary>
    /// <param name="concentration">Concentration field C(r).</param>
    /// <param name="velocity">Velocity (or wind) field v(r).</param>
    /// <param name="diffusionCoefficient">Diffusion coefficient D in m²/s.</param>
    /// <param name="source">Optional volumetric source/sink term S(r) in kg/(m³·s).</param>
    public static ScalarField AdvectionDiffusionRate(
        this ScalarField concentration,
        VectorField velocity,
        double diffusionCoefficient,
        ScalarField source = default)
    {
        var func = concentration.f;
        var hasSource = source.f != null;
        return new ScalarField(r =>
        {
            double laplacian = func.Laplacian((r.x, r.y, r.z));
            var grad = func.Gradient((r.x, r.y, r.z));
            double advection = velocity.fx(r) * grad.x + velocity.fy(r) * grad.y + velocity.fz(r) * grad.z;
            double result = diffusionCoefficient * laplacian - advection;

            if (hasSource) result += source.f(r);

            return result;
        });
    }

    /// <summary>
    /// Computes the Péclet number: Pe = u·L / D.
    /// <para>
    /// Pe ≫ 1 → advection-dominated transport.
    /// Pe ≪ 1 → diffusion-dominated transport.
    /// </para>
    /// </summary>
    /// <param name="velocity">Characteristic velocity in m/s.</param>
    /// <param name="characteristicLength">Length scale L in metres.</param>
    /// <param name="diffusionCoefficient">Diffusion coefficient D in m²/s.</param>
    public static double PecletNumber(this double velocity, double characteristicLength, double diffusionCoefficient)
    {
        if (diffusionCoefficient <= 0) throw new ArgumentException("Diffusion coefficient must be greater than zero.");
        return velocity * characteristicLength / diffusionCoefficient;
    }

    #endregion
}
