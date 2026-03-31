using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics;
using System;

namespace CSharpNumerics.Physics.Environmental;

/// <summary>
/// Extension methods for Fickian diffusion:
/// Fick's first law (flux), Fick's second law (rate of change),
/// and the analytical point-source solution.
/// </summary>
public static class DiffusionExtensions
{
    // ═══════════════════════════════════════════════════════════════
    //  Diffusion (Fick's laws)
    // ═══════════════════════════════════════════════════════════════

    #region Diffusion

    /// <summary>
    /// Fick's first law: computes the diffusive flux J = −D ∇C.
    /// Returns a <see cref="VectorField"/> representing the mass flux.
    /// </summary>
    /// <param name="concentration">Concentration field C(r).</param>
    /// <param name="diffusionCoefficient">Diffusion coefficient D in m²/s.</param>
    public static VectorField DiffusionFlux(
        this ScalarField concentration,
        double diffusionCoefficient)
    {
        var grad = concentration.GradientField();
        return new VectorField(
            r => -diffusionCoefficient * grad.fx(r),
            r => -diffusionCoefficient * grad.fy(r),
            r => -diffusionCoefficient * grad.fz(r));
    }

    /// <summary>
    /// Fick's second law: computes the rate of change of concentration due
    /// to isotropic diffusion: ∂C/∂t = D ∇²C.
    /// Returns a <see cref="ScalarField"/> representing the local time derivative.
    /// </summary>
    /// <param name="concentration">Concentration field C(r).</param>
    /// <param name="diffusionCoefficient">Diffusion coefficient D in m²/s.</param>
    public static ScalarField DiffusionRate(
        this ScalarField concentration,
        double diffusionCoefficient)
    {
        var func = concentration.f;
        return new ScalarField(r =>
            diffusionCoefficient * func.Laplacian((r.x, r.y, r.z)));
    }

    /// <summary>
    /// Analytical solution for isotropic 3D diffusion from an instantaneous
    /// point source of mass M released at time t = 0:
    /// C(r, t) = M / (4πDt)^(3/2) · exp(−|r − r₀|² / (4Dt)).
    /// </summary>
    /// <param name="mass">Total mass M released (kg).</param>
    /// <param name="diffusionCoefficient">Diffusion coefficient D in m²/s.</param>
    /// <param name="time">Time since release in seconds (must be &gt; 0).</param>
    /// <param name="sourcePosition">Position r₀ of the release.</param>
    public static ScalarField DiffusionPointSource(
        this double mass,
        double diffusionCoefficient,
        double time,
        Vector sourcePosition)
    {
        if (time <= 0) throw new ArgumentException("Time must be greater than zero.");
        if (diffusionCoefficient <= 0) throw new ArgumentException("Diffusion coefficient must be greater than zero.");

        double D = diffusionCoefficient;
        double t = time;
        double coeff = mass / Math.Pow(4 * Math.PI * D * t, 1.5);
        double denom = 4 * D * t;

        return new ScalarField(r =>
        {
            var d = r - sourcePosition;
            double r2 = d.x * d.x + d.y * d.y + d.z * d.z;
            return coeff * Math.Exp(-r2 / denom);
        });
    }

    #endregion
}
