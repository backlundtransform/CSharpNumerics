using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics;
using System;

namespace CSharpNumerics.Physics.FluidDynamics;

/// <summary>
/// Extension methods for vorticity analysis and stream functions:
/// vorticity (curl of velocity), enstrophy, helicity density,
/// and velocity field recovery from a 2-D stream function.
/// </summary>
public static class VorticityExtensions
{
    // ═══════════════════════════════════════════════════════════════
    //  Vorticity & circulation
    //
    //  ω = ∇ × v
    // ═══════════════════════════════════════════════════════════════

    #region Vorticity

    /// <summary>
    /// Computes the vorticity ω = ∇ × v at a point.
    /// </summary>
    /// <param name="velocity">The velocity vector field v(r).</param>
    /// <param name="point">The point at which to evaluate.</param>
    public static Vector Vorticity(
        this VectorField velocity, (double, double, double) point)
    {
        return velocity.Curl(point);
    }

    /// <summary>
    /// Computes the enstrophy density ½|ω|² at a point, a measure of
    /// rotational kinetic energy per unit volume (divided by ρ).
    /// </summary>
    /// <param name="velocity">The velocity vector field v(r).</param>
    /// <param name="point">The point at which to evaluate.</param>
    public static double Enstrophy(
        this VectorField velocity, (double, double, double) point)
    {
        var omega = velocity.Curl(point);
        return 0.5 * omega.Dot(omega);
    }

    /// <summary>
    /// Computes the helicity density v·ω at a point, quantifying the
    /// degree of linkage between velocity and vorticity.
    /// </summary>
    /// <param name="velocity">The velocity vector field v(r).</param>
    /// <param name="point">The point at which to evaluate.</param>
    public static double HelicityDensity(
        this VectorField velocity, (double, double, double) point)
    {
        var v = new Vector(velocity.fx(new Vector(point)),
                           velocity.fy(new Vector(point)),
                           velocity.fz(new Vector(point)));
        var omega = velocity.Curl(point);
        return v.Dot(omega);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Stream function (2-D incompressible)
    //
    //  vx = ∂ψ/∂y,  vy = −∂ψ/∂x
    // ═══════════════════════════════════════════════════════════════

    #region Stream Function

    /// <summary>
    /// Creates a 2-D velocity vector field from a stream function ψ(x, y):
    /// vx = ∂ψ/∂y, vy = −∂ψ/∂x, vz = 0.
    /// The resulting field is automatically divergence-free.
    /// </summary>
    /// <param name="streamFunction">The stream function ψ(r), depending on x and y.</param>
    public static VectorField VelocityFromStreamFunction(this ScalarField streamFunction)
    {
        // vx = ∂ψ/∂y
        Func<Vector, double> vx = r =>
        {
            var sf = streamFunction.f;
            return (sf(new Vector(r.x, r.y + 1e-6, r.z)) - sf(new Vector(r.x, r.y - 1e-6, r.z))) / 2e-6;
        };

        // vy = −∂ψ/∂x
        Func<Vector, double> vy = r =>
        {
            var sf = streamFunction.f;
            return -(sf(new Vector(r.x + 1e-6, r.y, r.z)) - sf(new Vector(r.x - 1e-6, r.y, r.z))) / 2e-6;
        };

        return new VectorField(vx, vy, _ => 0.0);
    }

    #endregion
}
