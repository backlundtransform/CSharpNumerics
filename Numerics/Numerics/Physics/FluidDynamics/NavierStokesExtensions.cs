using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics;
using System;

namespace CSharpNumerics.Physics.FluidDynamics;

/// <summary>
/// Extension methods for the Navier–Stokes and Euler equations:
/// convective acceleration, viscous term, pressure gradient,
/// full N-S residual, incompressibility check, and inviscid Euler residual.
/// </summary>
public static class NavierStokesExtensions
{
    // ═══════════════════════════════════════════════════════════════
    //  Navier–Stokes equations (incompressible)
    //
    //  ρ(∂v/∂t + (v·∇)v) = −∇p + μ∇²v + f
    //  ∇·v = 0  (incompressibility constraint)
    // ═══════════════════════════════════════════════════════════════

    #region Navier–Stokes

    /// <summary>
    /// Computes the convective acceleration term (v·∇)v of the Navier–Stokes
    /// equations at a point. This is the nonlinear advection term.
    /// </summary>
    /// <param name="velocity">The velocity vector field v(r).</param>
    /// <param name="point">The point at which to evaluate.</param>
    /// <param name="h">Step size for numerical differentiation (default 1×10⁻⁶).</param>
    public static Vector ConvectiveAcceleration(
        this VectorField velocity, (double, double, double) point, double h = 1e-6)
    {
        var p = new Vector(point);
        var vx = velocity.fx(p);
        var vy = velocity.fy(p);
        var vz = velocity.fz(p);

        // ∂vᵢ/∂xⱼ via central differences
        double Partial(Func<Vector, double> comp, int axis)
        {
            var pPlus = new Vector(
                point.Item1 + (axis == 0 ? h : 0),
                point.Item2 + (axis == 1 ? h : 0),
                point.Item3 + (axis == 2 ? h : 0));
            var pMinus = new Vector(
                point.Item1 - (axis == 0 ? h : 0),
                point.Item2 - (axis == 1 ? h : 0),
                point.Item3 - (axis == 2 ? h : 0));
            return (comp(pPlus) - comp(pMinus)) / (2.0 * h);
        }

        // (v·∇)vᵢ = vx ∂vᵢ/∂x + vy ∂vᵢ/∂y + vz ∂vᵢ/∂z
        double ax = vx * Partial(velocity.fx, 0) + vy * Partial(velocity.fx, 1) + vz * Partial(velocity.fx, 2);
        double ay = vx * Partial(velocity.fy, 0) + vy * Partial(velocity.fy, 1) + vz * Partial(velocity.fy, 2);
        double az = vx * Partial(velocity.fz, 0) + vy * Partial(velocity.fz, 1) + vz * Partial(velocity.fz, 2);

        return new Vector(ax, ay, az);
    }

    /// <summary>
    /// Computes the viscous diffusion term μ∇²v of the Navier–Stokes equations
    /// at a point, where ∇²v is the vector Laplacian applied component-wise.
    /// </summary>
    /// <param name="velocity">The velocity vector field v(r).</param>
    /// <param name="dynamicViscosity">Dynamic viscosity μ in Pa·s.</param>
    /// <param name="point">The point at which to evaluate.</param>
    /// <param name="h">Step size for numerical differentiation (default 1×10⁻⁶).</param>
    public static Vector ViscousTerm(
        this VectorField velocity, double dynamicViscosity,
        (double, double, double) point, double h = 1e-6)
    {
        double LaplacianComponent(Func<Vector, double> comp)
        {
            return new ScalarField(comp).Laplacian(point);
        }

        return dynamicViscosity * new Vector(
            LaplacianComponent(velocity.fx),
            LaplacianComponent(velocity.fy),
            LaplacianComponent(velocity.fz));
    }

    /// <summary>
    /// Computes the pressure gradient term −∇p of the Navier–Stokes equations.
    /// </summary>
    /// <param name="pressure">The pressure scalar field p(r).</param>
    /// <param name="point">The point at which to evaluate.</param>
    public static Vector PressureGradientForce(
        this ScalarField pressure, (double, double, double) point)
    {
        var grad = pressure.Gradient(point);
        return new Vector(-grad.x, -grad.y, -grad.z);
    }

    /// <summary>
    /// Computes the full Navier–Stokes acceleration (right-hand side / ρ)
    /// for an incompressible Newtonian fluid at a steady state:
    /// a = (−∇p + μ∇²v + f) / ρ − (v·∇)v.
    /// </summary>
    /// <param name="velocity">The velocity vector field v(r).</param>
    /// <param name="pressure">The pressure scalar field p(r).</param>
    /// <param name="density">Fluid density ρ in kg/m³.</param>
    /// <param name="dynamicViscosity">Dynamic viscosity μ in Pa·s.</param>
    /// <param name="point">The point at which to evaluate.</param>
    /// <param name="bodyForce">External body force per unit volume f (e.g. gravity). Pass zero vector if none.</param>
    public static Vector NavierStokesResidual(
        this VectorField velocity, ScalarField pressure,
        double density, double dynamicViscosity,
        (double, double, double) point, Vector bodyForce)
    {
        var pressureGrad = pressure.PressureGradientForce(point);
        var viscous = velocity.ViscousTerm(dynamicViscosity, point);
        var convective = velocity.ConvectiveAcceleration(point);

        return (1.0 / density) * (pressureGrad + viscous + bodyForce) - convective;
    }

    /// <summary>
    /// Computes the full Navier–Stokes acceleration with no external body force.
    /// </summary>
    public static Vector NavierStokesResidual(
        this VectorField velocity, ScalarField pressure,
        double density, double dynamicViscosity,
        (double, double, double) point)
    {
        return velocity.NavierStokesResidual(pressure, density, dynamicViscosity, point, new Vector(0, 0, 0));
    }

    /// <summary>
    /// Checks the incompressibility constraint ∇·v = 0.
    /// Returns the divergence of the velocity field, which should be
    /// zero for an incompressible flow.
    /// </summary>
    /// <param name="velocity">The velocity vector field v(r).</param>
    /// <param name="point">The point at which to evaluate.</param>
    public static double IncompressibilityResidual(
        this VectorField velocity, (double, double, double) point)
    {
        return velocity.Divergence(point);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Euler equations (inviscid flow)
    //
    //  ρ(∂v/∂t + (v·∇)v) = −∇p + f
    // ═══════════════════════════════════════════════════════════════

    #region Euler Equations

    /// <summary>
    /// Computes the Euler equation residual for inviscid flow:
    /// ∂v/∂t = (−∇p + f) / ρ − (v·∇)v.
    /// Equivalent to Navier–Stokes with μ = 0.
    /// </summary>
    /// <param name="velocity">The velocity vector field v(r).</param>
    /// <param name="pressure">The pressure scalar field p(r).</param>
    /// <param name="density">Fluid density ρ in kg/m³.</param>
    /// <param name="point">The point at which to evaluate.</param>
    /// <param name="bodyForce">External body force per unit volume f. Pass zero vector if none.</param>
    public static Vector EulerEquationResidual(
        this VectorField velocity, ScalarField pressure,
        double density, (double, double, double) point, Vector bodyForce)
    {
        var pressureGrad = pressure.PressureGradientForce(point);
        var convective = velocity.ConvectiveAcceleration(point);

        return (1.0 / density) * (pressureGrad + bodyForce) - convective;
    }

    /// <summary>
    /// Computes the Euler equation residual with no external body force.
    /// </summary>
    public static Vector EulerEquationResidual(
        this VectorField velocity, ScalarField pressure,
        double density, (double, double, double) point)
    {
        return velocity.EulerEquationResidual(pressure, density, point, new Vector(0, 0, 0));
    }

    #endregion
}
