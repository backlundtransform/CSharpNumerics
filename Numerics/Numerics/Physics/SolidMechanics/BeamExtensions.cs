using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.FiniteDifference;
using System;

namespace CSharpNumerics.Physics.SolidMechanics;

/// <summary>
/// Extension methods for the Euler–Bernoulli beam equation and
/// analytical beam deflection formulas for common configurations
/// (cantilever, simply supported).
/// </summary>
public static class BeamExtensions
{
    // ═══════════════════════════════════════════════════════════════
    //  Euler–Bernoulli Beam Equation: EIu⁗ = q
    // ═══════════════════════════════════════════════════════════════

    #region Euler-Bernoulli Beam Equation

    /// <summary>
    /// Bending moment from curvature: M = E · I · κ.
    /// </summary>
    /// <param name="youngsModulus">Young's modulus E in Pa.</param>
    /// <param name="secondMoment">Second moment of area I in m⁴.</param>
    /// <param name="curvature">Beam curvature κ = d²u/dx² in 1/m.</param>
    /// <returns>Bending moment M in N·m.</returns>
    public static double BendingMoment(this double youngsModulus, double secondMoment, double curvature)
    {
        return youngsModulus * secondMoment * curvature;
    }

    /// <summary>
    /// Bending stress from the flexure formula: σ = M · y / I.
    /// </summary>
    /// <param name="bendingMoment">Bending moment M in N·m.</param>
    /// <param name="distanceFromNeutralAxis">Distance y from neutral axis in metres.</param>
    /// <param name="secondMoment">Second moment of area I in m⁴.</param>
    /// <returns>Bending stress σ in Pa.</returns>
    public static double BendingStress(this double bendingMoment, double distanceFromNeutralAxis, double secondMoment)
    {
        if (secondMoment <= 0) throw new ArgumentException("Second moment of area must be greater than zero.");
        return bendingMoment * distanceFromNeutralAxis / secondMoment;
    }

    /// <summary>
    /// Beam shear force from the third derivative: V = E · I · d³u/dx³.
    /// </summary>
    /// <param name="youngsModulus">Young's modulus E in Pa.</param>
    /// <param name="secondMoment">Second moment of area I in m⁴.</param>
    /// <param name="thirdDerivative">Third derivative d³u/dx³ in 1/m².</param>
    /// <returns>Shear force V in newtons.</returns>
    public static double BeamShearForce(this double youngsModulus, double secondMoment, double thirdDerivative)
    {
        return youngsModulus * secondMoment * thirdDerivative;
    }

    /// <summary>
    /// Distributed load intensity from the Euler–Bernoulli beam equation: q = E · I · d⁴u/dx⁴.
    /// This is the core relationship EIu⁗ = q.
    /// </summary>
    /// <param name="youngsModulus">Young's modulus E in Pa.</param>
    /// <param name="secondMoment">Second moment of area I in m⁴.</param>
    /// <param name="fourthDerivative">Fourth derivative d⁴u/dx⁴ in 1/m³.</param>
    /// <returns>Load intensity q in N/m.</returns>
    public static double BeamLoadIntensity(this double youngsModulus, double secondMoment, double fourthDerivative)
    {
        return youngsModulus * secondMoment * fourthDerivative;
    }

    /// <summary>
    /// Computes the Euler–Bernoulli residual EI·d⁴u/dx⁴ − q at each node
    /// using finite differences (<see cref="GridOperators.Biharmonic1D"/>).
    /// A zero residual indicates the beam equation is satisfied.
    /// </summary>
    /// <param name="u">Nodal deflections along the beam.</param>
    /// <param name="EI">Flexural rigidity E·I in N·m².</param>
    /// <param name="dx">Nodal spacing in metres.</param>
    /// <param name="q">Distributed load at each node in N/m.</param>
    /// <param name="bc">Boundary condition type.</param>
    /// <returns>Residual vector (EI·d⁴u/dx⁴ − q) at each node.</returns>
    public static VectorN EulerBernoulliResidual(this VectorN u, double EI, double dx, VectorN q,
        BoundaryCondition bc = BoundaryCondition.Dirichlet)
    {
        var d4u = GridOperators.Biharmonic1D(u, dx, bc);
        var residual = new double[u.Length];
        for (int i = 0; i < u.Length; i++)
            residual[i] = EI * d4u[i] - q[i];
        return new VectorN(residual);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Analytical Beam Deflections
    // ═══════════════════════════════════════════════════════════════

    #region Analytical Beam Deflections

    /// <summary>
    /// Cantilever beam with point load P at the free end — deflection at position x
    /// (x measured from the fixed end): u(x) = Px²(3L − x) / (6EI).
    /// </summary>
    /// <param name="load">Point load P in newtons (positive downward).</param>
    /// <param name="length">Beam length L in metres.</param>
    /// <param name="EI">Flexural rigidity E·I in N·m².</param>
    /// <param name="x">Position from the fixed end in metres.</param>
    /// <returns>Deflection u(x) in metres (positive downward).</returns>
    public static double CantileverPointLoadDeflection(this double load, double length, double EI, double x)
    {
        if (EI <= 0) throw new ArgumentException("Flexural rigidity must be greater than zero.");
        if (x < 0 || x > length) throw new ArgumentException("Position x must be between 0 and L.");
        return load * x * x * (3.0 * length - x) / (6.0 * EI);
    }

    /// <summary>
    /// Cantilever beam with point load P at the free end — maximum deflection:
    /// δ_max = PL³ / (3EI).
    /// </summary>
    /// <param name="load">Point load P in newtons.</param>
    /// <param name="length">Beam length L in metres.</param>
    /// <param name="EI">Flexural rigidity E·I in N·m².</param>
    /// <returns>Maximum deflection at the free end in metres.</returns>
    public static double CantileverPointLoadMaxDeflection(this double load, double length, double EI)
    {
        if (EI <= 0) throw new ArgumentException("Flexural rigidity must be greater than zero.");
        return load * length * length * length / (3.0 * EI);
    }

    /// <summary>
    /// Cantilever beam with uniform distributed load q — deflection at position x:
    /// u(x) = q·x²(6L² − 4Lx + x²) / (24EI).
    /// </summary>
    /// <param name="loadPerLength">Distributed load q in N/m.</param>
    /// <param name="length">Beam length L in metres.</param>
    /// <param name="EI">Flexural rigidity E·I in N·m².</param>
    /// <param name="x">Position from the fixed end in metres.</param>
    /// <returns>Deflection u(x) in metres.</returns>
    public static double CantileverUniformLoadDeflection(this double loadPerLength, double length, double EI, double x)
    {
        if (EI <= 0) throw new ArgumentException("Flexural rigidity must be greater than zero.");
        if (x < 0 || x > length) throw new ArgumentException("Position x must be between 0 and L.");
        double L = length;
        return loadPerLength * x * x * (6.0 * L * L - 4.0 * L * x + x * x) / (24.0 * EI);
    }

    /// <summary>
    /// Cantilever beam with uniform distributed load q — maximum deflection:
    /// δ_max = qL⁴ / (8EI).
    /// </summary>
    /// <param name="loadPerLength">Distributed load q in N/m.</param>
    /// <param name="length">Beam length L in metres.</param>
    /// <param name="EI">Flexural rigidity E·I in N·m².</param>
    /// <returns>Maximum deflection at the free end in metres.</returns>
    public static double CantileverUniformLoadMaxDeflection(this double loadPerLength, double length, double EI)
    {
        if (EI <= 0) throw new ArgumentException("Flexural rigidity must be greater than zero.");
        return loadPerLength * length * length * length * length / (8.0 * EI);
    }

    /// <summary>
    /// Simply supported beam with uniform load q — deflection at position x:
    /// u(x) = qx(L³ − 2Lx² + x³) / (24EI).
    /// </summary>
    /// <param name="loadPerLength">Distributed load q in N/m.</param>
    /// <param name="length">Beam span L in metres.</param>
    /// <param name="EI">Flexural rigidity E·I in N·m².</param>
    /// <param name="x">Position from the left support in metres.</param>
    /// <returns>Deflection u(x) in metres.</returns>
    public static double SimplySupportedUniformLoadDeflection(this double loadPerLength, double length, double EI, double x)
    {
        if (EI <= 0) throw new ArgumentException("Flexural rigidity must be greater than zero.");
        if (x < 0 || x > length) throw new ArgumentException("Position x must be between 0 and L.");
        double L = length;
        return loadPerLength * x * (L * L * L - 2.0 * L * x * x + x * x * x) / (24.0 * EI);
    }

    /// <summary>
    /// Simply supported beam with uniform load q — maximum deflection at midspan:
    /// δ_max = 5qL⁴ / (384EI).
    /// </summary>
    /// <param name="loadPerLength">Distributed load q in N/m.</param>
    /// <param name="length">Beam span L in metres.</param>
    /// <param name="EI">Flexural rigidity E·I in N·m².</param>
    /// <returns>Maximum deflection at midspan in metres.</returns>
    public static double SimplySupportedUniformLoadMaxDeflection(this double loadPerLength, double length, double EI)
    {
        if (EI <= 0) throw new ArgumentException("Flexural rigidity must be greater than zero.");
        return 5.0 * loadPerLength * length * length * length * length / (384.0 * EI);
    }

    /// <summary>
    /// Simply supported beam with point load P at midspan — deflection at position x
    /// (for x ≤ L/2, symmetric): u(x) = Px(3L² − 4x²) / (48EI).
    /// </summary>
    /// <param name="load">Point load P in newtons at midspan.</param>
    /// <param name="length">Beam span L in metres.</param>
    /// <param name="EI">Flexural rigidity E·I in N·m².</param>
    /// <param name="x">Position from the left support in metres.</param>
    /// <returns>Deflection u(x) in metres.</returns>
    public static double SimplySupportedPointLoadDeflection(this double load, double length, double EI, double x)
    {
        if (EI <= 0) throw new ArgumentException("Flexural rigidity must be greater than zero.");
        if (x < 0 || x > length) throw new ArgumentException("Position x must be between 0 and L.");
        // Exploit symmetry: map to first half
        double xm = x <= length / 2.0 ? x : length - x;
        return load * xm * (3.0 * length * length - 4.0 * xm * xm) / (48.0 * EI);
    }

    /// <summary>
    /// Simply supported beam with point load P at midspan — maximum deflection:
    /// δ_max = PL³ / (48EI).
    /// </summary>
    /// <param name="load">Point load P in newtons.</param>
    /// <param name="length">Beam span L in metres.</param>
    /// <param name="EI">Flexural rigidity E·I in N·m².</param>
    /// <returns>Maximum deflection at midspan in metres.</returns>
    public static double SimplySupportedPointLoadMaxDeflection(this double load, double length, double EI)
    {
        if (EI <= 0) throw new ArgumentException("Flexural rigidity must be greater than zero.");
        return load * length * length * length / (48.0 * EI);
    }

    #endregion
}
