namespace CSharpNumerics.Physics.SolidMechanics.Interfaces;

/// <summary>
/// Abstraction for Euler–Bernoulli beam analysis.
/// Provides analytical deflection, bending moment, shear force,
/// and stress calculations for common support conditions.
/// </summary>
public interface IBeamModel
{
    /// <summary>
    /// Computes deflection, bending moment, and shear force for a cantilever
    /// beam with a point load P applied at distance <paramref name="a"/>
    /// from the fixed end.
    /// </summary>
    (double deflection, double moment, double shear) CantileverPointLoad(
        double P, double a, double L, double EI, double x);

    /// <summary>
    /// Computes deflection, bending moment, and shear force for a cantilever
    /// beam with a uniformly distributed load q.
    /// </summary>
    (double deflection, double moment, double shear) CantileverUniformLoad(
        double q, double L, double EI, double x);

    /// <summary>
    /// Computes deflection, bending moment, and shear force for a simply
    /// supported beam with a point load P at distance <paramref name="a"/>
    /// from the left support.
    /// </summary>
    (double deflection, double moment, double shear) SimplySupportedPointLoad(
        double P, double a, double L, double EI, double x);

    /// <summary>
    /// Computes deflection, bending moment, and shear force for a simply
    /// supported beam with a uniformly distributed load q.
    /// </summary>
    (double deflection, double moment, double shear) SimplySupportedUniformLoad(
        double q, double L, double EI, double x);

    /// <summary>
    /// Computes deflection, bending moment, and shear force for a fixed-fixed
    /// beam with a point load P at distance <paramref name="a"/>
    /// from the left support.
    /// </summary>
    (double deflection, double moment, double shear) FixedFixedPointLoad(
        double P, double a, double L, double EI, double x);

    /// <summary>
    /// Computes deflection, bending moment, and shear force for a fixed-fixed
    /// beam with a uniformly distributed load q.
    /// </summary>
    (double deflection, double moment, double shear) FixedFixedUniformLoad(
        double q, double L, double EI, double x);

    /// <summary>
    /// Computes the maximum bending stress σ = |M| · y_max / I.
    /// </summary>
    double BendingStress(double moment, double halfHeight, double secondMoment);
}
