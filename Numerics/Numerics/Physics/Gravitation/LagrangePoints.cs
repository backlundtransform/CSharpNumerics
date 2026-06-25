using System;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.Gravitation;

/// <summary>
/// The five Lagrange points of the circular restricted three-body problem for
/// two primaries of mass m₁ (larger) and m₂ (smaller) separated by a fixed
/// distance. Points are returned in the co-rotating frame with the larger mass
/// at the origin and the smaller mass on the +x axis at (separation, 0, 0).
///
/// L1–L3 are the collinear points (solved numerically); L4 and L5 are the
/// triangular points that form equilateral triangles with the two primaries.
/// </summary>
public static class LagrangePoints
{
    /// <summary>
    /// Computes all five Lagrange points (L1–L5) for the given primaries and
    /// separation. m1 sits at the origin, m2 at (separation, 0, 0).
    /// </summary>
    /// <param name="largeMass">Mass of the larger primary m₁ (kg).</param>
    /// <param name="smallMass">Mass of the smaller primary m₂ (kg), ≤ m₁.</param>
    /// <param name="separation">Distance between the primaries (m).</param>
    public static (Vector L1, Vector L2, Vector L3, Vector L4, Vector L5) All(
        double largeMass, double smallMass, double separation)
    {
        if (largeMass <= 0 || smallMass <= 0)
            throw new ArgumentOutOfRangeException(nameof(largeMass), "Masses must be positive.");
        if (smallMass > largeMass)
            throw new ArgumentException("smallMass must not exceed largeMass.", nameof(smallMass));
        if (separation <= 0)
            throw new ArgumentOutOfRangeException(nameof(separation), "Separation must be positive.");

        double mu = smallMass / (largeMass + smallMass); // mass parameter
        double cubeRoot = Math.Pow(mu / 3.0, 1.0 / 3.0);

        // Solve the collinear equilibrium equation in the barycentric, normalised
        // frame (R = 1): m1 at x = −mu, m2 at x = 1 − mu.
        double xL1 = SolveCollinear((1.0 - mu) - cubeRoot, mu);
        double xL2 = SolveCollinear((1.0 - mu) + cubeRoot, mu);
        double xL3 = SolveCollinear(-(1.0 + 5.0 * mu / 12.0), mu);

        // Shift (+mu) so m1 is at the origin, then scale by the separation.
        Vector OnAxis(double xNorm) => new Vector((xNorm + mu) * separation, 0, 0);

        // Triangular points: equilateral triangle with both primaries.
        double halfX = separation / 2.0;
        double height = separation * Math.Sqrt(3.0) / 2.0;

        return (
            OnAxis(xL1),
            OnAxis(xL2),
            OnAxis(xL3),
            new Vector(halfX, height, 0),
            new Vector(halfX, -height, 0));
    }

    // Newton's method on f(x) = x − (1−mu)(x+mu)/|x+mu|³ − mu(x−1+mu)/|x−1+mu|³.
    private static double SolveCollinear(double initialGuess, double mu)
    {
        double x = initialGuess;
        double m1Pos = -mu;             // larger mass position
        double m2Pos = 1.0 - mu;        // smaller mass position

        for (int iter = 0; iter < 100; iter++)
        {
            double d1 = x - m1Pos;
            double d2 = x - m2Pos;
            double a1 = Math.Abs(d1);
            double a2 = Math.Abs(d2);

            double f = x
                       - (1.0 - mu) * d1 / (a1 * a1 * a1)
                       - mu * d2 / (a2 * a2 * a2);

            // f'(x) = 1 + 2(1−mu)/|d1|³ + 2·mu/|d2|³  (always positive)
            double fp = 1.0
                        + 2.0 * (1.0 - mu) / (a1 * a1 * a1)
                        + 2.0 * mu / (a2 * a2 * a2);

            double dx = f / fp;
            x -= dx;
            if (Math.Abs(dx) < 1e-14) break;
        }

        return x;
    }
}
