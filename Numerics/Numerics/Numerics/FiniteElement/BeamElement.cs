namespace CSharpNumerics.Numerics.FiniteElement;

using CSharpNumerics.Numerics.FiniteElement.Interfaces;
using CSharpNumerics.Numerics.Objects;

/// <summary>
/// 2-node Hermite cubic beam element for Euler-Bernoulli bending analysis.
/// DOFs per node: (w, θ) — transverse displacement and rotation.
/// Total DOFs: 4 — (w1, θ1, w2, θ2).
/// </summary>
public class BeamElement : IElement1D
{
    /// <summary>Flexural rigidity EI (Young's modulus × second moment of area).</summary>
    public double EI { get; }

    /// <inheritdoc/>
    public double Length { get; }

    /// <inheritdoc/>
    public int DofsPerNode => 2;

    /// <inheritdoc/>
    public int TotalDofs => 4;

    /// <param name="ei">Flexural rigidity (E × I).</param>
    /// <param name="length">Element length.</param>
    public BeamElement(double ei, double length)
    {
        EI = ei;
        Length = length;
    }

    /// <summary>
    /// Returns the 4×4 Euler-Bernoulli element stiffness matrix.
    /// Standard form: (EI / L³) × symmetric 4×4 matrix.
    /// </summary>
    public Matrix LocalStiffness()
    {
        double L = Length;
        double L2 = L * L;
        double L3 = L * L * L;
        double c = EI / L3;

        var k = new double[4, 4];
        k[0, 0] = 12.0 * c;
        k[0, 1] = 6.0 * L * c;
        k[0, 2] = -12.0 * c;
        k[0, 3] = 6.0 * L * c;

        k[1, 0] = 6.0 * L * c;
        k[1, 1] = 4.0 * L2 * c;
        k[1, 2] = -6.0 * L * c;
        k[1, 3] = 2.0 * L2 * c;

        k[2, 0] = -12.0 * c;
        k[2, 1] = -6.0 * L * c;
        k[2, 2] = 12.0 * c;
        k[2, 3] = -6.0 * L * c;

        k[3, 0] = 6.0 * L * c;
        k[3, 1] = 2.0 * L2 * c;
        k[3, 2] = -6.0 * L * c;
        k[3, 3] = 4.0 * L2 * c;

        return new Matrix(k);
    }

    /// <summary>
    /// Consistent nodal load vector for a uniform distributed load q (force/length).
    /// Returns [qL/2, qL²/12, qL/2, −qL²/12].
    /// </summary>
    public VectorN LocalLoad(double q)
    {
        double L = Length;
        double qL = q * L;
        return new VectorN(new[]
        {
            qL / 2.0,
            qL * L / 12.0,
            qL / 2.0,
            -qL * L / 12.0
        });
    }

    /// <summary>
    /// Evaluates the four Hermite cubic shape functions at local coordinate xi in [0, L].
    /// N1 = 1 − 3ξ² + 2ξ³  (displacement at node 1)
    /// N2 = ξ − 2ξ² + ξ³    (rotation at node 1, scaled by L)
    /// N3 = 3ξ² − 2ξ³        (displacement at node 2)
    /// N4 = −ξ² + ξ³         (rotation at node 2, scaled by L)
    /// where ξ = xi / L (normalised coordinate).
    /// </summary>
    public VectorN ShapeFunctions(double xi)
    {
        double s = xi / Length;
        double s2 = s * s;
        double s3 = s2 * s;

        return new VectorN(new[]
        {
            1.0 - 3.0 * s2 + 2.0 * s3,
            Length * (s - 2.0 * s2 + s3),
            3.0 * s2 - 2.0 * s3,
            Length * (-s2 + s3)
        });
    }
}
