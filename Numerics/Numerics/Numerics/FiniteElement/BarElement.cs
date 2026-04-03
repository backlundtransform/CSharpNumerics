namespace CSharpNumerics.Numerics.FiniteElement;

using CSharpNumerics.Numerics.FiniteElement.Interfaces;
using CSharpNumerics.Numerics.Objects;

/// <summary>
/// 2-node linear bar element for axial stress/strain analysis.
/// DOFs: u1, u2 (axial displacement at each node).
/// Stiffness: (EA/L) * [[1, -1], [-1, 1]]
/// </summary>
public class BarElement : IElement1D
{
    /// <summary>Cross-sectional area times Young's modulus (EA).</summary>
    public double EA { get; }

    /// <inheritdoc/>
    public double Length { get; }

    /// <inheritdoc/>
    public int DofsPerNode => 1;

    /// <inheritdoc/>
    public int TotalDofs => 2;

    /// <param name="ea">Product of Young's modulus E and cross-sectional area A.</param>
    /// <param name="length">Element length.</param>
    public BarElement(double ea, double length)
    {
        EA = ea;
        Length = length;
    }

    /// <inheritdoc/>
    public Matrix LocalStiffness()
    {
        double k = EA / Length;
        var values = new double[2, 2];
        values[0, 0] = k;
        values[0, 1] = -k;
        values[1, 0] = -k;
        values[1, 1] = k;
        return new Matrix(values);
    }

    /// <inheritdoc/>
    public VectorN LocalLoad(double q)
    {
        double half = q * Length / 2.0;
        return new VectorN(new[] { half, half });
    }

    /// <inheritdoc/>
    public VectorN ShapeFunctions(double xi)
    {
        double n1 = 1.0 - xi / Length;
        double n2 = xi / Length;
        return new VectorN(new[] { n1, n2 });
    }
}
