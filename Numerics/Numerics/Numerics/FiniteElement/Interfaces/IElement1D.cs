namespace CSharpNumerics.Numerics.FiniteElement.Interfaces;

using CSharpNumerics.Numerics.Objects;

/// <summary>
/// Interface for a 1D finite element with two nodes.
/// </summary>
public interface IElement1D
{
    /// <summary>Number of degrees of freedom per node.</summary>
    int DofsPerNode { get; }

    /// <summary>Total degrees of freedom for the element (2 * DofsPerNode).</summary>
    int TotalDofs { get; }

    /// <summary>Physical length of the element.</summary>
    double Length { get; }

    /// <summary>
    /// Returns the element stiffness matrix in local coordinates.
    /// </summary>
    Matrix LocalStiffness();

    /// <summary>
    /// Returns the consistent nodal load vector for a uniform distributed load.
    /// </summary>
    /// <param name="q">Load intensity per unit length.</param>
    VectorN LocalLoad(double q);

    /// <summary>
    /// Evaluates the shape functions at a local coordinate xi in [0, L].
    /// </summary>
    /// <param name="xi">Local coordinate along the element.</param>
    VectorN ShapeFunctions(double xi);
}
