namespace CSharpNumerics.Numerics.FiniteElement.Interfaces;

using CSharpNumerics.Numerics.Objects;

/// <summary>
/// Interface for a 2D finite element used in plane stress/strain analysis.
/// </summary>
public interface IElement2D
{
    /// <summary>Number of nodes in this element.</summary>
    int NodesPerElement { get; }

    /// <summary>Degrees of freedom per node (always 2: ux, uy).</summary>
    int DofsPerNode { get; }

    /// <summary>Total DOFs for the element (NodesPerElement × DofsPerNode).</summary>
    int TotalDofs { get; }

    /// <summary>
    /// Computes the local element stiffness matrix.
    /// </summary>
    /// <param name="nodeCoords">Node coordinates [nodesPerElement × 2] — each row is (x, y).</param>
    /// <param name="thickness">Element thickness.</param>
    /// <param name="E">Young's modulus.</param>
    /// <param name="nu">Poisson's ratio.</param>
    /// <param name="planeStress">True for plane stress, false for plane strain.</param>
    Matrix LocalStiffness(double[,] nodeCoords, double thickness, double E, double nu, bool planeStress);

    /// <summary>
    /// Computes the consistent nodal load vector for a uniform body force.
    /// </summary>
    /// <param name="nodeCoords">Node coordinates [nodesPerElement × 2].</param>
    /// <param name="thickness">Element thickness.</param>
    /// <param name="qx">Body force in x-direction (force per unit volume).</param>
    /// <param name="qy">Body force in y-direction (force per unit volume).</param>
    VectorN LocalLoad(double[,] nodeCoords, double thickness, double qx, double qy);

    /// <summary>
    /// Evaluates shape functions at natural coordinates (xi, eta).
    /// </summary>
    /// <param name="xi">Natural coordinate ξ.</param>
    /// <param name="eta">Natural coordinate η.</param>
    double[] ShapeFunctions(double xi, double eta);
}
