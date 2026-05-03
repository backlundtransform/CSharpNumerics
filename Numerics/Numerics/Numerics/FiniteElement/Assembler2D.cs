namespace CSharpNumerics.Numerics.FiniteElement;

using System;
using System.Collections.Generic;
using CSharpNumerics.Numerics.FiniteElement.Interfaces;
using CSharpNumerics.Numerics.Objects;

/// <summary>
/// Assembles a global stiffness matrix and load vector from 2D finite elements,
/// applies Dirichlet boundary conditions, and solves via Preconditioned Conjugate Gradient.
/// </summary>
public class Assembler2D
{
    private readonly Mesh2D _mesh;
    private readonly IElement2D _element;
    private readonly double _E;
    private readonly double _nu;
    private readonly double _thickness;
    private readonly bool _planeStress;
    private readonly int _totalDofs;

    private List<(int row, int col, double val)> _triplets;
    private double[] _load;
    private SparseMatrix _globalStiffness;

    /// <summary>
    /// Creates a 2D assembler for the given mesh and element type.
    /// </summary>
    /// <param name="mesh">The 2D mesh.</param>
    /// <param name="element">Element template implementing the 2D element interface.</param>
    /// <param name="E">Young's modulus.</param>
    /// <param name="nu">Poisson's ratio.</param>
    /// <param name="thickness">Element thickness (typically 1.0 for plane stress unit-thickness problems).</param>
    /// <param name="planeStress">True for plane stress, false for plane strain.</param>
    public Assembler2D(Mesh2D mesh, IElement2D element, double E, double nu, double thickness, bool planeStress)
    {
        _mesh = mesh;
        _element = element;
        _E = E;
        _nu = nu;
        _thickness = thickness;
        _planeStress = planeStress;
        _totalDofs = mesh.NodeCount * 2;
        _load = new double[_totalDofs];
    }

    /// <summary>Total degrees of freedom in the system (2 × nodeCount).</summary>
    public int TotalDofs => _totalDofs;

    /// <summary>
    /// Assembles all element stiffness matrices into the global system.
    /// Must be called before Solve().
    /// </summary>
    public void Assemble()
    {
        _triplets = new List<(int, int, double)>();
        _load = new double[_totalDofs];

        int nodesPerElem = _element.NodesPerElement;
        int dofsPerElem = _element.TotalDofs;

        for (int e = 0; e < _mesh.ElementCount; e++)
        {
            var nodeCoords = _mesh.GetElementNodes(e);
            var ke = _element.LocalStiffness(nodeCoords, _thickness, _E, _nu, _planeStress);

            // Map local DOFs → global DOFs
            int[] globalDofs = GetElementDofs(e);

            for (int i = 0; i < dofsPerElem; i++)
                for (int j = 0; j < dofsPerElem; j++)
                {
                    double val = ke.values[i, j];
                    if (Math.Abs(val) > 1e-30)
                        _triplets.Add((globalDofs[i], globalDofs[j], val));
                }
        }

        _globalStiffness = SparseMatrix.FromTriplets(_totalDofs, _totalDofs, _triplets);
    }

    /// <summary>
    /// Applies a concentrated load at a specific node.
    /// </summary>
    /// <param name="nodeIndex">Node index.</param>
    /// <param name="direction">0 for x-direction, 1 for y-direction.</param>
    /// <param name="value">Force value.</param>
    public void ApplyNodalLoad(int nodeIndex, int direction, double value)
    {
        int dof = nodeIndex * 2 + direction;
        _load[dof] += value;
    }

    /// <summary>
    /// Applies a uniform body force to all elements via consistent load vectors.
    /// </summary>
    /// <param name="qx">Body force in x (force per unit volume).</param>
    /// <param name="qy">Body force in y (force per unit volume).</param>
    public void ApplyBodyForce(double qx, double qy)
    {
        for (int e = 0; e < _mesh.ElementCount; e++)
        {
            var nodeCoords = _mesh.GetElementNodes(e);
            var fe = _element.LocalLoad(nodeCoords, _thickness, qx, qy);
            int[] globalDofs = GetElementDofs(e);

            for (int i = 0; i < _element.TotalDofs; i++)
                _load[globalDofs[i]] += fe[i];
        }
    }

    /// <summary>
    /// Solves the assembled system with Dirichlet boundary conditions using PCG.
    /// </summary>
    /// <param name="fixedDofs">Dictionary mapping global DOF index → prescribed displacement value.</param>
    /// <returns>Solution vector of all DOFs (interleaved: ux0, uy0, ux1, uy1, ...).</returns>
    public VectorN Solve(Dictionary<int, double> fixedDofs)
    {
        if (_globalStiffness == null)
            throw new InvalidOperationException("Call Assemble() before Solve().");

        var rhs = new VectorN(_load);
        var K = _globalStiffness.ApplyDirichlet(fixedDofs, rhs);
        return K.SolvePCG(rhs);
    }

    /// <summary>
    /// Computes stresses at the centroid of each element from the solved displacement vector.
    /// Returns arrays of σxx, σyy, τxy for each element.
    /// </summary>
    /// <param name="displacement">Solution displacement vector from Solve().</param>
    /// <param name="stressXX">Output: normal stress σxx per element.</param>
    /// <param name="stressYY">Output: normal stress σyy per element.</param>
    /// <param name="stressXY">Output: shear stress τxy per element.</param>
    public void ComputeElementStresses(VectorN displacement, out double[] stressXX, out double[] stressYY, out double[] stressXY)
    {
        int nElem = _mesh.ElementCount;
        stressXX = new double[nElem];
        stressYY = new double[nElem];
        stressXY = new double[nElem];

        var D = TriElement.ConstitutiveMatrix(_E, _nu, _planeStress);

        for (int e = 0; e < nElem; e++)
        {
            var nodeCoords = _mesh.GetElementNodes(e);
            int[] globalDofs = GetElementDofs(e);

            // Extract element displacements
            int dofsPerElem = _element.TotalDofs;
            var ue = new double[dofsPerElem];
            for (int i = 0; i < dofsPerElem; i++)
                ue[i] = displacement[globalDofs[i]];

            // Compute B-matrix for this element (CST: constant B)
            var B = ComputeBMatrix(nodeCoords);

            // strain = B · ue
            var strain = new double[3];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < dofsPerElem; j++)
                    strain[i] += B[i, j] * ue[j];

            // stress = D · strain
            stressXX[e] = D[0, 0] * strain[0] + D[0, 1] * strain[1];
            stressYY[e] = D[1, 0] * strain[0] + D[1, 1] * strain[1];
            stressXY[e] = D[2, 2] * strain[2];
        }
    }

    private double[,] ComputeBMatrix(double[,] nodeCoords)
    {
        double x1 = nodeCoords[0, 0], y1 = nodeCoords[0, 1];
        double x2 = nodeCoords[1, 0], y2 = nodeCoords[1, 1];
        double x3 = nodeCoords[2, 0], y3 = nodeCoords[2, 1];

        double twoA = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);

        double dN1dx = (y2 - y3) / twoA;
        double dN2dx = (y3 - y1) / twoA;
        double dN3dx = (y1 - y2) / twoA;
        double dN1dy = (x3 - x2) / twoA;
        double dN2dy = (x1 - x3) / twoA;
        double dN3dy = (x2 - x1) / twoA;

        var B = new double[3, 6];
        B[0, 0] = dN1dx; B[0, 2] = dN2dx; B[0, 4] = dN3dx;
        B[1, 1] = dN1dy; B[1, 3] = dN2dy; B[1, 5] = dN3dy;
        B[2, 0] = dN1dy; B[2, 1] = dN1dx;
        B[2, 2] = dN2dy; B[2, 3] = dN2dx;
        B[2, 4] = dN3dy; B[2, 5] = dN3dx;

        return B;
    }

    private int[] GetElementDofs(int elemIndex)
    {
        int nodesPerElem = _element.NodesPerElement;
        var dofs = new int[nodesPerElem * 2];
        for (int i = 0; i < nodesPerElem; i++)
        {
            int nodeIdx = _mesh.Elements[elemIndex, i];
            dofs[i * 2] = nodeIdx * 2;
            dofs[i * 2 + 1] = nodeIdx * 2 + 1;
        }
        return dofs;
    }
}
