namespace CSharpNumerics.Numerics.FiniteElement;

using System;
using System.Collections.Generic;
using CSharpNumerics.Numerics.FiniteElement.Interfaces;
using CSharpNumerics.Numerics.Objects;

/// <summary>
/// Assembles a global stiffness matrix and load vector from 1D finite elements,
/// applies Dirichlet boundary conditions, and solves the resulting linear system.
/// </summary>
public class Assembler1D
{
    private readonly IElement1D[] _elements;
    private readonly Mesh1D _mesh;
    private readonly int _dofsPerNode;
    private readonly int _totalDofs;

    /// <summary>Global stiffness matrix after assembly.</summary>
    public Matrix GlobalStiffness;

    /// <summary>Global load vector after assembly.</summary>
    public VectorN GlobalLoad;

    /// <param name="mesh">The 1D mesh.</param>
    /// <param name="elements">Array of elements matching the mesh connectivity.</param>
    public Assembler1D(Mesh1D mesh, IElement1D[] elements)
    {
        _mesh = mesh;
        _elements = elements;
        _dofsPerNode = elements[0].DofsPerNode;
        _totalDofs = mesh.NodeCount * _dofsPerNode;
        GlobalStiffness = new Matrix(_totalDofs, _totalDofs);
        GlobalLoad = new VectorN(_totalDofs);
    }

    /// <summary>
    /// Assembles all element stiffness matrices and load vectors into the global system.
    /// </summary>
    /// <param name="distributedLoad">Uniform distributed load per unit length (applied to all elements).</param>
    public void Assemble(double distributedLoad = 0.0)
    {
        GlobalStiffness = new Matrix(_totalDofs, _totalDofs);
        GlobalLoad = new VectorN(_totalDofs);

        for (int e = 0; e < _elements.Length; e++)
        {
            var elem = _elements[e];
            var ke = elem.LocalStiffness();
            var fe = elem.LocalLoad(distributedLoad);

            var (nodeA, nodeB) = _mesh.Elements[e];
            int[] globalDofs = GetElementDofs(nodeA, nodeB);

            for (int i = 0; i < elem.TotalDofs; i++)
            {
                GlobalLoad[globalDofs[i]] += fe[i];

                for (int j = 0; j < elem.TotalDofs; j++)
                    GlobalStiffness.values[globalDofs[i], globalDofs[j]] += ke.values[i, j];
            }
        }
    }

    /// <summary>
    /// Applies a concentrated nodal load at a specific DOF.
    /// </summary>
    /// <param name="nodeIndex">Index of the node.</param>
    /// <param name="localDof">Local DOF index at the node (0 for displacement, 1 for rotation in beam elements).</param>
    /// <param name="value">Load value.</param>
    public void ApplyNodalLoad(int nodeIndex, int localDof, double value)
    {
        int globalDof = nodeIndex * _dofsPerNode + localDof;
        GlobalLoad[globalDof] += value;
    }

    /// <summary>
    /// Solves the global system Ku = F after applying Dirichlet boundary conditions
    /// by row/column elimination and Gaussian elimination with partial pivoting.
    /// </summary>
    /// <param name="fixedDofs">Dictionary mapping global DOF index → prescribed value.</param>
    /// <returns>Solution vector of all DOFs.</returns>
    public VectorN Solve(Dictionary<int, double> fixedDofs)
    {
        int n = _totalDofs;
        var A = new double[n, n];
        var b = new double[n];

        for (int i = 0; i < n; i++)
        {
            b[i] = GlobalLoad.Values[i];
            for (int j = 0; j < n; j++)
                A[i, j] = GlobalStiffness.values[i, j];
        }

        // Apply Dirichlet BCs by row/column elimination
        foreach (var (dof, val) in fixedDofs)
        {
            for (int i = 0; i < n; i++)
            {
                b[i] -= A[i, dof] * val;
                A[i, dof] = 0;
                A[dof, i] = 0;
            }
            A[dof, dof] = 1.0;
            b[dof] = val;
        }

        // Gaussian elimination with partial pivoting
        for (int col = 0; col < n; col++)
        {
            int maxRow = col;
            double maxVal = Math.Abs(A[col, col]);
            for (int row = col + 1; row < n; row++)
            {
                double v = Math.Abs(A[row, col]);
                if (v > maxVal) { maxVal = v; maxRow = row; }
            }

            if (maxRow != col)
            {
                for (int j = col; j < n; j++)
                    (A[col, j], A[maxRow, j]) = (A[maxRow, j], A[col, j]);
                (b[col], b[maxRow]) = (b[maxRow], b[col]);
            }

            double pivot = A[col, col];
            for (int row = col + 1; row < n; row++)
            {
                double factor = A[row, col] / pivot;
                for (int j = col; j < n; j++)
                    A[row, j] -= factor * A[col, j];
                b[row] -= factor * b[col];
            }
        }

        // Back substitution
        var x = new double[n];
        for (int i = n - 1; i >= 0; i--)
        {
            double sum = b[i];
            for (int j = i + 1; j < n; j++)
                sum -= A[i, j] * x[j];
            x[i] = sum / A[i, i];
        }

        return new VectorN(x);
    }

    private int[] GetElementDofs(int nodeA, int nodeB)
    {
        var dofs = new int[_dofsPerNode * 2];
        for (int d = 0; d < _dofsPerNode; d++)
        {
            dofs[d] = nodeA * _dofsPerNode + d;
            dofs[_dofsPerNode + d] = nodeB * _dofsPerNode + d;
        }
        return dofs;
    }
}
