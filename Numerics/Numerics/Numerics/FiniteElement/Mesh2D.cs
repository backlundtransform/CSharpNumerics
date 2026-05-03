namespace CSharpNumerics.Numerics.FiniteElement;

using System;
using CSharpNumerics.Numerics.FiniteElement.Enums;

/// <summary>
/// A structured 2D mesh over a rectangular domain, supporting triangular or quadrilateral elements.
/// Nodes are arranged in a regular grid with (nx+1) × (ny+1) nodes.
/// </summary>
public class Mesh2D
{
    /// <summary>Node coordinates [nodeCount × 2]. Each row is (x, y).</summary>
    public double[,] Nodes { get; }

    /// <summary>Element connectivity [elementCount × nodesPerElement]. Each row lists node indices.</summary>
    public int[,] Elements { get; }

    /// <summary>Total number of nodes.</summary>
    public int NodeCount { get; }

    /// <summary>Total number of elements.</summary>
    public int ElementCount { get; }

    /// <summary>Number of nodes per element (3 for Tri, 4 for Quad).</summary>
    public int NodesPerElement { get; }

    /// <summary>Element type used in this mesh.</summary>
    public ElementType Type { get; }

    /// <summary>Number of elements in the x-direction.</summary>
    public int Nx { get; }

    /// <summary>Number of elements in the y-direction.</summary>
    public int Ny { get; }

    /// <summary>Domain width.</summary>
    public double Width { get; }

    /// <summary>Domain height.</summary>
    public double Height { get; }

    /// <summary>
    /// Creates a structured 2D mesh over a rectangular domain [0, width] × [0, height].
    /// </summary>
    /// <param name="width">Domain width.</param>
    /// <param name="height">Domain height.</param>
    /// <param name="nx">Number of element divisions in x.</param>
    /// <param name="ny">Number of element divisions in y.</param>
    /// <param name="type">Element type (Tri or Quad).</param>
    public Mesh2D(double width, double height, int nx, int ny, ElementType type)
    {
        if (width <= 0 || height <= 0)
            throw new ArgumentException("Width and height must be positive.");
        if (nx <= 0 || ny <= 0)
            throw new ArgumentException("nx and ny must be positive.");

        Width = width;
        Height = height;
        Nx = nx;
        Ny = ny;
        Type = type;

        int nodesX = nx + 1;
        int nodesY = ny + 1;
        NodeCount = nodesX * nodesY;

        double dx = width / nx;
        double dy = height / ny;

        // Generate nodes: row-major (x varies fastest)
        Nodes = new double[NodeCount, 2];
        for (int iy = 0; iy < nodesY; iy++)
        {
            for (int ix = 0; ix < nodesX; ix++)
            {
                int nodeIdx = iy * nodesX + ix;
                Nodes[nodeIdx, 0] = ix * dx;
                Nodes[nodeIdx, 1] = iy * dy;
            }
        }

        // Generate element connectivity
        if (type == ElementType.Quad)
        {
            NodesPerElement = 4;
            ElementCount = nx * ny;
            Elements = new int[ElementCount, 4];

            for (int ey = 0; ey < ny; ey++)
            {
                for (int ex = 0; ex < nx; ex++)
                {
                    int elemIdx = ey * nx + ex;
                    int n0 = ey * nodesX + ex;         // bottom-left
                    int n1 = n0 + 1;                    // bottom-right
                    int n2 = n0 + nodesX + 1;           // top-right
                    int n3 = n0 + nodesX;               // top-left

                    Elements[elemIdx, 0] = n0;
                    Elements[elemIdx, 1] = n1;
                    Elements[elemIdx, 2] = n2;
                    Elements[elemIdx, 3] = n3;
                }
            }
        }
        else // Tri
        {
            NodesPerElement = 3;
            ElementCount = 2 * nx * ny; // each rectangle → 2 triangles
            Elements = new int[ElementCount, 3];

            int idx = 0;
            for (int ey = 0; ey < ny; ey++)
            {
                for (int ex = 0; ex < nx; ex++)
                {
                    int n0 = ey * nodesX + ex;         // bottom-left
                    int n1 = n0 + 1;                    // bottom-right
                    int n2 = n0 + nodesX + 1;           // top-right
                    int n3 = n0 + nodesX;               // top-left

                    // Triangle 1: n0, n1, n2 (lower-right)
                    Elements[idx, 0] = n0;
                    Elements[idx, 1] = n1;
                    Elements[idx, 2] = n2;
                    idx++;

                    // Triangle 2: n0, n2, n3 (upper-left)
                    Elements[idx, 0] = n0;
                    Elements[idx, 1] = n2;
                    Elements[idx, 2] = n3;
                    idx++;
                }
            }
        }
    }

    /// <summary>
    /// Extracts the node coordinates for a given element.
    /// </summary>
    /// <param name="elemIndex">Element index.</param>
    /// <returns>Array [nodesPerElement × 2] of node coordinates.</returns>
    public double[,] GetElementNodes(int elemIndex)
    {
        var coords = new double[NodesPerElement, 2];
        for (int i = 0; i < NodesPerElement; i++)
        {
            int nodeIdx = Elements[elemIndex, i];
            coords[i, 0] = Nodes[nodeIdx, 0];
            coords[i, 1] = Nodes[nodeIdx, 1];
        }
        return coords;
    }

    /// <summary>
    /// Gets the node index for a grid position (ix, iy) where ix ∈ [0, nx] and iy ∈ [0, ny].
    /// </summary>
    public int GetNodeIndex(int ix, int iy)
    {
        return iy * (Nx + 1) + ix;
    }
}
