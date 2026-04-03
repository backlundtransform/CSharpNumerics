namespace CSharpNumerics.Numerics.FiniteElement;

using System;
using CSharpNumerics.Numerics.FiniteElement.Interfaces;

/// <summary>
/// A uniform 1D mesh produced by subdividing an interval [start, end] into equal-length elements.
/// </summary>
public class Mesh1D
{
    /// <summary>Node positions along the 1D domain.</summary>
    public double[] Nodes { get; }

    /// <summary>Element connectivity: Elements[e] = (nodeA, nodeB).</summary>
    public (int NodeA, int NodeB)[] Elements { get; }

    /// <summary>Number of elements.</summary>
    public int ElementCount { get; }

    /// <summary>Number of nodes.</summary>
    public int NodeCount { get; }

    /// <summary>Length of each element (uniform).</summary>
    public double ElementLength { get; }

    /// <param name="start">Left endpoint of the domain.</param>
    /// <param name="end">Right endpoint of the domain.</param>
    /// <param name="numElements">Number of equal-length elements.</param>
    public Mesh1D(double start, double end, int numElements)
    {
        ElementCount = numElements;
        NodeCount = numElements + 1;
        ElementLength = (end - start) / numElements;

        Nodes = new double[NodeCount];
        for (int i = 0; i < NodeCount; i++)
            Nodes[i] = start + i * ElementLength;

        Elements = new (int, int)[ElementCount];
        for (int e = 0; e < ElementCount; e++)
            Elements[e] = (e, e + 1);
    }

    /// <summary>
    /// Creates an array of IElement1D instances from this mesh using a factory function.
    /// The factory receives (elementIndex, elementLength) and returns an element.
    /// </summary>
    public IElement1D[] CreateElements(Func<int, double, IElement1D> factory)
    {
        var elems = new IElement1D[ElementCount];
        for (int e = 0; e < ElementCount; e++)
            elems[e] = factory(e, ElementLength);
        return elems;
    }
}
