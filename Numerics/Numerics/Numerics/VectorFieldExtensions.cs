using Numerics.Objects;
using System.Collections.Generic;

namespace System;

/// <summary>
/// Provides extension methods for vector calculus operations on scalar fields:
/// gradient and Laplacian.
/// </summary>
public static class VectorFieldExtensions
{
    /// <summary>
    /// Creates a gradient vector field over a diagonal line of points starting at (xmin, ymin, zmin).
    /// Each step evaluates ∇f at (xmin+i, ymin+i, zmin+i).
    /// </summary>
    /// <param name="func">The scalar field.</param>
    /// <param name="xmin">Starting x coordinate.</param>
    /// <param name="ymin">Starting y coordinate.</param>
    /// <param name="zmin">Starting z coordinate.</param>
    /// <param name="stepSize">Increment per step.</param>
    /// <param name="maxSteps">Maximum value of the step parameter.</param>
    /// <returns>A dictionary mapping points to their gradient vectors.</returns>
    public static IDictionary<Vector, Vector> Gradient(this Func<Vector, double> func, double xmin, double ymin, double zmin, double stepSize, double maxSteps)
    {
        var vectorField = new Dictionary<Vector, Vector>();
        for (var i = 0.0; i < maxSteps; i += stepSize)
        {
            vectorField.Add(new Vector((xmin + i, ymin + i, zmin + i)), func.Gradient((xmin + i, ymin + i, zmin + i)));
        }
        return vectorField;
    }

    /// <summary>
    /// Computes the gradient ∇f at a point: (∂f/∂x, ∂f/∂y, ∂f/∂z).
    /// </summary>
    /// <param name="func">The scalar field.</param>
    /// <param name="points">Point at which to compute the gradient.</param>
    /// <returns>The gradient vector at <paramref name="points"/>.</returns>
    public static Vector Gradient(this Func<Vector, double> func, (double, double, double) points)
    {
        var dx = func.Derivate(new Vector(points), CSharpNumerics.Numerics.Enums.Cartesian.x);
        var dy = func.Derivate(new Vector(points), CSharpNumerics.Numerics.Enums.Cartesian.y);
        var dz = func.Derivate(new Vector(points), CSharpNumerics.Numerics.Enums.Cartesian.z);
        return new Vector(dx, dy, dz);
    }

    /// <summary>
    /// Computes the Laplacian ∇²f at a point: ∂²f/∂x² + ∂²f/∂y² + ∂²f/∂z².
    /// </summary>
    /// <param name="func">The scalar field.</param>
    /// <param name="points">Point at which to compute the Laplacian.</param>
    /// <returns>The Laplacian at <paramref name="points"/>.</returns>
    public static double Laplacian(this Func<Vector, double> func, (double, double, double) points)
    {
        var dx2 = func.Derivate(new Vector(points), CSharpNumerics.Numerics.Enums.Cartesian.x, 2);
        var dy2 = func.Derivate(new Vector(points), CSharpNumerics.Numerics.Enums.Cartesian.y, 2);
        var dz2 = func.Derivate(new Vector(points), CSharpNumerics.Numerics.Enums.Cartesian.z, 2);
        return dx2 + dy2 + dz2;
    }
}
