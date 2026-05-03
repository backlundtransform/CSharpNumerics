namespace CSharpNumerics.Numerics.FiniteElement;

using System;
using CSharpNumerics.Numerics.FiniteElement.Interfaces;
using CSharpNumerics.Numerics.Objects;

/// <summary>
/// 3-node Constant Strain Triangle (CST) element for 2D plane stress/strain analysis.
/// Linear displacement interpolation — strain is constant within the element.
/// DOFs: (ux1, uy1, ux2, uy2, ux3, uy3) → 6×6 stiffness matrix.
/// </summary>
public class TriElement : IElement2D
{
    /// <inheritdoc/>
    public int NodesPerElement => 3;

    /// <inheritdoc/>
    public int DofsPerNode => 2;

    /// <inheritdoc/>
    public int TotalDofs => 6;

    /// <inheritdoc/>
    public Matrix LocalStiffness(double[,] nodeCoords, double thickness, double E, double nu, bool planeStress)
    {
        double x1 = nodeCoords[0, 0], y1 = nodeCoords[0, 1];
        double x2 = nodeCoords[1, 0], y2 = nodeCoords[1, 1];
        double x3 = nodeCoords[2, 0], y3 = nodeCoords[2, 1];

        // Element area (signed)
        double area = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
        if (area <= 0)
            throw new InvalidOperationException("Element has non-positive area. Check node ordering (counter-clockwise required).");

        double twoA = 2.0 * area;

        // Shape function derivatives (constant for CST):
        // dN1/dx = (y2-y3)/(2A), dN2/dx = (y3-y1)/(2A), dN3/dx = (y1-y2)/(2A)
        // dN1/dy = (x3-x2)/(2A), dN2/dy = (x1-x3)/(2A), dN3/dy = (x2-x1)/(2A)
        double dN1dx = (y2 - y3) / twoA;
        double dN2dx = (y3 - y1) / twoA;
        double dN3dx = (y1 - y2) / twoA;
        double dN1dy = (x3 - x2) / twoA;
        double dN2dy = (x1 - x3) / twoA;
        double dN3dy = (x2 - x1) / twoA;

        // B-matrix (3×6): strain = B · u
        // [εxx]   [dN1/dx  0      dN2/dx  0      dN3/dx  0    ] [u1]
        // [εyy] = [0       dN1/dy 0       dN2/dy 0       dN3/dy] [v1]
        // [γxy]   [dN1/dy  dN1/dx dN2/dy  dN2/dx dN3/dy  dN3/dx] [...]
        var B = new double[3, 6];
        B[0, 0] = dN1dx; B[0, 2] = dN2dx; B[0, 4] = dN3dx;
        B[1, 1] = dN1dy; B[1, 3] = dN2dy; B[1, 5] = dN3dy;
        B[2, 0] = dN1dy; B[2, 1] = dN1dx;
        B[2, 2] = dN2dy; B[2, 3] = dN2dx;
        B[2, 4] = dN3dy; B[2, 5] = dN3dx;

        // Constitutive matrix D (3×3)
        var D = ConstitutiveMatrix(E, nu, planeStress);

        // Ke = t · A · Bᵀ · D · B (closed-form for CST — no numerical integration)
        // First compute DB = D · B (3×6)
        var DB = new double[3, 6];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 6; j++)
            {
                double sum = 0;
                for (int k = 0; k < 3; k++)
                    sum += D[i, k] * B[k, j];
                DB[i, j] = sum;
            }

        // Ke = t · A · Bᵀ · DB (6×6)
        double factor = thickness * area;
        var ke = new double[6, 6];
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
            {
                double sum = 0;
                for (int k = 0; k < 3; k++)
                    sum += B[k, i] * DB[k, j];
                ke[i, j] = factor * sum;
            }

        return new Matrix(ke);
    }

    /// <inheritdoc/>
    public VectorN LocalLoad(double[,] nodeCoords, double thickness, double qx, double qy)
    {
        double x1 = nodeCoords[0, 0], y1 = nodeCoords[0, 1];
        double x2 = nodeCoords[1, 0], y2 = nodeCoords[1, 1];
        double x3 = nodeCoords[2, 0], y3 = nodeCoords[2, 1];

        double area = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
        double factor = thickness * area / 3.0;

        // Consistent load: distribute equally to all 3 nodes
        return new VectorN(new[]
        {
            factor * qx, factor * qy,
            factor * qx, factor * qy,
            factor * qx, factor * qy
        });
    }

    /// <inheritdoc/>
    public double[] ShapeFunctions(double xi, double eta)
    {
        // Area coordinates: N1 = 1 - xi - eta, N2 = xi, N3 = eta
        return new[] { 1.0 - xi - eta, xi, eta };
    }

    /// <summary>
    /// Computes the constitutive (material) matrix D for plane stress or plane strain.
    /// </summary>
    internal static double[,] ConstitutiveMatrix(double E, double nu, bool planeStress)
    {
        var D = new double[3, 3];

        if (planeStress)
        {
            double c = E / (1.0 - nu * nu);
            D[0, 0] = c;
            D[0, 1] = c * nu;
            D[1, 0] = c * nu;
            D[1, 1] = c;
            D[2, 2] = c * (1.0 - nu) / 2.0;
        }
        else
        {
            // Plane strain
            double c = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
            D[0, 0] = c * (1.0 - nu);
            D[0, 1] = c * nu;
            D[1, 0] = c * nu;
            D[1, 1] = c * (1.0 - nu);
            D[2, 2] = c * (1.0 - 2.0 * nu) / 2.0;
        }

        return D;
    }
}
