namespace CSharpNumerics.Numerics.FiniteElement;

using System;
using CSharpNumerics.Numerics.FiniteElement.Interfaces;
using CSharpNumerics.Numerics.Objects;

/// <summary>
/// 4-node bilinear quadrilateral (Q4) element for 2D plane stress/strain analysis.
/// DOFs: (ux1, uy1, ux2, uy2, ux3, uy3, ux4, uy4).
/// Stiffness is integrated with 2×2 Gauss quadrature.
/// </summary>
public class QuadElement : IElement2D
{
    /// <inheritdoc/>
    public int NodesPerElement => 4;

    /// <inheritdoc/>
    public int DofsPerNode => 2;

    /// <inheritdoc/>
    public int TotalDofs => 8;

    /// <inheritdoc/>
    public Matrix LocalStiffness(double[,] nodeCoords, double thickness, double E, double nu, bool planeStress)
    {
        var d = TriElement.ConstitutiveMatrix(E, nu, planeStress);
        var stiffness = new double[8, 8];

        double invSqrt3 = 1.0 / Math.Sqrt(3.0);
        double[] gaussPoints = { -invSqrt3, invSqrt3 };

        foreach (double xi in gaussPoints)
        {
            foreach (double eta in gaussPoints)
            {
                ComputeKinematics(nodeCoords, xi, eta, out var b, out double detJ);

                var db = new double[3, 8];
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 8; j++)
                    {
                        double sum = 0.0;
                        for (int k = 0; k < 3; k++)
                            sum += d[i, k] * b[k, j];
                        db[i, j] = sum;
                    }

                for (int i = 0; i < 8; i++)
                    for (int j = 0; j < 8; j++)
                    {
                        double sum = 0.0;
                        for (int k = 0; k < 3; k++)
                            sum += b[k, i] * db[k, j];
                        stiffness[i, j] += thickness * detJ * sum;
                    }
            }
        }

        return new Matrix(stiffness);
    }

    /// <inheritdoc/>
    public VectorN LocalLoad(double[,] nodeCoords, double thickness, double qx, double qy)
    {
        var load = new double[8];
        double invSqrt3 = 1.0 / Math.Sqrt(3.0);
        double[] gaussPoints = { -invSqrt3, invSqrt3 };

        foreach (double xi in gaussPoints)
        {
            foreach (double eta in gaussPoints)
            {
                var shape = ShapeFunctions(xi, eta);
                ComputeJacobian(nodeCoords, xi, eta, out _, out double detJ);

                for (int node = 0; node < 4; node++)
                {
                    load[node * 2] += thickness * detJ * shape[node] * qx;
                    load[node * 2 + 1] += thickness * detJ * shape[node] * qy;
                }
            }
        }

        return new VectorN(load);
    }

    /// <inheritdoc/>
    public double[] ShapeFunctions(double xi, double eta)
    {
        return
        [
            0.25 * (1.0 - xi) * (1.0 - eta),
            0.25 * (1.0 + xi) * (1.0 - eta),
            0.25 * (1.0 + xi) * (1.0 + eta),
            0.25 * (1.0 - xi) * (1.0 + eta)
        ];
    }

    private static double[] ShapeFunctionDerivativesXi(double eta)
    {
        return
        [
            -0.25 * (1.0 - eta),
             0.25 * (1.0 - eta),
             0.25 * (1.0 + eta),
            -0.25 * (1.0 + eta)
        ];
    }

    private static double[] ShapeFunctionDerivativesEta(double xi)
    {
        return
        [
            -0.25 * (1.0 - xi),
            -0.25 * (1.0 + xi),
             0.25 * (1.0 + xi),
             0.25 * (1.0 - xi)
        ];
    }

    private static void ComputeKinematics(double[,] nodeCoords, double xi, double eta, out double[,] b, out double detJ)
    {
        ComputeJacobian(nodeCoords, xi, eta, out var inverseJ, out detJ);

        var dNdXi = ShapeFunctionDerivativesXi(eta);
        var dNdEta = ShapeFunctionDerivativesEta(xi);
        var dNdx = new double[4];
        var dNdy = new double[4];

        for (int i = 0; i < 4; i++)
        {
            dNdx[i] = inverseJ[0, 0] * dNdXi[i] + inverseJ[0, 1] * dNdEta[i];
            dNdy[i] = inverseJ[1, 0] * dNdXi[i] + inverseJ[1, 1] * dNdEta[i];
        }

        b = new double[3, 8];
        for (int i = 0; i < 4; i++)
        {
            int col = i * 2;
            b[0, col] = dNdx[i];
            b[1, col + 1] = dNdy[i];
            b[2, col] = dNdy[i];
            b[2, col + 1] = dNdx[i];
        }
    }

    private static void ComputeJacobian(double[,] nodeCoords, double xi, double eta, out double[,] inverseJ, out double detJ)
    {
        var dNdXi = ShapeFunctionDerivativesXi(eta);
        var dNdEta = ShapeFunctionDerivativesEta(xi);

        double j11 = 0.0;
        double j12 = 0.0;
        double j21 = 0.0;
        double j22 = 0.0;

        for (int i = 0; i < 4; i++)
        {
            double x = nodeCoords[i, 0];
            double y = nodeCoords[i, 1];
            j11 += dNdXi[i] * x;
            j12 += dNdEta[i] * x;
            j21 += dNdXi[i] * y;
            j22 += dNdEta[i] * y;
        }

        detJ = j11 * j22 - j12 * j21;
        if (detJ <= 0.0)
            throw new InvalidOperationException("Quadrilateral element has non-positive Jacobian determinant.");

        double invDet = 1.0 / detJ;
        inverseJ = new double[2, 2];
        inverseJ[0, 0] = invDet * j22;
        inverseJ[0, 1] = -invDet * j12;
        inverseJ[1, 0] = -invDet * j21;
        inverseJ[1, 1] = invDet * j11;
    }
}