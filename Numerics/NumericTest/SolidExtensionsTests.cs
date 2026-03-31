using CSharpNumerics.Physics.SolidMechanics;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.FiniteDifference;

namespace NumericsTests;

[TestClass]
public class SolidExtensionsTests
{
    // ════════════════════════════════════════════
    //  Stress & Strain
    // ════════════════════════════════════════════

    [TestMethod]
    public void NormalStress_ForceAndArea_ReturnsCorrectValue()
    {
        double force = 1000.0;   // N
        double area = 0.01;      // m²
        double sigma = force.NormalStress(area);
        Assert.AreEqual(100_000.0, sigma, 1e-6); // 100 kPa
    }

    [TestMethod]
    public void NormalStrain_StressAndModulus_ReturnsCorrectValue()
    {
        double E = 200e9;        // 200 GPa (steel)
        double sigma = 100e6;    // 100 MPa
        double epsilon = sigma.NormalStrain(E);
        Assert.AreEqual(0.0005, epsilon, 1e-10);
    }

    [TestMethod]
    public void HookesLaw_RoundTrip()
    {
        double E = 200e9;
        double strain = 0.001;
        double stress = E.HookesLaw(strain);  // σ = Eε
        double strainBack = stress.NormalStrain(E);  // ε = σ/E
        Assert.AreEqual(strain, strainBack, 1e-12);
    }

    [TestMethod]
    public void ShearStress_ForceAndArea()
    {
        double V = 5000.0;
        double A = 0.02;
        Assert.AreEqual(250_000.0, V.ShearStress(A), 1e-6);
    }

    [TestMethod]
    public void ShearModulus_SteelApproximation()
    {
        double E = 200e9;
        double nu = 0.3;
        double G = E.ShearModulus(nu);
        // G ≈ 76.92 GPa
        Assert.AreEqual(200e9 / 2.6, G, 1e3);
    }

    // ════════════════════════════════════════════
    //  Second Moment of Area
    // ════════════════════════════════════════════

    [TestMethod]
    public void RectangularSecondMoment_10x10cm()
    {
        double b = 0.10;  // 10 cm
        double h = 0.10;
        double I = b.RectangularSecondMoment(h);
        // I = 0.1 * 0.1³ / 12 = 8.333...e-6 m⁴
        Assert.AreEqual(8.333333e-6, I, 1e-12);
    }

    [TestMethod]
    public void CircularSecondMoment_5cmRadius()
    {
        double r = 0.05;
        double I = r.CircularSecondMoment();
        // I = π * 0.05⁴ / 4
        double expected = Math.PI * Math.Pow(0.05, 4) / 4.0;
        Assert.AreEqual(expected, I, 1e-14);
    }

    [TestMethod]
    public void TubularSecondMoment_HollowTube()
    {
        double R = 0.05;
        double r = 0.04;
        double I = R.TubularSecondMoment(r);
        double expected = Math.PI * (Math.Pow(R, 4) - Math.Pow(r, 4)) / 4.0;
        Assert.AreEqual(expected, I, 1e-14);
    }

    // ════════════════════════════════════════════
    //  Euler-Bernoulli Beam Equation
    // ════════════════════════════════════════════

    [TestMethod]
    public void BendingMoment_EI_Times_Curvature()
    {
        double E = 200e9;
        double I = 8.33e-6;
        double kappa = 0.01;  // 1/m
        double M = E.BendingMoment(I, kappa);
        Assert.AreEqual(200e9 * 8.33e-6 * 0.01, M, 1e-3);
    }

    [TestMethod]
    public void BendingStress_FlexureFormula()
    {
        double M = 10_000.0;  // N·m
        double y = 0.05;      // 5 cm from neutral axis
        double I = 8.33e-6;   // m⁴
        double sigma = M.BendingStress(y, I);
        // σ = My/I = 10000 * 0.05 / 8.33e-6 ≈ 60.024 MPa
        Assert.AreEqual(M * y / I, sigma, 1e-3);
    }

    [TestMethod]
    public void BeamLoadIntensity_EI_Times_FourthDerivative()
    {
        double E = 200e9;
        double I = 8.33e-6;
        double d4u = 0.5;  // 1/m³
        double q = E.BeamLoadIntensity(I, d4u);
        Assert.AreEqual(E * I * d4u, q, 1e-3);
    }

    // ════════════════════════════════════════════
    //  Analytical Beam Deflections
    // ════════════════════════════════════════════

    [TestMethod]
    public void CantileverPointLoad_MaxDeflection_SteelBeam()
    {
        // Steel beam: 10x10 cm cross-section, 1 m long, 1 kN tip load
        double E = 200e9;
        double I = 0.10.RectangularSecondMoment(0.10);  // 8.333e-6 m⁴
        double EI = E * I;
        double P = 1000.0;  // 1 kN
        double L = 1.0;     // 1 m

        double dmax = P.CantileverPointLoadMaxDeflection(L, EI);
        double expected = P * L * L * L / (3.0 * EI);
        Assert.AreEqual(expected, dmax, 1e-12);
    }

    [TestMethod]
    public void CantileverPointLoad_DeflectionAtFreeEnd_EqualsMax()
    {
        double P = 500.0;
        double L = 2.0;
        double EI = 1e6;

        double atEnd = P.CantileverPointLoadDeflection(L, EI, L);
        double max = P.CantileverPointLoadMaxDeflection(L, EI);
        Assert.AreEqual(max, atEnd, 1e-10);
    }

    [TestMethod]
    public void CantileverPointLoad_ZeroAtFixedEnd()
    {
        double P = 1000.0;
        double L = 1.0;
        double EI = 1e6;
        double atFixed = P.CantileverPointLoadDeflection(L, EI, 0.0);
        Assert.AreEqual(0.0, atFixed, 1e-15);
    }

    [TestMethod]
    public void CantileverUniformLoad_MaxDeflection()
    {
        double q = 500.0;   // N/m
        double L = 2.0;
        double EI = 1e6;

        double dmax = q.CantileverUniformLoadMaxDeflection(L, EI);
        double expected = q * L * L * L * L / (8.0 * EI);
        Assert.AreEqual(expected, dmax, 1e-10);
    }

    [TestMethod]
    public void CantileverUniformLoad_AtFreeEnd_EqualsMax()
    {
        double q = 500.0;
        double L = 2.0;
        double EI = 1e6;

        double atEnd = q.CantileverUniformLoadDeflection(L, EI, L);
        double max = q.CantileverUniformLoadMaxDeflection(L, EI);
        Assert.AreEqual(max, atEnd, 1e-10);
    }

    [TestMethod]
    public void SimplySupportedUniformLoad_MaxDeflection()
    {
        double q = 1000.0;
        double L = 3.0;
        double EI = 2e6;

        double dmax = q.SimplySupportedUniformLoadMaxDeflection(L, EI);
        double expected = 5.0 * q * L * L * L * L / (384.0 * EI);
        Assert.AreEqual(expected, dmax, 1e-10);
    }

    [TestMethod]
    public void SimplySupportedUniformLoad_ZeroAtSupports()
    {
        double q = 1000.0;
        double L = 3.0;
        double EI = 2e6;

        Assert.AreEqual(0.0, q.SimplySupportedUniformLoadDeflection(L, EI, 0.0), 1e-15);
        Assert.AreEqual(0.0, q.SimplySupportedUniformLoadDeflection(L, EI, L), 1e-12);
    }

    [TestMethod]
    public void SimplySupportedPointLoad_MaxAtMidspan()
    {
        double P = 1000.0;
        double L = 2.0;
        double EI = 1e6;

        double atMid = P.SimplySupportedPointLoadDeflection(L, EI, L / 2.0);
        double max = P.SimplySupportedPointLoadMaxDeflection(L, EI);
        Assert.AreEqual(max, atMid, 1e-10);
    }

    [TestMethod]
    public void SimplySupportedPointLoad_Symmetry()
    {
        double P = 1000.0;
        double L = 4.0;
        double EI = 1e6;

        double at1 = P.SimplySupportedPointLoadDeflection(L, EI, 1.0);
        double at3 = P.SimplySupportedPointLoadDeflection(L, EI, 3.0);
        Assert.AreEqual(at1, at3, 1e-10);
    }

    // ════════════════════════════════════════════
    //  Biharmonic1D (GridOperators)
    // ════════════════════════════════════════════

    [TestMethod]
    public void Biharmonic1D_OnQuarticPolynomial_Returns24()
    {
        // u(x) = x⁴ → d⁴u/dx⁴ = 24 everywhere
        int n = 101;
        double dx = 0.01;
        var values = new double[n];
        for (int i = 0; i < n; i++)
        {
            double x = i * dx;
            values[i] = x * x * x * x;
        }
        var u = new VectorN(values);

        var d4u = GridOperators.Biharmonic1D(u, dx, BoundaryCondition.Neumann);

        // Check interior points (away from boundaries where stencil is exact)
        for (int i = 4; i < n - 4; i++)
            Assert.AreEqual(24.0, d4u[i], 0.5, $"Failed at node {i}");
    }

    // ════════════════════════════════════════════
    //  Euler-Bernoulli Residual
    // ════════════════════════════════════════════

    [TestMethod]
    public void EulerBernoulliResidual_UniformLoad_CantileverDeflection()
    {
        // Verify that the analytical cantilever deflection under uniform load
        // satisfies EIu⁗ = q (residual ≈ 0 at interior points)
        double q = 1000.0;
        double L = 1.0;
        double EI = 1e6;
        int n = 201;
        double dx = L / (n - 1);

        var uValues = new double[n];
        var qValues = new double[n];
        for (int i = 0; i < n; i++)
        {
            double x = i * dx;
            uValues[i] = q.CantileverUniformLoadDeflection(L, EI, x);
            qValues[i] = q;
        }

        var u = new VectorN(uValues);
        var qVec = new VectorN(qValues);
        var residual = u.EulerBernoulliResidual(EI, dx, qVec, BoundaryCondition.Neumann);

        // Interior points should have small residual (FD approximation)
        for (int i = 10; i < n - 10; i++)
            Assert.AreEqual(0.0, residual[i], q * 0.05,
                $"Residual too large at node {i}: {residual[i]}");
    }
}
