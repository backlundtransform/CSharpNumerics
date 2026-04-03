using CSharpNumerics.Numerics.FiniteElement;
using CSharpNumerics.Numerics.FiniteElement.Interfaces;

namespace NumericTest;

[TestClass]
public class FiniteElementTests
{
    #region BarElement Tests

    [TestMethod]
    public void BarElement_Stiffness_IsSymmetric()
    {
        double EA = 200e9 * 0.01; // steel, 100 cm²
        double L = 1.0;
        var bar = new BarElement(EA, L);
        var K = bar.LocalStiffness();

        Assert.AreEqual(K.values[0, 1], K.values[1, 0], 1e-10);
    }

    [TestMethod]
    public void BarElement_Stiffness_DiagonalValues()
    {
        double EA = 1000.0;
        double L = 2.0;
        var bar = new BarElement(EA, L);
        var K = bar.LocalStiffness();

        Assert.AreEqual(500.0, K.values[0, 0], 1e-10); // EA/L
        Assert.AreEqual(500.0, K.values[1, 1], 1e-10);
        Assert.AreEqual(-500.0, K.values[0, 1], 1e-10);
    }

    [TestMethod]
    public void BarElement_ConsistentLoad_Uniform()
    {
        double EA = 1000.0;
        double L = 4.0;
        double q = 10.0;
        var bar = new BarElement(EA, L);
        var f = bar.LocalLoad(q);

        Assert.AreEqual(20.0, f[0], 1e-10); // qL/2
        Assert.AreEqual(20.0, f[1], 1e-10);
    }

    [TestMethod]
    public void BarElement_ShapeFunctions_AtEndpoints()
    {
        var bar = new BarElement(100, 2.0);
        var n0 = bar.ShapeFunctions(0.0);
        var nL = bar.ShapeFunctions(2.0);

        Assert.AreEqual(1.0, n0[0], 1e-10);
        Assert.AreEqual(0.0, n0[1], 1e-10);
        Assert.AreEqual(0.0, nL[0], 1e-10);
        Assert.AreEqual(1.0, nL[1], 1e-10);
    }

    [TestMethod]
    public void BarElement_TwoElement_AxialBar()
    {
        // Two-element bar: fixed at left end, force P at right end
        // Exact solution: u = Px / (EA)
        double EA = 1000.0;
        double L = 2.0; // total length
        int nElem = 2;
        double P = 100.0;

        var mesh = new Mesh1D(0, L, nElem);
        var elements = mesh.CreateElements((_, len) => new BarElement(EA, len));
        var asm = new Assembler1D(mesh, elements);
        asm.Assemble();
        asm.ApplyNodalLoad(nElem, 0, P); // force at right tip

        var bc = new Dictionary<int, double> { { 0, 0.0 } }; // fixed at left
        var u = asm.Solve(bc);

        // Analytical: u(x) = Px/(EA)
        for (int i = 0; i < mesh.NodeCount; i++)
        {
            double x = mesh.Nodes[i];
            double exact = P * x / EA;
            Assert.AreEqual(exact, u[i], 1e-6, $"Node {i}: expected {exact}, got {u[i]}");
        }
    }

    #endregion

    #region BeamElement Tests

    [TestMethod]
    public void BeamElement_Stiffness_IsSymmetric()
    {
        double EI = 1e6;
        double L = 1.0;
        var beam = new BeamElement(EI, L);
        var K = beam.LocalStiffness();

        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                Assert.AreEqual(K.values[i, j], K.values[j, i], 1e-6,
                    $"K[{i},{j}] != K[{j},{i}]");
    }

    [TestMethod]
    public void BeamElement_Stiffness_KnownValues()
    {
        double EI = 1.0;
        double L = 1.0;
        var beam = new BeamElement(EI, L);
        var K = beam.LocalStiffness();

        // Standard textbook values for L=1, EI=1
        Assert.AreEqual(12.0, K.values[0, 0], 1e-10);
        Assert.AreEqual(6.0, K.values[0, 1], 1e-10);
        Assert.AreEqual(-12.0, K.values[0, 2], 1e-10);
        Assert.AreEqual(6.0, K.values[0, 3], 1e-10);
        Assert.AreEqual(4.0, K.values[1, 1], 1e-10);
        Assert.AreEqual(-6.0, K.values[1, 2], 1e-10);
        Assert.AreEqual(2.0, K.values[1, 3], 1e-10);
        Assert.AreEqual(12.0, K.values[2, 2], 1e-10);
        Assert.AreEqual(-6.0, K.values[2, 3], 1e-10);
        Assert.AreEqual(4.0, K.values[3, 3], 1e-10);
    }

    [TestMethod]
    public void BeamElement_ShapeFunctions_Partition()
    {
        // Shape functions should sum to [1, 0] at any point (displacement sums to 1)
        var beam = new BeamElement(1.0, 2.0);
        double L = 2.0;

        foreach (double xi in new[] { 0.0, 0.5, 1.0, 1.5, 2.0 })
        {
            var N = beam.ShapeFunctions(xi);
            // N1 + N3 should equal 1 (displacement partition of unity)
            Assert.AreEqual(1.0, N[0] + N[2], 1e-10, $"Displacement partition at xi={xi}");
        }
    }

    [TestMethod]
    public void BeamElement_ConsistentLoad_Uniform()
    {
        double EI = 1.0;
        double L = 2.0;
        double q = 6.0;
        var beam = new BeamElement(EI, L);
        var f = beam.LocalLoad(q);

        Assert.AreEqual(6.0, f[0], 1e-10);  // qL/2
        Assert.AreEqual(2.0, f[1], 1e-10);  // qL²/12
        Assert.AreEqual(6.0, f[2], 1e-10);  // qL/2
        Assert.AreEqual(-2.0, f[3], 1e-10); // -qL²/12
    }

    [TestMethod]
    public void BeamElement_Cantilever_TipLoad()
    {
        // Cantilever beam with point load P at free end
        // Exact: w_tip = PL³ / (3EI), θ_tip = PL² / (2EI)
        double EI = 1e4;
        double L = 1.0;
        double P = 100.0;
        int nElem = 10;

        var mesh = new Mesh1D(0, L, nElem);
        var elements = mesh.CreateElements((_, len) => new BeamElement(EI, len));
        var asm = new Assembler1D(mesh, elements);
        asm.Assemble();

        // Tip load (downward) at last node, DOF 0 = transverse displacement
        asm.ApplyNodalLoad(nElem, 0, P);

        // Fixed end: w=0, θ=0 at node 0
        var bc = new Dictionary<int, double>
        {
            { 0, 0.0 }, // w1 = 0
            { 1, 0.0 }  // θ1 = 0
        };

        var u = asm.Solve(bc);

        double w_exact = P * L * L * L / (3.0 * EI);
        double theta_exact = P * L * L / (2.0 * EI);

        int tipW = nElem * 2;
        int tipTheta = nElem * 2 + 1;

        Assert.AreEqual(w_exact, u[tipW], w_exact * 0.01,
            $"Tip deflection: expected {w_exact:E4}, got {u[tipW]:E4}");
        Assert.AreEqual(theta_exact, u[tipTheta], theta_exact * 0.01,
            $"Tip rotation: expected {theta_exact:E4}, got {u[tipTheta]:E4}");
    }

    [TestMethod]
    public void BeamElement_SimplySupported_UniformLoad()
    {
        // Simply supported beam with uniform load q
        // Exact mid-span deflection: w_max = 5qL⁴ / (384 EI)
        double EI = 1e4;
        double L = 2.0;
        double q = 100.0;
        int nElem = 10;

        var mesh = new Mesh1D(0, L, nElem);
        var elements = mesh.CreateElements((_, len) => new BeamElement(EI, len));
        var asm = new Assembler1D(mesh, elements);
        asm.Assemble(q);

        // Pin supports: w=0 at both ends
        var bc = new Dictionary<int, double>
        {
            { 0, 0.0 },            // w at left end
            { nElem * 2, 0.0 }     // w at right end
        };

        var u = asm.Solve(bc);

        double w_exact = 5.0 * q * L * L * L * L / (384.0 * EI);
        int midNode = nElem / 2;
        int midDof = midNode * 2;

        Assert.AreEqual(w_exact, u[midDof], w_exact * 0.02,
            $"Mid-span deflection: expected {w_exact:E4}, got {u[midDof]:E4}");
    }

    #endregion

    #region Mesh1D Tests

    [TestMethod]
    public void Mesh1D_UniformSubdivision()
    {
        var mesh = new Mesh1D(0, 10, 5);

        Assert.AreEqual(6, mesh.NodeCount);
        Assert.AreEqual(5, mesh.ElementCount);
        Assert.AreEqual(2.0, mesh.ElementLength, 1e-10);

        Assert.AreEqual(0.0, mesh.Nodes[0], 1e-10);
        Assert.AreEqual(10.0, mesh.Nodes[5], 1e-10);
        Assert.AreEqual(4.0, mesh.Nodes[2], 1e-10);
    }

    [TestMethod]
    public void Mesh1D_ElementConnectivity()
    {
        var mesh = new Mesh1D(0, 1, 3);

        Assert.AreEqual((0, 1), mesh.Elements[0]);
        Assert.AreEqual((1, 2), mesh.Elements[1]);
        Assert.AreEqual((2, 3), mesh.Elements[2]);
    }

    [TestMethod]
    public void Mesh1D_CreateElements_Factory()
    {
        var mesh = new Mesh1D(0, 6, 3);
        var elements = mesh.CreateElements((_, len) => new BarElement(100, len));

        Assert.AreEqual(3, elements.Length);
        Assert.AreEqual(2.0, elements[0].Length, 1e-10);
    }

    #endregion

    #region Assembler1D Tests

    [TestMethod]
    public void Assembler1D_GlobalStiffness_Dimensions()
    {
        var mesh = new Mesh1D(0, 1, 4);
        var elements = mesh.CreateElements((_, len) => new BarElement(100, len));
        var asm = new Assembler1D(mesh, elements);
        asm.Assemble();

        Assert.AreEqual(5, asm.GlobalStiffness.rowLength);
        Assert.AreEqual(5, asm.GlobalStiffness.columnLength);
    }

    [TestMethod]
    public void Assembler1D_BeamGlobalStiffness_Dimensions()
    {
        var mesh = new Mesh1D(0, 1, 3);
        var elements = mesh.CreateElements((_, len) => new BeamElement(1e4, len));
        var asm = new Assembler1D(mesh, elements);
        asm.Assemble();

        // 4 nodes × 2 DOFs/node = 8
        Assert.AreEqual(8, asm.GlobalStiffness.rowLength);
        Assert.AreEqual(8, asm.GlobalStiffness.columnLength);
    }

    #endregion
}
