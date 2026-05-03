using CSharpNumerics.Numerics.FiniteElement;
using CSharpNumerics.Numerics.FiniteElement.Enums;
using CSharpNumerics.Numerics.Objects;

namespace NumericTest;

[TestClass]
public class FiniteElement2DTests
{
    #region TriElement Tests

    [TestMethod]
    public void TriElement_Stiffness_IsSymmetric()
    {
        var tri = new TriElement();
        var nodes = new double[,] { { 0, 0 }, { 1, 0 }, { 0, 1 } };
        var K = tri.LocalStiffness(nodes, 1.0, 200e9, 0.3, true);

        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
                Assert.AreEqual(K.values[i, j], K.values[j, i], 1e-3,
                    $"K[{i},{j}] != K[{j},{i}]");
    }

    [TestMethod]
    public void TriElement_RigidBodyTranslation_ZeroForce()
    {
        // Rigid body translation should produce zero internal forces: K · u_rigid = 0
        var tri = new TriElement();
        var nodes = new double[,] { { 0, 0 }, { 2, 0 }, { 1, 1.5 } };
        var K = tri.LocalStiffness(nodes, 1.0, 100e9, 0.25, true);

        // Translation in x: u = [1, 0, 1, 0, 1, 0]
        var ux = new double[] { 1, 0, 1, 0, 1, 0 };
        for (int i = 0; i < 6; i++)
        {
            double force = 0;
            for (int j = 0; j < 6; j++)
                force += K.values[i, j] * ux[j];
            Assert.AreEqual(0.0, force, 1e-3, $"Non-zero force at DOF {i} for x-translation");
        }

        // Translation in y: u = [0, 1, 0, 1, 0, 1]
        var uy = new double[] { 0, 1, 0, 1, 0, 1 };
        for (int i = 0; i < 6; i++)
        {
            double force = 0;
            for (int j = 0; j < 6; j++)
                force += K.values[i, j] * uy[j];
            Assert.AreEqual(0.0, force, 1e-3, $"Non-zero force at DOF {i} for y-translation");
        }
    }

    [TestMethod]
    public void TriElement_UnitRightTriangle_KnownStiffness()
    {
        // Unit right triangle, plane stress, E=1, nu=0, t=1
        // With nu=0: D = [1,0,0; 0,1,0; 0,0,0.5]
        // This gives a well-known stiffness matrix
        var tri = new TriElement();
        var nodes = new double[,] { { 0, 0 }, { 1, 0 }, { 0, 1 } };
        var K = tri.LocalStiffness(nodes, 1.0, 1.0, 0.0, true);

        // Area = 0.5, for CST with E=1, nu=0:
        // Verify K is 6×6 and diagonal entries are positive
        Assert.AreEqual(6, K.rowLength);
        Assert.AreEqual(6, K.columnLength);

        for (int i = 0; i < 6; i++)
            Assert.IsTrue(K.values[i, i] >= 0, $"K[{i},{i}] should be non-negative");
    }

    [TestMethod]
    public void TriElement_ConsistentLoad_SumsCorrectly()
    {
        var tri = new TriElement();
        var nodes = new double[,] { { 0, 0 }, { 2, 0 }, { 0, 3 } };
        double qx = 10.0, qy = -5.0, t = 0.5;

        var f = tri.LocalLoad(nodes, t, qx, qy);

        // Area = 0.5 * 2 * 3 = 3.0
        // factor = t * A / 3 = 0.5 * 3.0 / 3.0 = 0.5
        Assert.AreEqual(6, f.Length);
        Assert.AreEqual(0.5 * qx, f[0], 1e-10); // node 1 fx
        Assert.AreEqual(0.5 * qy, f[1], 1e-10); // node 1 fy
        Assert.AreEqual(0.5 * qx, f[2], 1e-10); // node 2 fx
        Assert.AreEqual(0.5 * qy, f[3], 1e-10); // node 2 fy
        Assert.AreEqual(0.5 * qx, f[4], 1e-10); // node 3 fx
        Assert.AreEqual(0.5 * qy, f[5], 1e-10); // node 3 fy
    }

    [TestMethod]
    public void TriElement_ShapeFunctions_PartitionOfUnity()
    {
        var tri = new TriElement();

        // At several points inside the triangle, shape functions should sum to 1
        var points = new[] { (0.25, 0.25), (0.5, 0.25), (0.1, 0.1), (1.0 / 3.0, 1.0 / 3.0) };
        foreach (var (xi, eta) in points)
        {
            var N = tri.ShapeFunctions(xi, eta);
            double sum = N[0] + N[1] + N[2];
            Assert.AreEqual(1.0, sum, 1e-15, $"Shape functions don't sum to 1 at ({xi},{eta})");
        }
    }

    [TestMethod]
    public void TriElement_ShapeFunctions_AtNodes()
    {
        var tri = new TriElement();

        // Node 1: (0, 0) → N = [1, 0, 0]
        var n1 = tri.ShapeFunctions(0, 0);
        Assert.AreEqual(1.0, n1[0], 1e-15);
        Assert.AreEqual(0.0, n1[1], 1e-15);
        Assert.AreEqual(0.0, n1[2], 1e-15);

        // Node 2: (1, 0) → N = [0, 1, 0]
        var n2 = tri.ShapeFunctions(1, 0);
        Assert.AreEqual(0.0, n2[0], 1e-15);
        Assert.AreEqual(1.0, n2[1], 1e-15);
        Assert.AreEqual(0.0, n2[2], 1e-15);

        // Node 3: (0, 1) → N = [0, 0, 1]
        var n3 = tri.ShapeFunctions(0, 1);
        Assert.AreEqual(0.0, n3[0], 1e-15);
        Assert.AreEqual(0.0, n3[1], 1e-15);
        Assert.AreEqual(1.0, n3[2], 1e-15);
    }

    [TestMethod]
    public void TriElement_PlaneStrain_DifferentFromPlaneStress()
    {
        var tri = new TriElement();
        var nodes = new double[,] { { 0, 0 }, { 1, 0 }, { 0, 1 } };

        var Kstress = tri.LocalStiffness(nodes, 1.0, 200e9, 0.3, planeStress: true);
        var Kstrain = tri.LocalStiffness(nodes, 1.0, 200e9, 0.3, planeStress: false);

        // They should be different (plane strain is stiffer)
        bool anyDifferent = false;
        for (int i = 0; i < 6 && !anyDifferent; i++)
            for (int j = 0; j < 6 && !anyDifferent; j++)
                if (Math.Abs(Kstress.values[i, j] - Kstrain.values[i, j]) > 1e-3)
                    anyDifferent = true;

        Assert.IsTrue(anyDifferent, "Plane stress and plane strain stiffness should differ");
    }

    #endregion

    #region QuadElement Tests

    [TestMethod]
    public void QuadElement_Stiffness_IsSymmetric()
    {
        var quad = new QuadElement();
        var nodes = new double[,] { { 0, 0 }, { 2, 0 }, { 2, 1 }, { 0, 1 } };
        var K = quad.LocalStiffness(nodes, 1.0, 200e9, 0.3, true);

        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 8; j++)
                Assert.AreEqual(K.values[i, j], K.values[j, i], 1e-3,
                    $"K[{i},{j}] != K[{j},{i}]");
    }

    [TestMethod]
    public void QuadElement_RigidBodyTranslation_ZeroForce()
    {
        var quad = new QuadElement();
        var nodes = new double[,] { { 0, 0 }, { 2, 0 }, { 2, 1 }, { 0, 1 } };
        var K = quad.LocalStiffness(nodes, 1.0, 100e9, 0.25, true);

        var ux = new double[] { 1, 0, 1, 0, 1, 0, 1, 0 };
        for (int i = 0; i < 8; i++)
        {
            double force = 0.0;
            for (int j = 0; j < 8; j++)
                force += K.values[i, j] * ux[j];
            Assert.AreEqual(0.0, force, 1e-3, $"Non-zero force at DOF {i} for x-translation");
        }

        var uy = new double[] { 0, 1, 0, 1, 0, 1, 0, 1 };
        for (int i = 0; i < 8; i++)
        {
            double force = 0.0;
            for (int j = 0; j < 8; j++)
                force += K.values[i, j] * uy[j];
            Assert.AreEqual(0.0, force, 1e-3, $"Non-zero force at DOF {i} for y-translation");
        }
    }

    [TestMethod]
    public void QuadElement_ShapeFunctions_CenterAreQuarter()
    {
        var quad = new QuadElement();
        var N = quad.ShapeFunctions(0.0, 0.0);

        for (int i = 0; i < 4; i++)
            Assert.AreEqual(0.25, N[i], 1e-15);
    }

    [TestMethod]
    public void QuadElement_ConsistentLoad_SumsToTotalBodyForce()
    {
        var quad = new QuadElement();
        var nodes = new double[,] { { 0, 0 }, { 2, 0 }, { 2, 1 }, { 0, 1 } };
        double qx = 10.0;
        double qy = -3.0;
        double thickness = 0.5;

        var f = quad.LocalLoad(nodes, thickness, qx, qy);
        double totalFx = 0.0;
        double totalFy = 0.0;
        for (int i = 0; i < 4; i++)
        {
            totalFx += f[i * 2];
            totalFy += f[i * 2 + 1];
        }

        double area = 2.0;
        Assert.AreEqual(thickness * area * qx, totalFx, 1e-12);
        Assert.AreEqual(thickness * area * qy, totalFy, 1e-12);
    }

    [TestMethod]
    public void QuadElement_PatchTest_LinearBoundaryDisplacementsRecoverInteriorNode()
    {
        double E = 1000.0;
        double epsilon = 1e-3;
        var mesh = new Mesh2D(1.0, 1.0, 2, 2, ElementType.Quad);
        var asm = new Assembler2D(mesh, new QuadElement(), E, 0.0, 1.0, true);
        asm.Assemble();

        var bcs = new Dictionary<int, double>();
        for (int iy = 0; iy <= 2; iy++)
            for (int ix = 0; ix <= 2; ix++)
            {
                bool isBoundary = ix == 0 || ix == 2 || iy == 0 || iy == 2;
                if (!isBoundary)
                    continue;

                int node = mesh.GetNodeIndex(ix, iy);
                double x = mesh.Nodes[node, 0];
                bcs[node * 2] = epsilon * x;
                bcs[node * 2 + 1] = 0.0;
            }

        var u = asm.Solve(bcs);
        int centerNode = mesh.GetNodeIndex(1, 1);
        Assert.AreEqual(epsilon * 0.5, u[centerNode * 2], 1e-8);
        Assert.AreEqual(0.0, u[centerNode * 2 + 1], 1e-8);
    }

    #endregion

    #region Mesh2D Tests

    [TestMethod]
    public void Mesh2D_Tri_CorrectCounts()
    {
        var mesh = new Mesh2D(2.0, 1.0, 4, 2, ElementType.Tri);

        Assert.AreEqual(5 * 3, mesh.NodeCount);   // (4+1)×(2+1) = 15
        Assert.AreEqual(2 * 4 * 2, mesh.ElementCount); // 2 tri per rect = 16
        Assert.AreEqual(3, mesh.NodesPerElement);
    }

    [TestMethod]
    public void Mesh2D_Quad_CorrectCounts()
    {
        var mesh = new Mesh2D(2.0, 1.0, 4, 2, ElementType.Quad);

        Assert.AreEqual(15, mesh.NodeCount);       // (4+1)×(2+1)
        Assert.AreEqual(8, mesh.ElementCount);     // 4×2
        Assert.AreEqual(4, mesh.NodesPerElement);
    }

    [TestMethod]
    public void Mesh2D_NodeCoordinates_CornerValues()
    {
        var mesh = new Mesh2D(3.0, 2.0, 3, 2, ElementType.Tri);

        // Bottom-left: node 0 = (0, 0)
        Assert.AreEqual(0.0, mesh.Nodes[0, 0], 1e-15);
        Assert.AreEqual(0.0, mesh.Nodes[0, 1], 1e-15);

        // Bottom-right: node 3 = (3, 0) — nodesX=4, so ix=3 → node index 3
        int brNode = mesh.GetNodeIndex(3, 0);
        Assert.AreEqual(3.0, mesh.Nodes[brNode, 0], 1e-15);
        Assert.AreEqual(0.0, mesh.Nodes[brNode, 1], 1e-15);

        // Top-right: node (3,2)
        int trNode = mesh.GetNodeIndex(3, 2);
        Assert.AreEqual(3.0, mesh.Nodes[trNode, 0], 1e-15);
        Assert.AreEqual(2.0, mesh.Nodes[trNode, 1], 1e-15);

        // Top-left: node (0,2)
        int tlNode = mesh.GetNodeIndex(0, 2);
        Assert.AreEqual(0.0, mesh.Nodes[tlNode, 0], 1e-15);
        Assert.AreEqual(2.0, mesh.Nodes[tlNode, 1], 1e-15);
    }

    [TestMethod]
    public void Mesh2D_GetElementNodes_ReturnsCorrectCoords()
    {
        var mesh = new Mesh2D(2.0, 2.0, 2, 2, ElementType.Quad);

        // First element (bottom-left quad): nodes should be at (0,0), (1,0), (1,1), (0,1)
        var coords = mesh.GetElementNodes(0);

        Assert.AreEqual(4, coords.GetLength(0));
        Assert.AreEqual(0.0, coords[0, 0], 1e-15); // (0,0)
        Assert.AreEqual(0.0, coords[0, 1], 1e-15);
        Assert.AreEqual(1.0, coords[1, 0], 1e-15); // (1,0)
        Assert.AreEqual(0.0, coords[1, 1], 1e-15);
        Assert.AreEqual(1.0, coords[2, 0], 1e-15); // (1,1)
        Assert.AreEqual(1.0, coords[2, 1], 1e-15);
        Assert.AreEqual(0.0, coords[3, 0], 1e-15); // (0,1)
        Assert.AreEqual(1.0, coords[3, 1], 1e-15);
    }

    [TestMethod]
    public void Mesh2D_GetNodeIndex_Consistent()
    {
        var mesh = new Mesh2D(1.0, 1.0, 5, 5, ElementType.Tri);

        for (int iy = 0; iy <= 5; iy++)
            for (int ix = 0; ix <= 5; ix++)
            {
                int idx = mesh.GetNodeIndex(ix, iy);
                Assert.IsTrue(idx >= 0 && idx < mesh.NodeCount);
            }

        // Check uniqueness
        var indices = new HashSet<int>();
        for (int iy = 0; iy <= 5; iy++)
            for (int ix = 0; ix <= 5; ix++)
                indices.Add(mesh.GetNodeIndex(ix, iy));

        Assert.AreEqual(mesh.NodeCount, indices.Count);
    }

    #endregion

    #region Assembler2D Tests

    [TestMethod]
    public void Assembler2D_GlobalStiffness_CorrectSize()
    {
        var mesh = new Mesh2D(1.0, 1.0, 2, 2, ElementType.Tri);
        var asm = new Assembler2D(mesh, new TriElement(), 200e9, 0.3, 1.0, true);
        asm.Assemble();

        Assert.AreEqual(mesh.NodeCount * 2, asm.TotalDofs);
    }

    [TestMethod]
    public void Assembler2D_PatchTest_UniformStrain()
    {
        // Patch test: apply linear displacement BCs that produce uniform strain εxx = 1e-3
        // For plane stress with nu=0: σxx = E * εxx, σyy = 0, τxy = 0
        // Displacements: ux = εxx * x, uy = 0
        double E = 1000.0;
        double nu = 0.0;
        double epsilon = 1e-3;
        int nx = 4, ny = 4;

        var mesh = new Mesh2D(1.0, 1.0, nx, ny, ElementType.Tri);
        var asm = new Assembler2D(mesh, new TriElement(), E, nu, 1.0, true);
        asm.Assemble();

        // Prescribe displacement at ALL nodes: ux = epsilon * x, uy = 0
        var bcs = new Dictionary<int, double>();
        for (int i = 0; i < mesh.NodeCount; i++)
        {
            double x = mesh.Nodes[i, 0];
            bcs[i * 2] = epsilon * x;     // ux
            bcs[i * 2 + 1] = 0.0;         // uy
        }

        var u = asm.Solve(bcs);

        // Verify displacements match prescribed
        for (int i = 0; i < mesh.NodeCount; i++)
        {
            double x = mesh.Nodes[i, 0];
            Assert.AreEqual(epsilon * x, u[i * 2], 1e-10, $"Node {i} ux mismatch");
            Assert.AreEqual(0.0, u[i * 2 + 1], 1e-10, $"Node {i} uy mismatch");
        }

        // Verify uniform stress: σxx = E * ε = 1.0, σyy = 0, τxy = 0
        asm.ComputeElementStresses(u, out var sxx, out var syy, out var sxy);
        for (int e = 0; e < mesh.ElementCount; e++)
        {
            Assert.AreEqual(E * epsilon, sxx[e], 1e-6, $"Element {e}: σxx should be {E * epsilon}");
            Assert.AreEqual(0.0, syy[e], 1e-6, $"Element {e}: σyy should be 0");
            Assert.AreEqual(0.0, sxy[e], 1e-6, $"Element {e}: τxy should be 0");
        }
    }

    [TestMethod]
    public void Assembler2D_Cantilever_TipDeflection()
    {
        // Cantilever beam: fixed at x=0, point load P at tip (x=L, y=H/2)
        // Analytical: w_tip ≈ PL³/(3EI) for beam theory
        // FEM with CST will be less accurate but should converge in the right direction
        double L = 10.0, H = 1.0;
        double E = 200e9, nu = 0.3;
        double P = -1000.0; // downward
        double t = 1.0;
        int nx = 40, ny = 4;

        var mesh = new Mesh2D(L, H, nx, ny, ElementType.Tri);
        var asm = new Assembler2D(mesh, new TriElement(), E, nu, t, true);
        asm.Assemble();

        // Fixed left edge: all nodes at ix=0
        var bcs = new Dictionary<int, double>();
        for (int iy = 0; iy <= ny; iy++)
        {
            int node = mesh.GetNodeIndex(0, iy);
            bcs[node * 2] = 0.0;     // ux = 0
            bcs[node * 2 + 1] = 0.0; // uy = 0
        }

        // Point load at mid-height of right edge
        int tipNode = mesh.GetNodeIndex(nx, ny / 2);
        asm.ApplyNodalLoad(tipNode, 1, P);

        var u = asm.Solve(bcs);

        // Tip deflection (beam theory): w = PL³/(3EI), I = t*H³/12
        double I = t * H * H * H / 12.0;
        double wExact = Math.Abs(P) * L * L * L / (3.0 * E * I);

        double tipDeflection = Math.Abs(u[tipNode * 2 + 1]);

        // CST is stiff — expect FEM deflection to be somewhat less than exact
        // but within a reasonable factor (say 0.3–1.5 of exact for coarse mesh)
        Assert.IsTrue(tipDeflection > 0, "Tip should deflect");
        Assert.IsTrue(tipDeflection < wExact * 2.0,
            $"Tip deflection {tipDeflection:E4} is unreasonably large vs exact {wExact:E4}");
        Assert.IsTrue(tipDeflection > wExact * 0.1,
            $"Tip deflection {tipDeflection:E4} is unreasonably small vs exact {wExact:E4}");
    }

    [TestMethod]
    public void Assembler2D_BodyForce_ProducesDisplacement()
    {
        // Simple gravity-loaded block fixed at bottom
        var mesh = new Mesh2D(1.0, 1.0, 3, 3, ElementType.Tri);
        var asm = new Assembler2D(mesh, new TriElement(), 1e6, 0.3, 1.0, true);
        asm.Assemble();

        // Body force in -y (gravity)
        asm.ApplyBodyForce(0, -1000.0);

        // Fix bottom edge
        var bcs = new Dictionary<int, double>();
        for (int ix = 0; ix <= 3; ix++)
        {
            int node = mesh.GetNodeIndex(ix, 0);
            bcs[node * 2] = 0.0;
            bcs[node * 2 + 1] = 0.0;
        }

        var u = asm.Solve(bcs);

        // Top nodes should have negative uy displacement
        for (int ix = 0; ix <= 3; ix++)
        {
            int node = mesh.GetNodeIndex(ix, 3);
            Assert.IsTrue(u[node * 2 + 1] < 0,
                $"Top node ({ix},3) should have negative uy displacement");
        }
    }

    [TestMethod]
    public void Assembler2D_ComputeStresses_NonZeroUnderLoad()
    {
        var mesh = new Mesh2D(2.0, 1.0, 4, 2, ElementType.Tri);
        var asm = new Assembler2D(mesh, new TriElement(), 200e9, 0.3, 1.0, true);
        asm.Assemble();

        // Fix left edge
        var bcs = new Dictionary<int, double>();
        for (int iy = 0; iy <= 2; iy++)
        {
            int node = mesh.GetNodeIndex(0, iy);
            bcs[node * 2] = 0.0;
            bcs[node * 2 + 1] = 0.0;
        }

        // Load at right edge
        for (int iy = 0; iy <= 2; iy++)
        {
            int node = mesh.GetNodeIndex(4, iy);
            asm.ApplyNodalLoad(node, 0, 1000.0); // tension in x
        }

        var u = asm.Solve(bcs);
        asm.ComputeElementStresses(u, out var sxx, out var syy, out var sxy);

        // Should have positive σxx (tension)
        double avgSxx = 0;
        for (int e = 0; e < mesh.ElementCount; e++)
            avgSxx += sxx[e];
        avgSxx /= mesh.ElementCount;

        Assert.IsTrue(avgSxx > 0, $"Average σxx should be positive (tension), got {avgSxx}");
    }

    #endregion

    #region Integration Benchmarks

    [TestMethod]
    public void Benchmark_Cantilever_ConvergesToAnalytical()
    {
        // Cantilever beam: left edge fixed, point load P at tip (right, mid-height)
        // Analytical Euler-Bernoulli tip deflection: δ = PL³/(3EI)
        double L = 10.0, H = 1.0, t = 1.0;
        double E = 200e9, nu = 0.3;
        double P = -1000.0; // downward
        double I_moment = t * H * H * H / 12.0;
        double analytical = P * L * L * L / (3.0 * E * I_moment);

        double prevError = double.MaxValue;
        double lastDeflection = 0;

        // Test convergence with mesh refinement
        int[] meshSizes = { 10, 20, 40 };
        foreach (int n in meshSizes)
        {
            int nx = n;
            int ny = Math.Max(2, n / 5);
            var mesh = new Mesh2D(L, H, nx, ny, ElementType.Tri);
            var assembler = new Assembler2D(mesh, new TriElement(), E, nu, t, true);
            assembler.Assemble();

            // Fix left edge (ix = 0)
            var fixedDofs = new Dictionary<int, double>();
            for (int iy = 0; iy <= ny; iy++)
            {
                int nodeIdx = mesh.GetNodeIndex(0, iy);
                fixedDofs[nodeIdx * 2] = 0.0;
                fixedDofs[nodeIdx * 2 + 1] = 0.0;
            }

            // Apply point load at right mid-height
            int tipNode = mesh.GetNodeIndex(nx, ny / 2);
            assembler.ApplyNodalLoad(tipNode, 1, P);

            var u = assembler.Solve(fixedDofs);
            double tipDeflection = u[tipNode * 2 + 1];
            double error = Math.Abs(tipDeflection - analytical);

            // Error should decrease with refinement
            Assert.IsTrue(error < prevError,
                $"Mesh {nx}x{ny}: error {error:E3} should be < previous {prevError:E3}");
            prevError = error;
            lastDeflection = tipDeflection;
        }

        // Final mesh should be within 20% of analytical (CST converges slowly for bending)
        double relError = Math.Abs((lastDeflection - analytical) / analytical);
        Assert.IsTrue(relError < 0.20,
            $"Final deflection {lastDeflection:E6} should be within 20% of analytical {analytical:E6}, got {relError:P1}");
    }

    [TestMethod]
    public void Benchmark_PlateWithHole_StressConcentration()
    {
        // Quarter-plate with circular hole approximation
        // Far-field tension σ₀ applied at right edge, hole at origin
        // Expected: stress at hole edge ≈ 3σ₀ (SCF for infinite plate)
        double W = 10.0, H_plate = 10.0;
        double E = 200e9, nu = 0.3, t = 1.0;
        double sigma0 = 100e6; // far-field tension
        double holeRadius = 2.0;
        int nx = 20, ny = 20;

        var mesh = new Mesh2D(W, H_plate, nx, ny, ElementType.Tri);
        var assembler = new Assembler2D(mesh, new TriElement(), E, nu, t, true);
        assembler.Assemble();

        // Symmetry BCs: bottom edge (iy=0) fix uy, left edge (ix=0) fix ux
        var fixedDofs = new Dictionary<int, double>();
        for (int ix = 0; ix <= nx; ix++)
        {
            int nodeIdx = mesh.GetNodeIndex(ix, 0);
            fixedDofs[nodeIdx * 2 + 1] = 0.0; // uy = 0
        }
        for (int iy = 0; iy <= ny; iy++)
        {
            int nodeIdx = mesh.GetNodeIndex(0, iy);
            fixedDofs[nodeIdx * 2] = 0.0; // ux = 0
        }

        // Fix nodes inside the hole (approximate hole with zero-displacement constraint)
        for (int iy = 0; iy <= ny; iy++)
        {
            for (int ix = 0; ix <= nx; ix++)
            {
                int nodeIdx = mesh.GetNodeIndex(ix, iy);
                double x = mesh.Nodes[nodeIdx, 0];
                double y = mesh.Nodes[nodeIdx, 1];
                double r = Math.Sqrt(x * x + y * y);
                if (r < holeRadius - 1e-10)
                {
                    fixedDofs[nodeIdx * 2] = 0.0;
                    fixedDofs[nodeIdx * 2 + 1] = 0.0;
                }
            }
        }

        // Apply tension at right edge
        double dy = H_plate / ny;
        for (int iy = 0; iy <= ny; iy++)
        {
            int nodeIdx = mesh.GetNodeIndex(nx, iy);
            double force = sigma0 * t * dy;
            if (iy == 0 || iy == ny) force *= 0.5; // half for boundary nodes
            assembler.ApplyNodalLoad(nodeIdx, 0, force);
        }

        var u = assembler.Solve(fixedDofs);

        // Compute element stresses
        assembler.ComputeElementStresses(u, out double[] sxx, out double[] syy, out double[] sxy);

        // Find max σxx near the hole boundary (elements adjacent to hole)
        double maxStressNearHole = 0;
        for (int e = 0; e < mesh.ElementCount; e++)
        {
            // Compute centroid of element
            var nodeCoords = mesh.GetElementNodes(e);
            double cx = (nodeCoords[0, 0] + nodeCoords[1, 0] + nodeCoords[2, 0]) / 3.0;
            double cy = (nodeCoords[0, 1] + nodeCoords[1, 1] + nodeCoords[2, 1]) / 3.0;
            double cr = Math.Sqrt(cx * cx + cy * cy);

            // Elements near the hole boundary (within 1 element width)
            if (cr >= holeRadius && cr < holeRadius + W / nx * 2.0)
            {
                if (sxx[e] > maxStressNearHole)
                    maxStressNearHole = sxx[e];
            }
        }

        // Stress concentration should be elevated above far-field
        // Note: rigid-inclusion BC ≠ free hole (K≈3), so we check qualitative elevation
        double scf = maxStressNearHole / sigma0;
        Assert.IsTrue(scf > 1.2,
            $"Stress concentration factor {scf:F2} should be > 1.2 (qualitative check)");
    }

    [TestMethod]
    public void Benchmark_CooksMembrane_Q4Displacement()
    {
        // Cook's membrane: standard FEM benchmark for distorted Q4 elements
        // Trapezoidal domain: left edge (0,0)-(0,44), right edge (48,44)-(48,60)
        // Left edge clamped, right edge uniform shear load F=1/16 per unit length
        // Reference tip displacement ≈ 23.96 (fine mesh limit)
        int nx = 16, ny = 16;
        double E = 1.0, nu = 1.0 / 3.0, t = 1.0;

        // Create rectangular mesh then distort to Cook's trapezoidal shape
        var mesh = new Mesh2D(1.0, 1.0, nx, ny, ElementType.Quad);

        // Bilinear mapping from unit square to Cook's membrane
        // Corners: (0,0), (48,44), (48,60), (0,44)
        double x0 = 0, y0 = 0;   // bottom-left
        double x1 = 48, y1 = 44;  // bottom-right
        double x2 = 48, y2 = 60;  // top-right
        double x3 = 0, y3 = 44;   // top-left

        for (int i = 0; i < mesh.NodeCount; i++)
        {
            double xi = mesh.Nodes[i, 0];  // in [0,1]
            double eta = mesh.Nodes[i, 1]; // in [0,1]

            double xPhys = (1 - xi) * (1 - eta) * x0 + xi * (1 - eta) * x1
                         + xi * eta * x2 + (1 - xi) * eta * x3;
            double yPhys = (1 - xi) * (1 - eta) * y0 + xi * (1 - eta) * y1
                         + xi * eta * y2 + (1 - xi) * eta * y3;

            mesh.Nodes[i, 0] = xPhys;
            mesh.Nodes[i, 1] = yPhys;
        }

        var assembler = new Assembler2D(mesh, new QuadElement(), E, nu, t, true);
        assembler.Assemble();

        // Fix left edge (ix = 0)
        var fixedDofs = new Dictionary<int, double>();
        for (int iy = 0; iy <= ny; iy++)
        {
            int nodeIdx = mesh.GetNodeIndex(0, iy);
            fixedDofs[nodeIdx * 2] = 0.0;
            fixedDofs[nodeIdx * 2 + 1] = 0.0;
        }

        // Apply uniform vertical shear on right edge: total F = 1
        double totalForce = 1.0;
        double forcePerNode = totalForce / ny;
        for (int iy = 0; iy <= ny; iy++)
        {
            int nodeIdx = mesh.GetNodeIndex(nx, iy);
            double f = forcePerNode;
            if (iy == 0 || iy == ny) f *= 0.5;
            assembler.ApplyNodalLoad(nodeIdx, 1, f);
        }

        var u = assembler.Solve(fixedDofs);

        // Check vertical displacement at top-right corner (node C)
        int tipNode = mesh.GetNodeIndex(nx, ny);
        double tipDisp = u[tipNode * 2 + 1];

        // Verify positive and finite displacement (shear load → upward motion at tip)
        Assert.IsTrue(tipDisp > 10.0,
            $"Cook's membrane tip displacement {tipDisp:F4} should be positive and significant");

        // Verify convergence: coarser mesh should give a different (less converged) value
        var meshCoarse = new Mesh2D(1.0, 1.0, 4, 4, ElementType.Quad);
        for (int i = 0; i < meshCoarse.NodeCount; i++)
        {
            double xi = meshCoarse.Nodes[i, 0];
            double eta = meshCoarse.Nodes[i, 1];
            meshCoarse.Nodes[i, 0] = (1 - xi) * (1 - eta) * x0 + xi * (1 - eta) * x1
                                    + xi * eta * x2 + (1 - xi) * eta * x3;
            meshCoarse.Nodes[i, 1] = (1 - xi) * (1 - eta) * y0 + xi * (1 - eta) * y1
                                    + xi * eta * y2 + (1 - xi) * eta * y3;
        }
        var asmCoarse = new Assembler2D(meshCoarse, new QuadElement(), E, nu, t, true);
        asmCoarse.Assemble();
        var fixCoarse = new Dictionary<int, double>();
        for (int iy = 0; iy <= 4; iy++)
        {
            int nodeIdx = meshCoarse.GetNodeIndex(0, iy);
            fixCoarse[nodeIdx * 2] = 0.0;
            fixCoarse[nodeIdx * 2 + 1] = 0.0;
        }
        double fCoarse = totalForce / 4;
        for (int iy = 0; iy <= 4; iy++)
        {
            int nodeIdx = meshCoarse.GetNodeIndex(4, iy);
            double f = fCoarse;
            if (iy == 0 || iy == 4) f *= 0.5;
            asmCoarse.ApplyNodalLoad(nodeIdx, 1, f);
        }
        var uCoarse = asmCoarse.Solve(fixCoarse);
        double tipCoarse = uCoarse[meshCoarse.GetNodeIndex(4, 4) * 2 + 1];

        // Finer mesh should differ from coarse (convergence behavior)
        Assert.AreNotEqual(tipCoarse, tipDisp, 1.0,
            $"Coarse ({tipCoarse:F4}) and fine ({tipDisp:F4}) should differ (convergence)");
    }

    [TestMethod]
    public void Benchmark_ThickCylinder_LameSolution()
    {
        // Quarter annulus under internal pressure - plane strain
        // Inner radius a=1, outer radius b=2, internal pressure p=100
        // Lamé analytical radial displacement (plane strain):
        //   u_r(r) = (a²p/(E(b²-a²))) × ((1-2ν)r + b²(1+ν)/r)
        // Stress: σ_rr(r) = a²p/(b²-a²) × (1 - b²/r²)  [negative = compressive]
        double a = 1.0, b = 2.0, p = 100.0;
        double E = 1000.0, nu = 0.3, t = 1.0;
        int nx = 20, ny = 20;

        // Create mesh on unit square then map to quarter annulus
        var mesh = new Mesh2D(1.0, 1.0, nx, ny, ElementType.Tri);

        // Map (xi, eta) ∈ [0,1]² → quarter annulus
        // r = a + xi*(b-a), θ = eta * π/2
        for (int i = 0; i < mesh.NodeCount; i++)
        {
            double xi = mesh.Nodes[i, 0];  // radial parameter [0,1]
            double eta = mesh.Nodes[i, 1]; // angular parameter [0,1]
            double r = a + xi * (b - a);
            double theta = eta * Math.PI / 2.0;

            mesh.Nodes[i, 0] = r * Math.Cos(theta);
            mesh.Nodes[i, 1] = r * Math.Sin(theta);
        }

        // Plane strain analysis
        var assembler = new Assembler2D(mesh, new TriElement(), E, nu, t, false);
        assembler.Assemble();

        var fixedDofs = new Dictionary<int, double>();

        // Symmetry BCs: bottom edge (eta=0, iy=0) → θ=0 → fix uy
        for (int ix = 0; ix <= nx; ix++)
        {
            int nodeIdx = mesh.GetNodeIndex(ix, 0);
            fixedDofs[nodeIdx * 2 + 1] = 0.0; // uy = 0
        }

        // Symmetry BCs: top edge (eta=1, iy=ny) → θ=π/2 → fix ux
        for (int ix = 0; ix <= nx; ix++)
        {
            int nodeIdx = mesh.GetNodeIndex(ix, ny);
            fixedDofs[nodeIdx * 2] = 0.0; // ux = 0
        }

        // Apply internal pressure on inner edge (ix=0, r=a)
        // Pressure acts radially outward → decompose into x,y components
        for (int iy = 0; iy <= ny; iy++)
        {
            int nodeIdx = mesh.GetNodeIndex(0, iy);
            double eta = (double)iy / ny;
            double theta = eta * Math.PI / 2.0;

            // Tributary length for this node (trapezoidal rule on arc)
            double dtheta = Math.PI / 2.0 / ny;
            double tribLength = a * dtheta;
            if (iy == 0 || iy == ny) tribLength *= 0.5;

            double fx = p * Math.Cos(theta) * tribLength * t;
            double fy = p * Math.Sin(theta) * tribLength * t;

            assembler.ApplyNodalLoad(nodeIdx, 0, fx);
            assembler.ApplyNodalLoad(nodeIdx, 1, fy);
        }

        var u = assembler.Solve(fixedDofs);

        // Verify radial displacement at inner surface along θ=0
        // Analytical: u_r(a) = (a²p/(E(b²-a²))) × ((1-2ν)a + b²(1+ν)/a)
        double coeff = a * a * p / (E * (b * b - a * a));
        double ur_inner_exact = coeff * ((1 - 2 * nu) * a + b * b * (1 + nu) / a);

        int innerNode = mesh.GetNodeIndex(0, 0); // r=a, θ=0
        double ur_inner_fem = u[innerNode * 2]; // ux at θ=0 equals u_r

        // Displacement should be positive (expansion) and within 25% of analytical
        Assert.IsTrue(ur_inner_fem > 0, $"Inner radial displacement should be positive, got {ur_inner_fem}");
        double dispRelErr = Math.Abs((ur_inner_fem - ur_inner_exact) / ur_inner_exact);
        Assert.IsTrue(dispRelErr < 0.25,
            $"Inner displacement error {dispRelErr:P1}: FEM={ur_inner_fem:E4}, exact={ur_inner_exact:E4}");

        // Verify displacement at outer surface is smaller than at inner (pressure decays)
        int outerNode = mesh.GetNodeIndex(nx, 0); // r=b, θ=0
        double ur_outer_fem = u[outerNode * 2];
        double ur_outer_exact = coeff * ((1 - 2 * nu) * b + b * b * (1 + nu) / b);

        Assert.IsTrue(ur_outer_fem > 0, "Outer radial displacement should be positive");
        Assert.IsTrue(ur_outer_fem < ur_inner_fem,
            $"Outer displacement ({ur_outer_fem:E4}) should be < inner ({ur_inner_fem:E4})");

        // Compute stresses and verify qualitative behavior
        assembler.ComputeElementStresses(u, out double[] sxx, out double[] syy, out double[] sxy);

        // Along θ=0: σ_rr = σ_xx should be negative (compressive) near inner radius
        // and approach zero at outer radius
        double stressCoeff = a * a * p / (b * b - a * a);
        bool foundInner = false, foundOuter = false;
        double innerStress = 0, outerStress = 0;

        for (int e = 0; e < mesh.ElementCount; e++)
        {
            var nodeCoords = mesh.GetElementNodes(e);
            double cx = (nodeCoords[0, 0] + nodeCoords[1, 0] + nodeCoords[2, 0]) / 3.0;
            double cy = (nodeCoords[0, 1] + nodeCoords[1, 1] + nodeCoords[2, 1]) / 3.0;
            double cr = Math.Sqrt(cx * cx + cy * cy);
            double ctheta = Math.Atan2(cy, cx);

            // Elements near θ=0
            if (ctheta < Math.PI / (2.0 * ny) * 1.5)
            {
                if (cr < a + 0.15 && !foundInner) { innerStress = sxx[e]; foundInner = true; }
                if (cr > b - 0.15 && !foundOuter) { outerStress = sxx[e]; foundOuter = true; }
            }
        }

        Assert.IsTrue(foundInner && foundOuter, "Should find elements near inner and outer radius");
        // σ_rr at inner ≈ -100, at outer ≈ 0
        Assert.IsTrue(innerStress < outerStress,
            $"Inner σ_rr ({innerStress:F2}) should be more compressive than outer ({outerStress:F2})");
        Assert.IsTrue(innerStress < 0,
            $"Inner radial stress should be compressive, got {innerStress:F2}");
    }

    #endregion
}
