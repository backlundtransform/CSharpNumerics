using CSharpNumerics.Numerics;
using CSharpNumerics.Numerics.FiniteDifference;
using CSharpNumerics.Numerics.Objects;

namespace NumericsTests
{
    [TestClass]
    public class FiniteDifferenceTests
    {
        // ════════════════════════════════════════════
        //  Grid2D
        // ════════════════════════════════════════════

        [TestMethod]
        public void Grid2D_Index_RoundTrips()
        {
            var grid = new Grid2D(5, 4, 1.0);
            for (int iy = 0; iy < 4; iy++)
                for (int ix = 0; ix < 5; ix++)
                {
                    int flat = grid.Index(ix, iy);
                    var (rx, ry) = grid.Index2D(flat);
                    Assert.AreEqual(ix, rx);
                    Assert.AreEqual(iy, ry);
                }
        }

        [TestMethod]
        public void Grid2D_ToVector_ToArray_RoundTrips()
        {
            var grid = new Grid2D(3, 3, 0.5);
            var field = new double[3, 3];
            field[1, 1] = 42.0;
            field[0, 2] = -7.0;

            var v = grid.ToVector(field);
            var back = grid.ToArray(v);

            Assert.AreEqual(42.0, back[1, 1]);
            Assert.AreEqual(-7.0, back[0, 2]);
            Assert.AreEqual(0.0, back[0, 0]);
        }

        [TestMethod]
        public void Grid2D_Initialize_SetsValues()
        {
            var grid = new Grid2D(4, 4, 0.25);
            var v = grid.Initialize((x, y) => x + y);

            // Cell (2,3) → x=0.5, y=0.75 → 1.25
            Assert.AreEqual(1.25, v[grid.Index(2, 3)], 1e-12);
        }

        // ════════════════════════════════════════════
        //  1D Laplacian
        // ════════════════════════════════════════════

        [TestMethod]
        public void Laplacian1D_Quadratic_IsConstant()
        {
            // u = x², d²u/dx² = 2 everywhere
            int n = 50;
            double dx = 0.1;
            var values = new double[n];
            for (int i = 0; i < n; i++)
                values[i] = (i * dx) * (i * dx);

            var u = new VectorN(values);
            var lap = GridOperators.Laplacian1D(u, dx, BoundaryCondition.Dirichlet);

            // Check interior points (skip boundaries which are affected by BC)
            for (int i = 2; i < n - 2; i++)
                Assert.AreEqual(2.0, lap[i], 0.01, $"Laplacian at i={i}");
        }

        [TestMethod]
        public void Laplacian1D_Periodic_Wraps()
        {
            // Sine wave: u = sin(2πx/L), d²u/dx² = -(2π/L)² sin(2πx/L)
            int n = 100;
            double L = n * 0.1;
            double dx = 0.1;
            double k = 2 * Math.PI / L;
            var values = new double[n];
            for (int i = 0; i < n; i++)
                values[i] = Math.Sin(k * i * dx);

            var u = new VectorN(values);
            var lap = GridOperators.Laplacian1D(u, dx, BoundaryCondition.Periodic);

            // Should be ≈ -k² sin(kx)
            for (int i = 0; i < n; i++)
            {
                double expected = -k * k * Math.Sin(k * i * dx);
                Assert.AreEqual(expected, lap[i], 0.05, $"Periodic Laplacian at i={i}");
            }
        }

        // ════════════════════════════════════════════
        //  1D Gradient
        // ════════════════════════════════════════════

        [TestMethod]
        public void Gradient1D_Linear_IsConstant()
        {
            // u = 3x, du/dx = 3
            int n = 20;
            double dx = 0.5;
            var values = new double[n];
            for (int i = 0; i < n; i++)
                values[i] = 3.0 * i * dx;

            var u = new VectorN(values);
            var grad = GridOperators.Gradient1D(u, dx, BoundaryCondition.Neumann);

            for (int i = 1; i < n - 1; i++)
                Assert.AreEqual(3.0, grad[i], 1e-10, $"Gradient at i={i}");
        }

        // ════════════════════════════════════════════
        //  2D Laplacian
        // ════════════════════════════════════════════

        [TestMethod]
        public void Laplacian2D_Quadratic_IsConstant()
        {
            // u = x² + y², ∇²u = 4
            int nx = 20, ny = 20;
            double d = 0.1;
            var grid = new Grid2D(nx, ny, d);
            var u = grid.Initialize((x, y) => x * x + y * y);

            var lap = GridOperators.Laplacian2D(u, grid, BoundaryCondition.Dirichlet);

            // Check interior (skip 2-cell border)
            for (int iy = 2; iy < ny - 2; iy++)
                for (int ix = 2; ix < nx - 2; ix++)
                    Assert.AreEqual(4.0, lap[grid.Index(ix, iy)], 0.01,
                        $"Laplacian at ({ix},{iy})");
        }

        [TestMethod]
        public void Laplacian2D_Neumann_NoFlux()
        {
            // Constant field: ∇²u = 0
            int nx = 10, ny = 10;
            var grid = new Grid2D(nx, ny, 1.0);
            var u = grid.Initialize((x, y) => 5.0);

            var lap = GridOperators.Laplacian2D(u, grid, BoundaryCondition.Neumann);

            for (int i = 0; i < grid.Length; i++)
                Assert.AreEqual(0.0, lap[i], 1e-12, $"Laplacian at flat index {i}");
        }

        // ════════════════════════════════════════════
        //  2D Gradient
        // ════════════════════════════════════════════

        [TestMethod]
        public void Gradient2D_Linear_ReturnsConstants()
        {
            // u = 2x + 3y → ∂u/∂x = 2, ∂u/∂y = 3
            int nx = 20, ny = 20;
            double d = 0.1;
            var grid = new Grid2D(nx, ny, d);
            var u = grid.Initialize((x, y) => 2.0 * x + 3.0 * y);

            var (dux, duy) = GridOperators.Gradient2D(u, grid, BoundaryCondition.Neumann);

            for (int iy = 1; iy < ny - 1; iy++)
                for (int ix = 1; ix < nx - 1; ix++)
                {
                    int idx = grid.Index(ix, iy);
                    Assert.AreEqual(2.0, dux[idx], 0.01, $"du/dx at ({ix},{iy})");
                    Assert.AreEqual(3.0, duy[idx], 0.01, $"du/dy at ({ix},{iy})");
                }
        }

        // ════════════════════════════════════════════
        //  2D Divergence
        // ════════════════════════════════════════════

        [TestMethod]
        public void Divergence2D_LinearField_IsConstant()
        {
            // F = (2x, 3y) → ∇·F = 2 + 3 = 5
            int nx = 20, ny = 20;
            double d = 0.1;
            var grid = new Grid2D(nx, ny, d);
            var fx = grid.Initialize((x, y) => 2.0 * x);
            var fy = grid.Initialize((x, y) => 3.0 * y);

            var div = GridOperators.Divergence2D(fx, fy, grid, BoundaryCondition.Neumann);

            for (int iy = 1; iy < ny - 1; iy++)
                for (int ix = 1; ix < nx - 1; ix++)
                    Assert.AreEqual(5.0, div[grid.Index(ix, iy)], 0.01,
                        $"Divergence at ({ix},{iy})");
        }

        // ════════════════════════════════════════════
        //  Advection
        // ════════════════════════════════════════════

        [TestMethod]
        public void Advection2D_UniformVelocity_ShiftsField()
        {
            // Constant velocity (1,0) advecting u = x → v·∇u = 1
            int nx = 20, ny = 10;
            double d = 0.1;
            var grid = new Grid2D(nx, ny, d);
            var u = grid.Initialize((x, y) => x);
            var vx = grid.Initialize((x, y) => 1.0);
            var vy = grid.Zeros();

            var adv = GridOperators.Advection2D(u, vx, vy, grid, BoundaryCondition.Neumann);

            // Interior: upwind gives (u[i] - u[i-1])/dx = 1
            for (int iy = 1; iy < ny - 1; iy++)
                for (int ix = 2; ix < nx - 1; ix++)
                    Assert.AreEqual(1.0, adv[grid.Index(ix, iy)], 0.01,
                        $"Advection at ({ix},{iy})");
        }

        // ════════════════════════════════════════════
        //  Integration test: Heat equation with RK4
        // ════════════════════════════════════════════

        [TestMethod]
        public void HeatEquation_1D_Dirichlet_Decays()
        {
            // ∂u/∂t = α·∂²u/∂x², u(0)=u(L)=0
            // Initial hot spot in the middle → should decay toward zero
            double alpha = 0.01;
            int n = 51; // odd so centre is exactly at n/2
            double dx = 0.02;
            var values = new double[n];
            values[n / 2] = 100.0; // hot spot

            var u0 = new VectorN(values);

            Func<(double t, VectorN y), VectorN> rhs = args =>
                alpha * GridOperators.Laplacian1D(args.y, dx, BoundaryCondition.Dirichlet);

            // CFL: dt < dx²/(2·α) = 0.02
            var result = rhs.RungeKutta(0, 1.0, 0.005, u0);

            // Energy should decrease (peak < initial)
            double maxVal = 0;
            for (int i = 0; i < n; i++)
                maxVal = Math.Max(maxVal, result[i]);

            Assert.IsTrue(maxVal < 100.0, "Peak should have decreased from diffusion");
            Assert.IsTrue(maxVal > 0.0, "Solution should still be positive");

            // u should be symmetric around centre
            for (int i = 0; i < n / 2; i++)
                Assert.AreEqual(result[i], result[n - 1 - i], 0.05,
                    $"Symmetry at i={i}");
        }

        [TestMethod]
        public void HeatEquation_2D_Dirichlet_Decays()
        {
            // 2D heat equation: ∂u/∂t = α·∇²u
            double alpha = 0.01;
            int nx = 20, ny = 20;
            double d = 0.1;
            var grid = new Grid2D(nx, ny, d);

            // Hot spot in centre
            var field = new double[nx, ny];
            field[nx / 2, ny / 2] = 100.0;
            var u0 = grid.ToVector(field);

            Func<(double t, VectorN y), VectorN> rhs = args =>
                alpha * GridOperators.Laplacian2D(args.y, grid, BoundaryCondition.Dirichlet);

            // dt < dx²/(4·α) = 0.025 for 2D stability
            var result = rhs.RungeKutta(0, 0.5, 0.005, u0);

            double maxVal = 0;
            for (int i = 0; i < grid.Length; i++)
                maxVal = Math.Max(maxVal, result[i]);

            Assert.IsTrue(maxVal < 100.0, "2D peak should decay");
            Assert.IsTrue(maxVal > 0.0, "Solution should still be positive");
        }

        // ════════════════════════════════════════════
        //  Grid3D
        // ════════════════════════════════════════════

        [TestMethod]
        public void Grid3D_Index_RoundTrips()
        {
            var grid = new Grid3D(5, 4, 3, 1.0);
            for (int iz = 0; iz < 3; iz++)
                for (int iy = 0; iy < 4; iy++)
                    for (int ix = 0; ix < 5; ix++)
                    {
                        int flat = grid.Index(ix, iy, iz);
                        var (rx, ry, rz) = grid.Index3D(flat);
                        Assert.AreEqual(ix, rx, $"ix at flat={flat}");
                        Assert.AreEqual(iy, ry, $"iy at flat={flat}");
                        Assert.AreEqual(iz, rz, $"iz at flat={flat}");
                    }
        }

        [TestMethod]
        public void Grid3D_ToVector_ToArray_RoundTrips()
        {
            var grid = new Grid3D(3, 3, 3, 0.5);
            var field = new double[3, 3, 3];
            field[1, 1, 1] = 42.0;
            field[0, 2, 1] = -7.0;

            var v = grid.ToVector(field);
            var back = grid.ToArray(v);

            Assert.AreEqual(42.0, back[1, 1, 1]);
            Assert.AreEqual(-7.0, back[0, 2, 1]);
            Assert.AreEqual(0.0, back[0, 0, 0]);
        }

        [TestMethod]
        public void Grid3D_Initialize_SetsCorrectPositions()
        {
            var grid = new Grid3D(4, 4, 4, 0.25);
            var v = grid.Initialize((x, y, z) => x + y + z);

            // Cell (2,3,1) → x=0.5, y=0.75, z=0.25 → 1.5
            Assert.AreEqual(1.5, v[grid.Index(2, 3, 1)], 1e-12);
        }

        [TestMethod]
        public void Grid3D_Length_IsNxTimesNyTimesNz()
        {
            var grid = new Grid3D(5, 4, 3, 1.0);
            Assert.AreEqual(60, grid.Length);
        }

        // ════════════════════════════════════════════
        //  3D Laplacian
        // ════════════════════════════════════════════

        [TestMethod]
        public void Laplacian3D_Quadratic_IsConstant()
        {
            // u = x² + y² + z², ∇²u = 2 + 2 + 2 = 6
            int n = 15;
            double d = 0.1;
            var grid = new Grid3D(n, n, n, d);
            var u = grid.Initialize((x, y, z) => x * x + y * y + z * z);

            var lap = GridOperators3D.Laplacian3D(u, grid, BoundaryCondition.Dirichlet);

            // Check interior (skip 2-cell border)
            for (int iz = 2; iz < n - 2; iz++)
                for (int iy = 2; iy < n - 2; iy++)
                    for (int ix = 2; ix < n - 2; ix++)
                        Assert.AreEqual(6.0, lap[grid.Index(ix, iy, iz)], 0.01,
                            $"Laplacian at ({ix},{iy},{iz})");
        }

        [TestMethod]
        public void Laplacian3D_Linear_IsZero()
        {
            // u = 2x + 3y + 5z, ∇²u = 0 at interior cells
            int n = 10;
            double d = 0.1;
            var grid = new Grid3D(n, n, n, d);
            var u = grid.Initialize((x, y, z) => 2 * x + 3 * y + 5 * z);

            var lap = GridOperators3D.Laplacian3D(u, grid, BoundaryCondition.Neumann);

            // Neumann ghost cells produce one-sided differences at boundaries,
            // so only check interior (skip first/last in each dimension).
            for (int iz = 1; iz < n - 1; iz++)
                for (int iy = 1; iy < n - 1; iy++)
                    for (int ix = 1; ix < n - 1; ix++)
                        Assert.AreEqual(0.0, lap[grid.Index(ix, iy, iz)], 1e-10,
                            $"Laplacian at ({ix},{iy},{iz});");
        }

        [TestMethod]
        public void Laplacian3D_Periodic_Wraps()
        {
            // Constant field → ∇²u = 0 with periodic BCs
            int n = 8;
            var grid = new Grid3D(n, n, n, 1.0);
            var u = grid.Initialize((x, y, z) => 5.0);

            var lap = GridOperators3D.Laplacian3D(u, grid, BoundaryCondition.Periodic);

            for (int i = 0; i < grid.Length; i++)
                Assert.AreEqual(0.0, lap[i], 1e-12, $"Periodic Laplacian at flat={i}");
        }

        [TestMethod]
        public void Laplacian3D_Neumann_MirrorsInterior()
        {
            // Constant field → ∇²u = 0 with Neumann
            int n = 8;
            var grid = new Grid3D(n, n, n, 1.0);
            var u = grid.Initialize((x, y, z) => 7.0);

            var lap = GridOperators3D.Laplacian3D(u, grid, BoundaryCondition.Neumann);

            for (int i = 0; i < grid.Length; i++)
                Assert.AreEqual(0.0, lap[i], 1e-12);
        }

        // ════════════════════════════════════════════
        //  3D Gradient
        // ════════════════════════════════════════════

        [TestMethod]
        public void Gradient3D_Linear_ReturnsCoefficients()
        {
            // u = 2x + 3y + 5z → ∇u = (2, 3, 5)
            int n = 15;
            double d = 0.1;
            var grid = new Grid3D(n, n, n, d);
            var u = grid.Initialize((x, y, z) => 2 * x + 3 * y + 5 * z);

            var (dux, duy, duz) = GridOperators3D.Gradient3D(u, grid, BoundaryCondition.Neumann);

            for (int iz = 1; iz < n - 1; iz++)
                for (int iy = 1; iy < n - 1; iy++)
                    for (int ix = 1; ix < n - 1; ix++)
                    {
                        int idx = grid.Index(ix, iy, iz);
                        Assert.AreEqual(2.0, dux[idx], 0.01, $"du/dx at ({ix},{iy},{iz})");
                        Assert.AreEqual(3.0, duy[idx], 0.01, $"du/dy at ({ix},{iy},{iz})");
                        Assert.AreEqual(5.0, duz[idx], 0.01, $"du/dz at ({ix},{iy},{iz})");
                    }
        }

        // ════════════════════════════════════════════
        //  3D Divergence
        // ════════════════════════════════════════════

        [TestMethod]
        public void Divergence3D_IdentityField_IsThree()
        {
            // F = (x, y, z) → ∇·F = 1 + 1 + 1 = 3
            int n = 15;
            double d = 0.1;
            var grid = new Grid3D(n, n, n, d);
            var fx = grid.Initialize((x, y, z) => x);
            var fy = grid.Initialize((x, y, z) => y);
            var fz = grid.Initialize((x, y, z) => z);

            var div = GridOperators3D.Divergence3D(fx, fy, fz, grid, BoundaryCondition.Neumann);

            for (int iz = 1; iz < n - 1; iz++)
                for (int iy = 1; iy < n - 1; iy++)
                    for (int ix = 1; ix < n - 1; ix++)
                        Assert.AreEqual(3.0, div[grid.Index(ix, iy, iz)], 0.01,
                            $"Divergence at ({ix},{iy},{iz})");
        }

        // ════════════════════════════════════════════
        //  3D Advection
        // ════════════════════════════════════════════

        [TestMethod]
        public void Advection3D_UniformFlow_ShiftsProfile()
        {
            // v = (1,0,0), u = x → v·∇u = 1
            int n = 15;
            double d = 0.1;
            var grid = new Grid3D(n, n, n, d);
            var u = grid.Initialize((x, y, z) => x);
            var vx = grid.Initialize((x, y, z) => 1.0);
            var vy = grid.Zeros();
            var vz = grid.Zeros();

            var adv = GridOperators3D.Advection3D(u, vx, vy, vz, grid, BoundaryCondition.Neumann);

            for (int iz = 1; iz < n - 1; iz++)
                for (int iy = 1; iy < n - 1; iy++)
                    for (int ix = 2; ix < n - 1; ix++)
                        Assert.AreEqual(1.0, adv[grid.Index(ix, iy, iz)], 0.01,
                            $"Advection at ({ix},{iy},{iz})");
        }

        // ════════════════════════════════════════════
        //  3D Poisson solver
        // ════════════════════════════════════════════

        [TestMethod]
        public void SolvePoisson3D_KnownSolution_Converges()
        {
            // ∇²u = 0 with u=0 on all faces except top (iz=nz-1)=100
            // → linear gradient in z: u ≈ 100 * iz / (nz-1)
            int n = 10;
            double d = 1.0 / (n - 1);
            var grid = new Grid3D(n, n, n, d);

            var rhs = grid.Zeros();
            var mask = new bool[grid.Length];
            var bcVals = new double[grid.Length];

            // All 6 faces are Dirichlet with values consistent with
            // the expected linear solution u = 100 * z.
            for (int iz = 0; iz < n; iz++)
                for (int iy = 0; iy < n; iy++)
                    for (int ix = 0; ix < n; ix++)
                    {
                        if (ix == 0 || ix == n - 1 || iy == 0 || iy == n - 1 || iz == 0 || iz == n - 1)
                        {
                            int idx = grid.Index(ix, iy, iz);
                            mask[idx] = true;
                            bcVals[idx] = 100.0 * iz * d; // linear in z
                        }
                    }

            var (sol, iters) = GridOperators3D.SolvePoisson3D(rhs, grid, mask, bcVals);

            Assert.IsTrue(iters < 10000, "Should converge.");

            // Check interior: u ≈ 100 * z (linear gradient)
            for (int iz = 1; iz < n - 1; iz++)
            {
                double expected = 100.0 * iz * d;
                // Check centre column
                int idx = grid.Index(n / 2, n / 2, iz);
                Assert.AreEqual(expected, sol[idx], 2.0,
                    $"Poisson solution at iz={iz}");
            }
        }

        // ════════════════════════════════════════════
        //  3D Heat equation integration
        // ════════════════════════════════════════════

        [TestMethod]
        public void HeatEquation_3D_Dirichlet_Decays()
        {
            // ∂u/∂t = α·∇²u, hot spot in centre → should decay
            double alpha = 0.01;
            int n = 10;
            double d = 0.1;
            var grid = new Grid3D(n, n, n, d);

            var field = new double[n, n, n];
            field[n / 2, n / 2, n / 2] = 100.0;
            var u0 = grid.ToVector(field);

            // Forward Euler: dt < dx²/(6α) ≈ 0.0167 for 3D
            double dt = 0.01;
            int steps = 50;
            var u = u0;

            for (int s = 0; s < steps; s++)
            {
                var lap = GridOperators3D.Laplacian3D(u, grid, BoundaryCondition.Dirichlet);
                var next = new double[grid.Length];
                for (int i = 0; i < grid.Length; i++)
                    next[i] = u[i] + dt * alpha * lap[i];
                u = new VectorN(next);
            }

            double maxVal = 0;
            for (int i = 0; i < grid.Length; i++)
                maxVal = Math.Max(maxVal, u[i]);

            Assert.IsTrue(maxVal < 100.0, "3D peak should decay from diffusion");
            Assert.IsTrue(maxVal > 0.0, "Solution should still be positive");
        }
    }
}

