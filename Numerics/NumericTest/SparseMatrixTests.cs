using CSharpNumerics.Numerics.Objects;

namespace NumericTest;

[TestClass]
public class SparseMatrixTests
{
    #region Construction & Properties

    [TestMethod]
    public void FromTriplets_BasicConstruction_CorrectDimensions()
    {
        var triplets = new List<(int, int, double)>
        {
            (0, 0, 2.0), (0, 1, -1.0),
            (1, 0, -1.0), (1, 1, 2.0), (1, 2, -1.0),
            (2, 1, -1.0), (2, 2, 2.0)
        };

        var A = SparseMatrix.FromTriplets(3, 3, triplets);

        Assert.AreEqual(3, A.Rows);
        Assert.AreEqual(3, A.Cols);
        Assert.AreEqual(7, A.NonZeroCount);
    }

    [TestMethod]
    public void FromTriplets_AccumulatesDuplicates()
    {
        var triplets = new List<(int, int, double)>
        {
            (0, 0, 1.0),
            (0, 0, 3.0),
            (1, 1, 5.0)
        };

        var A = SparseMatrix.FromTriplets(2, 2, triplets);

        Assert.AreEqual(4.0, A.Get(0, 0), 1e-15);
        Assert.AreEqual(5.0, A.Get(1, 1), 1e-15);
        Assert.AreEqual(2, A.NonZeroCount);
    }

    [TestMethod]
    public void FromTriplets_EmptyRows_Handled()
    {
        // Row 1 has no entries
        var triplets = new List<(int, int, double)>
        {
            (0, 0, 1.0),
            (2, 2, 3.0)
        };

        var A = SparseMatrix.FromTriplets(3, 3, triplets);

        Assert.AreEqual(1.0, A.Get(0, 0), 1e-15);
        Assert.AreEqual(0.0, A.Get(1, 0), 1e-15);
        Assert.AreEqual(0.0, A.Get(1, 1), 1e-15);
        Assert.AreEqual(3.0, A.Get(2, 2), 1e-15);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentOutOfRangeException))]
    public void FromTriplets_OutOfBounds_Throws()
    {
        var triplets = new List<(int, int, double)> { (3, 0, 1.0) };
        SparseMatrix.FromTriplets(3, 3, triplets);
    }

    #endregion

    #region Multiply

    [TestMethod]
    public void Multiply_MatchesDenseResult()
    {
        // Dense:
        // [2 -1  0]   [1]   [ 1]
        // [-1 2 -1] × [1] = [ 0]
        // [0 -1  2]   [1]   [ 1]
        var triplets = new List<(int, int, double)>
        {
            (0, 0, 2.0), (0, 1, -1.0),
            (1, 0, -1.0), (1, 1, 2.0), (1, 2, -1.0),
            (2, 1, -1.0), (2, 2, 2.0)
        };
        var A = SparseMatrix.FromTriplets(3, 3, triplets);
        var x = new VectorN(new[] { 1.0, 1.0, 1.0 });

        var y = A.Multiply(x);

        Assert.AreEqual(1.0, y[0], 1e-15);
        Assert.AreEqual(0.0, y[1], 1e-15);
        Assert.AreEqual(1.0, y[2], 1e-15);
    }

    [TestMethod]
    public void Multiply_IdentityMatrix_ReturnsInput()
    {
        var triplets = new List<(int, int, double)>
        {
            (0, 0, 1.0), (1, 1, 1.0), (2, 2, 1.0), (3, 3, 1.0)
        };
        var I = SparseMatrix.FromTriplets(4, 4, triplets);
        var x = new VectorN(new[] { 3.0, -1.0, 7.5, 0.2 });

        var y = I.Multiply(x);

        for (int i = 0; i < 4; i++)
            Assert.AreEqual(x[i], y[i], 1e-15);
    }

    [TestMethod]
    public void Multiply_DenseRoundTrip_Matches()
    {
        // Build a 4×4 dense matrix, convert to sparse via triplets, compare Ax
        double[,] dense = {
            { 4, -1, 0, 0 },
            { -1, 4, -1, 0 },
            { 0, -1, 4, -1 },
            { 0, 0, -1, 4 }
        };

        var triplets = new List<(int, int, double)>();
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                if (Math.Abs(dense[i, j]) > 1e-30)
                    triplets.Add((i, j, dense[i, j]));

        var A = SparseMatrix.FromTriplets(4, 4, triplets);
        var x = new VectorN(new[] { 1.0, 2.0, 3.0, 4.0 });

        var ySparse = A.Multiply(x);

        // Dense multiply
        for (int i = 0; i < 4; i++)
        {
            double sum = 0;
            for (int j = 0; j < 4; j++)
                sum += dense[i, j] * x[j];
            Assert.AreEqual(sum, ySparse[i], 1e-12, $"Row {i} mismatch");
        }
    }

    #endregion

    #region Diagonal

    [TestMethod]
    public void Diagonal_ExtractsCorrectValues()
    {
        var triplets = new List<(int, int, double)>
        {
            (0, 0, 10.0), (0, 1, -1.0),
            (1, 0, -1.0), (1, 1, 20.0),
            (2, 2, 30.0)
        };
        var A = SparseMatrix.FromTriplets(3, 3, triplets);

        var d = A.Diagonal();

        Assert.AreEqual(10.0, d[0], 1e-15);
        Assert.AreEqual(20.0, d[1], 1e-15);
        Assert.AreEqual(30.0, d[2], 1e-15);
    }

    [TestMethod]
    public void Diagonal_EmptyRow_ReturnsZero()
    {
        var triplets = new List<(int, int, double)>
        {
            (0, 0, 5.0),
            (2, 2, 7.0)
        };
        var A = SparseMatrix.FromTriplets(3, 3, triplets);

        var d = A.Diagonal();

        Assert.AreEqual(5.0, d[0], 1e-15);
        Assert.AreEqual(0.0, d[1], 1e-15);
        Assert.AreEqual(7.0, d[2], 1e-15);
    }

    #endregion

    #region ApplyDirichlet

    [TestMethod]
    public void ApplyDirichlet_PreservesSymmetry()
    {
        // Symmetric 3×3 tridiagonal
        var triplets = new List<(int, int, double)>
        {
            (0, 0, 4.0), (0, 1, -1.0),
            (1, 0, -1.0), (1, 1, 4.0), (1, 2, -1.0),
            (2, 1, -1.0), (2, 2, 4.0)
        };
        var A = SparseMatrix.FromTriplets(3, 3, triplets);
        var rhs = new VectorN(new[] { 1.0, 2.0, 3.0 });
        var fixedDofs = new Dictionary<int, double> { { 0, 0.5 } };

        var B = A.ApplyDirichlet(fixedDofs, rhs);

        // Check symmetry
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Assert.AreEqual(B.Get(i, j), B.Get(j, i), 1e-15,
                    $"B[{i},{j}] != B[{j},{i}]");

        // Fixed row: diagonal=1, off-diag=0
        Assert.AreEqual(1.0, B.Get(0, 0), 1e-15);
        Assert.AreEqual(0.0, B.Get(0, 1), 1e-15);

        // Fixed column: zeroed out for non-fixed rows
        Assert.AreEqual(0.0, B.Get(1, 0), 1e-15);

        // RHS adjusted: rhs[1] -= A[1,0]*0.5 = 2.0 - (-1.0)*0.5 = 2.5
        Assert.AreEqual(0.5, rhs[0], 1e-15);
        Assert.AreEqual(2.5, rhs[1], 1e-15);
        Assert.AreEqual(3.0, rhs[2], 1e-15);
    }

    [TestMethod]
    public void ApplyDirichlet_MultipleDofs()
    {
        var triplets = new List<(int, int, double)>
        {
            (0, 0, 2.0), (0, 1, -1.0),
            (1, 0, -1.0), (1, 1, 2.0), (1, 2, -1.0),
            (2, 1, -1.0), (2, 2, 2.0)
        };
        var A = SparseMatrix.FromTriplets(3, 3, triplets);
        var rhs = new VectorN(new[] { 0.0, 0.0, 0.0 });
        var fixedDofs = new Dictionary<int, double> { { 0, 1.0 }, { 2, 2.0 } };

        var B = A.ApplyDirichlet(fixedDofs, rhs);

        // Row 0 and Row 2 should be identity rows
        Assert.AreEqual(1.0, B.Get(0, 0), 1e-15);
        Assert.AreEqual(0.0, B.Get(0, 1), 1e-15);
        Assert.AreEqual(1.0, B.Get(2, 2), 1e-15);
        Assert.AreEqual(0.0, B.Get(2, 1), 1e-15);

        // Row 1: columns 0 and 2 zeroed
        Assert.AreEqual(0.0, B.Get(1, 0), 1e-15);
        Assert.AreEqual(2.0, B.Get(1, 1), 1e-15);
        Assert.AreEqual(0.0, B.Get(1, 2), 1e-15);

        // RHS: row 0 = 1.0, row 2 = 2.0
        // row 1: 0 - (-1.0)*1.0 - (-1.0)*2.0 = 3.0
        Assert.AreEqual(1.0, rhs[0], 1e-15);
        Assert.AreEqual(3.0, rhs[1], 1e-15);
        Assert.AreEqual(2.0, rhs[2], 1e-15);
    }

    #endregion

    #region SolvePCG

    [TestMethod]
    public void SolvePCG_TridiagonalSystem_MatchesExact()
    {
        // 4×4 SPD tridiagonal: classic 1D Laplacian
        var triplets = new List<(int, int, double)>
        {
            (0, 0, 2.0), (0, 1, -1.0),
            (1, 0, -1.0), (1, 1, 2.0), (1, 2, -1.0),
            (2, 1, -1.0), (2, 2, 2.0), (2, 3, -1.0),
            (3, 2, -1.0), (3, 3, 2.0)
        };
        var A = SparseMatrix.FromTriplets(4, 4, triplets);
        var b = new VectorN(new[] { 1.0, 0.0, 0.0, 1.0 });

        var x = A.SolvePCG(b);

        // Verify A*x ≈ b
        var Ax = A.Multiply(x);
        for (int i = 0; i < 4; i++)
            Assert.AreEqual(b[i], Ax[i], 1e-8, $"Residual at row {i}");
    }

    [TestMethod]
    public void SolvePCG_WithDirichlet_CorrectSolution()
    {
        // 3×3 system with DOF 0 fixed to 1.0
        var triplets = new List<(int, int, double)>
        {
            (0, 0, 2.0), (0, 1, -1.0),
            (1, 0, -1.0), (1, 1, 2.0), (1, 2, -1.0),
            (2, 1, -1.0), (2, 2, 2.0)
        };
        var A = SparseMatrix.FromTriplets(3, 3, triplets);
        var rhs = new VectorN(new[] { 0.0, 0.0, 1.0 });
        var fixedDofs = new Dictionary<int, double> { { 0, 1.0 } };

        var B = A.ApplyDirichlet(fixedDofs, rhs);
        var x = B.SolvePCG(rhs);

        Assert.AreEqual(1.0, x[0], 1e-8, "Fixed DOF should be 1.0");

        // Verify B*x ≈ rhs
        var Bx = B.Multiply(x);
        for (int i = 0; i < 3; i++)
            Assert.AreEqual(rhs[i], Bx[i], 1e-8, $"Residual at row {i}");
    }

    [TestMethod]
    public void SolvePCG_LargerSystem_Converges()
    {
        // 50×50 1D Laplacian: K = tridiag(-1, 2, -1)
        int n = 50;
        var triplets = new List<(int, int, double)>();
        for (int i = 0; i < n; i++)
        {
            triplets.Add((i, i, 2.0));
            if (i > 0) triplets.Add((i, i - 1, -1.0));
            if (i < n - 1) triplets.Add((i, i + 1, -1.0));
        }
        var A = SparseMatrix.FromTriplets(n, n, triplets);

        // RHS: uniform load
        var bArr = new double[n];
        for (int i = 0; i < n; i++) bArr[i] = 1.0;
        var b = new VectorN(bArr);

        var x = A.SolvePCG(b);
        var Ax = A.Multiply(x);

        double maxResidual = 0;
        for (int i = 0; i < n; i++)
            maxResidual = Math.Max(maxResidual, Math.Abs(Ax[i] - b[i]));

        Assert.IsTrue(maxResidual < 1e-8, $"Max residual {maxResidual} exceeds tolerance");
    }

    [TestMethod]
    public void SolvePCG_ZeroRHS_ReturnsZero()
    {
        var triplets = new List<(int, int, double)>
        {
            (0, 0, 2.0), (0, 1, -1.0),
            (1, 0, -1.0), (1, 1, 2.0)
        };
        var A = SparseMatrix.FromTriplets(2, 2, triplets);
        var b = new VectorN(new[] { 0.0, 0.0 });

        var x = A.SolvePCG(b);

        Assert.AreEqual(0.0, x[0], 1e-15);
        Assert.AreEqual(0.0, x[1], 1e-15);
    }

    #endregion

    #region Get

    [TestMethod]
    public void Get_ExistingEntry_ReturnsValue()
    {
        var triplets = new List<(int, int, double)>
        {
            (0, 0, 5.0), (0, 2, 3.0), (1, 1, 7.0)
        };
        var A = SparseMatrix.FromTriplets(3, 3, triplets);

        Assert.AreEqual(5.0, A.Get(0, 0), 1e-15);
        Assert.AreEqual(3.0, A.Get(0, 2), 1e-15);
        Assert.AreEqual(7.0, A.Get(1, 1), 1e-15);
    }

    [TestMethod]
    public void Get_MissingEntry_ReturnsZero()
    {
        var triplets = new List<(int, int, double)>
        {
            (0, 0, 5.0), (2, 2, 3.0)
        };
        var A = SparseMatrix.FromTriplets(3, 3, triplets);

        Assert.AreEqual(0.0, A.Get(0, 1), 1e-15);
        Assert.AreEqual(0.0, A.Get(1, 1), 1e-15);
        Assert.AreEqual(0.0, A.Get(2, 0), 1e-15);
    }

    #endregion
}
