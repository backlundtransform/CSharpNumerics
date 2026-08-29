using CSharpNumerics.Numerics.LinearAlgebra;
using CSharpNumerics.Numerics.Objects;

namespace NumericsTests
{
    [TestClass]
    public class CholeskyDecompositionTest
    {
        private const double Tolerance = 1e-10;

        private static readonly double[,] SpdValues =
        {
            { 4, 12, -16 },
            { 12, 37, -43 },
            { -16, -43, 98 }
        };

        private static void AssertMatricesEqual(Matrix expected, Matrix actual, double tolerance = Tolerance)
        {
            Assert.AreEqual(expected.rowLength, actual.rowLength);
            Assert.AreEqual(expected.columnLength, actual.columnLength);

            for (var i = 0; i < expected.rowLength; i++)
            {
                for (var j = 0; j < expected.columnLength; j++)
                {
                    Assert.IsTrue(Math.Abs(expected.values[i, j] - actual.values[i, j]) < tolerance,
                        $"Mismatch at ({i},{j}): expected {expected.values[i, j]}, actual {actual.values[i, j]}");
                }
            }
        }

        [TestMethod]
        public void TestKnownFactorization()
        {
            var matrix = new Matrix(SpdValues);
            var cholesky = matrix.Cholesky();

            Assert.IsTrue(cholesky.IsPositiveDefinite);

            var expectedLower = new Matrix(new double[,]
            {
                { 2, 0, 0 },
                { 6, 1, 0 },
                { -8, 5, 3 }
            });

            AssertMatricesEqual(expectedLower, cholesky.Lower);
        }

        [TestMethod]
        public void TestRoundTrip()
        {
            var matrix = new Matrix(SpdValues);
            var lower = matrix.Cholesky().Lower;

            AssertMatricesEqual(matrix, lower * lower.Transpose());
        }

        [TestMethod]
        public void TestSolveKnownSystem()
        {
            var matrix = new Matrix(SpdValues);
            var x = new VectorN(new double[] { 1, 2, 3 });
            var b = matrix * x;

            var solution = matrix.Cholesky().Solve(b);

            Assert.IsTrue(Math.Abs(solution[0] - 1) < Tolerance);
            Assert.IsTrue(Math.Abs(solution[1] - 2) < Tolerance);
            Assert.IsTrue(Math.Abs(solution[2] - 3) < Tolerance);
        }

        [TestMethod]
        public void TestSolveList()
        {
            var matrix = new Matrix(SpdValues);
            var b = matrix * new List<double> { 1, 2, 3 };

            var solution = matrix.Cholesky().Solve(b);

            Assert.IsTrue(Math.Abs(solution[0] - 1) < Tolerance);
            Assert.IsTrue(Math.Abs(solution[1] - 2) < Tolerance);
            Assert.IsTrue(Math.Abs(solution[2] - 3) < Tolerance);
        }

        [TestMethod]
        public void TestSolveMatchesLu()
        {
            var matrix = new Matrix(SpdValues);
            var b = new VectorN(new double[] { 1, -2, 5 });

            var choleskySolution = matrix.Cholesky().Solve(b);
            var luSolution = matrix.Lu().Solve(b);

            for (var i = 0; i < 3; i++)
            {
                Assert.IsTrue(Math.Abs(choleskySolution[i] - luSolution[i]) < Tolerance);
            }
        }

        [TestMethod]
        public void TestDeterminant()
        {
            var matrix = new Matrix(SpdValues);

            Assert.IsTrue(Math.Abs(matrix.Cholesky().Determinant() - 36) < Tolerance);
        }

        [TestMethod]
        public void TestInverseRoundTrip()
        {
            var matrix = new Matrix(SpdValues);
            var inverse = matrix.Cholesky().Inverse();

            AssertMatricesEqual(new Matrix(matrix.identity), matrix * inverse, 1e-9);
        }

        [TestMethod]
        public void TestSymmetricIndefiniteIsRejected()
        {
            var matrix = new Matrix(new double[,] { { 1, 2 }, { 2, 1 } });
            var cholesky = matrix.Cholesky();

            Assert.IsFalse(cholesky.IsPositiveDefinite);
            Assert.ThrowsException<Exception>(() => cholesky.Solve(new VectorN(new double[] { 1, 1 })));
        }

        [TestMethod]
        public void TestNonSymmetricIsRejected()
        {
            var matrix = new Matrix(new double[,] { { 4, 1 }, { 2, 3 } });
            var cholesky = matrix.Cholesky();

            Assert.IsFalse(cholesky.IsPositiveDefinite);
            Assert.ThrowsException<Exception>(() => cholesky.Lower);
        }

        [TestMethod]
        public void TestNonSquareMatrixThrows()
        {
            var matrix = new Matrix(new double[,] { { 1, 2, 3 }, { 4, 5, 6 } });

            Assert.ThrowsException<Exception>(() => matrix.Cholesky());
        }
    }
}
