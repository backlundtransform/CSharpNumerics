using CSharpNumerics.Numerics.LinearAlgebra;
using CSharpNumerics.Numerics.Objects;

namespace NumericsTests
{
    [TestClass]
    public class QrDecompositionTest
    {
        private const double Tolerance = 1e-10;

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
        public void TestRoundTripSquare()
        {
            var matrix = new Matrix(new double[,] { { 1, 5, 2 }, { 0, 3, 7 }, { 2, -1, 4 } });
            var qr = matrix.Qr();

            AssertMatricesEqual(matrix, qr.Q * qr.R);
        }

        [TestMethod]
        public void TestRoundTripOverdetermined()
        {
            var matrix = new Matrix(new double[,] { { 1, 1 }, { 1, 2 }, { 1, 3 }, { 1, 4 } });
            var qr = matrix.Qr();

            AssertMatricesEqual(matrix, qr.Q * qr.R);
        }

        [TestMethod]
        public void TestQHasOrthonormalColumns()
        {
            var matrix = new Matrix(new double[,] { { 1, 1 }, { 1, 2 }, { 1, 3 }, { 1, 4 } });
            var q = matrix.Qr().Q;
            var identity = q.Transpose() * q;

            for (var i = 0; i < 2; i++)
            {
                for (var j = 0; j < 2; j++)
                {
                    Assert.IsTrue(Math.Abs(identity.values[i, j] - (i == j ? 1 : 0)) < Tolerance);
                }
            }
        }

        [TestMethod]
        public void TestRIsUpperTriangular()
        {
            var matrix = new Matrix(new double[,] { { 1, 5, 2 }, { 0, 3, 7 }, { 2, -1, 4 } });
            var r = matrix.Qr().R;

            for (var i = 1; i < 3; i++)
            {
                for (var j = 0; j < i; j++)
                {
                    Assert.IsTrue(r.values[i, j] == 0);
                }
            }
        }

        [TestMethod]
        public void TestLeastSquaresSolution()
        {
            var matrix = new Matrix(new double[,] { { 1, 1 }, { 1, 2 }, { 1, 3 } });
            var solution = matrix.Qr().Solve(new VectorN(new double[] { 6, 0, 0 }));

            Assert.IsTrue(Math.Abs(solution[0] - 8) < Tolerance);
            Assert.IsTrue(Math.Abs(solution[1] - (-3)) < Tolerance);
        }

        [TestMethod]
        public void TestSquareSolveMatchesLu()
        {
            var matrix = new Matrix(new double[,] { { 2, 1, -1 }, { -3, -1, 2 }, { -2, 1, 2 } });
            var b = new VectorN(new double[] { 8, -11, -3 });

            var qrSolution = matrix.Qr().Solve(b);
            var luSolution = matrix.Lu().Solve(b);

            for (var i = 0; i < 3; i++)
            {
                Assert.IsTrue(Math.Abs(qrSolution[i] - luSolution[i]) < Tolerance);
            }
        }

        [TestMethod]
        public void TestSolveList()
        {
            var matrix = new Matrix(new double[,] { { 1, 1 }, { 1, 2 }, { 1, 3 } });
            var solution = matrix.Qr().Solve(new List<double> { 6, 0, 0 });

            Assert.IsTrue(Math.Abs(solution[0] - 8) < Tolerance);
            Assert.IsTrue(Math.Abs(solution[1] - (-3)) < Tolerance);
        }

        [TestMethod]
        public void TestSolveMultipleRightHandSides()
        {
            var matrix = new Matrix(new double[,] { { 2, 1, -1 }, { -3, -1, 2 }, { -2, 1, 2 } });
            var rhs = new Matrix(new double[,] { { 8, 1 }, { -11, 0 }, { -3, 2 } });

            var x = matrix.Qr().Solve(rhs);

            AssertMatricesEqual(rhs, matrix * x, 1e-9);
        }

        [TestMethod]
        public void TestRankDeficientIsDetected()
        {
            var matrix = new Matrix(new double[,] { { 1, 2 }, { 2, 4 }, { 3, 6 } });
            var qr = matrix.Qr();

            Assert.IsFalse(qr.IsFullRank);
            Assert.ThrowsException<Exception>(() => qr.Solve(new VectorN(new double[] { 1, 1, 1 })));
        }

        [TestMethod]
        public void TestUnderdeterminedThrows()
        {
            var matrix = new Matrix(new double[,] { { 1, 2, 3 }, { 4, 5, 6 } });

            Assert.ThrowsException<Exception>(() => matrix.Qr());
        }
    }
}
