using CSharpNumerics.Numerics;
using CSharpNumerics.Numerics.LinearAlgebra;
using CSharpNumerics.Numerics.LinearAlgebra.Decompositions;
using CSharpNumerics.Numerics.Objects;

namespace NumericsTests
{
    [TestClass]
    public class LuDecompositionTest
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
        public void TestKnownFactorization()
        {
            var matrix = new Matrix(new double[,] { { 4, 3 }, { 6, 3 } });
            var lu = matrix.Lu();

            var lower = lu.Lower;
            var upper = lu.Upper;

            Assert.IsTrue(Math.Abs(lower.values[0, 0] - 1) < Tolerance);
            Assert.IsTrue(Math.Abs(lower.values[1, 0] - 2.0 / 3.0) < Tolerance);
            Assert.IsTrue(Math.Abs(lower.values[1, 1] - 1) < Tolerance);

            Assert.IsTrue(Math.Abs(upper.values[0, 0] - 6) < Tolerance);
            Assert.IsTrue(Math.Abs(upper.values[0, 1] - 3) < Tolerance);
            Assert.IsTrue(Math.Abs(upper.values[1, 1] - 1) < Tolerance);
        }

        [TestMethod]
        public void TestRoundTrip()
        {
            var matrix = new Matrix(new double[,]
            {
                { 2, 1, 1, 0 },
                { 4, 3, 3, 1 },
                { 8, 7, 9, 5 },
                { 6, 7, 9, 8 }
            });

            var lu = matrix.Lu();
            AssertMatricesEqual(lu.PermutationMatrix * matrix, lu.Lower * lu.Upper);
        }

        [TestMethod]
        public void TestLowerIsUnitLowerAndUpperIsUpper()
        {
            var matrix = new Matrix(new double[,] { { 1, 5, 2 }, { 0, 3, 7 }, { 2, -1, 4 } });
            var lu = matrix.Lu();

            for (var i = 0; i < 3; i++)
            {
                Assert.IsTrue(lu.Lower.values[i, i] == 1);
                for (var j = i + 1; j < 3; j++)
                {
                    Assert.IsTrue(lu.Lower.values[i, j] == 0);
                    Assert.IsTrue(lu.Upper.values[j, i] == 0);
                }
            }
        }

        [TestMethod]
        public void TestSolveKnownSystem()
        {
            var matrix = new Matrix(new double[,] { { 2, 1, -1 }, { -3, -1, 2 }, { -2, 1, 2 } });
            var solution = matrix.Lu().Solve(new VectorN(new double[] { 8, -11, -3 }));

            Assert.IsTrue(Math.Abs(solution[0] - 2) < Tolerance);
            Assert.IsTrue(Math.Abs(solution[1] - 3) < Tolerance);
            Assert.IsTrue(Math.Abs(solution[2] - (-1)) < Tolerance);
        }

        [TestMethod]
        public void TestSolveList()
        {
            var matrix = new Matrix(new double[,] { { 2, 1, -1 }, { -3, -1, 2 }, { -2, 1, 2 } });
            var solution = matrix.Lu().Solve(new List<double> { 8, -11, -3 });

            Assert.IsTrue(Math.Abs(solution[0] - 2) < Tolerance);
            Assert.IsTrue(Math.Abs(solution[1] - 3) < Tolerance);
            Assert.IsTrue(Math.Abs(solution[2] - (-1)) < Tolerance);
        }

        [TestMethod]
        public void TestSolveReusesFactorization()
        {
            var matrix = new Matrix(new double[,] { { 4, 3 }, { 6, 3 } });
            var lu = matrix.Lu();

            var x1 = lu.Solve(new VectorN(new double[] { 10, 12 }));
            var x2 = lu.Solve(new VectorN(new double[] { 7, 9 }));

            var b1 = matrix * x1;
            var b2 = matrix * x2;

            Assert.IsTrue(Math.Abs(b1[0] - 10) < Tolerance);
            Assert.IsTrue(Math.Abs(b1[1] - 12) < Tolerance);
            Assert.IsTrue(Math.Abs(b2[0] - 7) < Tolerance);
            Assert.IsTrue(Math.Abs(b2[1] - 9) < Tolerance);
        }

        [TestMethod]
        public void TestSolveMultipleRightHandSides()
        {
            var matrix = new Matrix(new double[,] { { 2, 1, -1 }, { -3, -1, 2 }, { -2, 1, 2 } });
            var rhs = new Matrix(new double[,] { { 8, 1 }, { -11, 0 }, { -3, 2 } });

            var x = matrix.Lu().Solve(rhs);

            AssertMatricesEqual(rhs, matrix * x);
        }

        [TestMethod]
        public void TestDeterminantMatchesCofactorExpansion()
        {
            var matrix = new Matrix(new double[,] { { 1, 5, 2 }, { 0, 3, 7 }, { 2, -1, 4 } });

            Assert.IsTrue(Math.Abs(matrix.Lu().Determinant() - 77) < Tolerance);
            Assert.IsTrue(Math.Abs(matrix.Lu().Determinant() - matrix.Determinant()) < Tolerance);
        }

        [TestMethod]
        public void TestInverseRoundTrip()
        {
            var matrix = new Matrix(new double[,] { { 1, 5, 2 }, { 0, 3, 7 }, { 2, -1, 4 } });
            var inverse = matrix.Lu().Inverse();

            AssertMatricesEqual(new Matrix(matrix.identity), matrix * inverse);
            AssertMatricesEqual(new Matrix(matrix.identity), inverse * matrix);
        }

        [TestMethod]
        public void TestMatrixInverseUsesLuAndMatchesAdjugate()
        {
            var matrix = new Matrix(new double[,] { { 1, 5, 2 }, { 0, 3, 7 }, { 2, -1, 4 } });
            var inverse = matrix.Inverse();
            var expected = matrix.Adjugate() / matrix.Determinant();

            AssertMatricesEqual(expected, inverse);
        }

        [TestMethod]
        public void TestSingularMatrixIsDetected()
        {
            var matrix = new Matrix(new double[,] { { 1, 2 }, { 2, 4 } });
            var lu = matrix.Lu();

            Assert.IsTrue(lu.IsSingular);
            Assert.ThrowsException<Exception>(() => lu.Solve(new VectorN(new double[] { 1, 1 })));
            Assert.ThrowsException<Exception>(() => matrix.Inverse());
        }

        [TestMethod]
        public void TestNonSquareMatrixThrows()
        {
            var matrix = new Matrix(new double[,] { { 1, 2, 3 }, { 4, 5, 6 } });

            Assert.ThrowsException<Exception>(() => matrix.Lu());
        }

        [TestMethod]
        public void TestPivotingHandlesZeroLeadingElement()
        {
            var matrix = new Matrix(new double[,] { { 0, 1 }, { 1, 0 } });
            var lu = matrix.Lu();

            Assert.IsFalse(lu.IsSingular);

            var solution = lu.Solve(new VectorN(new double[] { 3, 5 }));

            Assert.IsTrue(Math.Abs(solution[0] - 5) < Tolerance);
            Assert.IsTrue(Math.Abs(solution[1] - 3) < Tolerance);
        }

        [TestMethod]
        public void TestLinearSystemSolverStillWorks()
        {
            var matrix = new Matrix(new double[,] { { 2, 1, -1 }, { -3, -1, 2 }, { -2, 1, 2 } });
            var solution = matrix.LinearSystemSolver(new VectorN(new double[] { 8, -11, -3 }));

            Assert.IsTrue(Math.Abs(solution[0] - 2) < Tolerance);
            Assert.IsTrue(Math.Abs(solution[1] - 3) < Tolerance);
            Assert.IsTrue(Math.Abs(solution[2] - (-1)) < Tolerance);
        }

        [TestMethod]
        public void TestLinearSystemSolverVector()
        {
            var matrix = new Matrix(new double[,] { { 2, 1, -1 }, { -3, -1, 2 }, { -2, 1, 2 } });
            var solution = matrix.LinearSystemSolver(new Vector(8, -11, -3));

            Assert.IsTrue(Math.Abs(solution.x - 2) < Tolerance);
            Assert.IsTrue(Math.Abs(solution.y - 3) < Tolerance);
            Assert.IsTrue(Math.Abs(solution.z - (-1)) < Tolerance);
        }
    }
}
