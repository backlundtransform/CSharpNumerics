using CSharpNumerics.Numerics.LinearAlgebra;
using CSharpNumerics.Numerics.LinearAlgebra.Decompositions;
using CSharpNumerics.Numerics.Objects;

namespace NumericsTests
{
    [TestClass]
    public class EigenDecompositionTest
    {
        private const double Tolerance = 1e-9;

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
        public void TestSymmetricKnownEigenvalues()
        {
            var matrix = new Matrix(new double[,] { { 2, 1 }, { 1, 2 } });
            var eigen = matrix.Eigen();

            Assert.IsTrue(eigen.IsSymmetric);
            Assert.IsTrue(Math.Abs(eigen.RealEigenvalues[0] - 1) < Tolerance);
            Assert.IsTrue(Math.Abs(eigen.RealEigenvalues[1] - 3) < Tolerance);
            Assert.IsTrue(eigen.ImaginaryEigenvalues[0] == 0);
            Assert.IsTrue(eigen.ImaginaryEigenvalues[1] == 0);
        }

        [TestMethod]
        public void TestSymmetricEigenvaluesAreAscending()
        {
            var matrix = new Matrix(new double[,]
            {
                { 4, 1, -2, 2 },
                { 1, 2, 0, 1 },
                { -2, 0, 3, -2 },
                { 2, 1, -2, -1 }
            });
            var eigen = matrix.Eigen();

            for (var i = 1; i < 4; i++)
            {
                Assert.IsTrue(eigen.RealEigenvalues[i] >= eigen.RealEigenvalues[i - 1]);
            }
        }

        [TestMethod]
        public void TestSymmetricEigenpairResiduals()
        {
            var matrix = new Matrix(new double[,]
            {
                { 4, 1, -2, 2 },
                { 1, 2, 0, 1 },
                { -2, 0, 3, -2 },
                { 2, 1, -2, -1 }
            });
            var eigen = matrix.Eigen();
            var v = eigen.EigenVectors;

            for (var k = 0; k < 4; k++)
            {
                var vector = v.ColumnSlice(k);
                var av = matrix * vector;

                for (var i = 0; i < 4; i++)
                {
                    Assert.IsTrue(Math.Abs(av[i] - eigen.RealEigenvalues[k] * vector[i]) < Tolerance,
                        $"Residual too large for eigenpair {k}, component {i}");
                }
            }
        }

        [TestMethod]
        public void TestSymmetricEigenvectorsAreOrthonormal()
        {
            var matrix = new Matrix(new double[,]
            {
                { 4, 1, -2, 2 },
                { 1, 2, 0, 1 },
                { -2, 0, 3, -2 },
                { 2, 1, -2, -1 }
            });
            var v = matrix.Eigen().EigenVectors;
            var identity = v.Transpose() * v;

            for (var i = 0; i < 4; i++)
            {
                for (var j = 0; j < 4; j++)
                {
                    Assert.IsTrue(Math.Abs(identity.values[i, j] - (i == j ? 1 : 0)) < Tolerance);
                }
            }
        }

        [TestMethod]
        public void TestSymmetricRoundTrip()
        {
            var matrix = new Matrix(new double[,] { { 4, 12, -16 }, { 12, 37, -43 }, { -16, -43, 98 } });
            var eigen = matrix.Eigen();

            AssertMatricesEqual(matrix * eigen.EigenVectors, eigen.EigenVectors * eigen.DiagonalMatrix, 1e-8);
        }

        [TestMethod]
        public void TestNonSymmetricRealEigenvalues()
        {
            var matrix = new Matrix(new double[,] { { 4, 1 }, { 2, 3 } });
            var eigen = matrix.Eigen();

            Assert.IsFalse(eigen.IsSymmetric);

            var values = new List<double>(eigen.RealEigenvalues);
            values.Sort();

            Assert.IsTrue(Math.Abs(values[0] - 2) < Tolerance);
            Assert.IsTrue(Math.Abs(values[1] - 5) < Tolerance);
            Assert.IsTrue(Math.Abs(eigen.ImaginaryEigenvalues[0]) < Tolerance);
            Assert.IsTrue(Math.Abs(eigen.ImaginaryEigenvalues[1]) < Tolerance);
        }

        [TestMethod]
        public void TestNonSymmetricRoundTrip()
        {
            var matrix = new Matrix(new double[,] { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 10 } });
            var eigen = matrix.Eigen();

            AssertMatricesEqual(matrix * eigen.EigenVectors, eigen.EigenVectors * eigen.DiagonalMatrix, 1e-8);
        }

        [TestMethod]
        public void TestComplexConjugatePair()
        {
            var matrix = new Matrix(new double[,] { { 0, -1 }, { 1, 0 } });
            var eigen = matrix.Eigen();

            Assert.IsTrue(Math.Abs(eigen.RealEigenvalues[0]) < Tolerance);
            Assert.IsTrue(Math.Abs(eigen.RealEigenvalues[1]) < Tolerance);
            Assert.IsTrue(Math.Abs(Math.Abs(eigen.ImaginaryEigenvalues[0]) - 1) < Tolerance);
            Assert.IsTrue(Math.Abs(eigen.ImaginaryEigenvalues[0] + eigen.ImaginaryEigenvalues[1]) < Tolerance);

            AssertMatricesEqual(matrix * eigen.EigenVectors, eigen.EigenVectors * eigen.DiagonalMatrix);
        }

        [TestMethod]
        public void TestDeterminantMatchesEigenvalueProduct()
        {
            var matrix = new Matrix(new double[,] { { 4, 12, -16 }, { 12, 37, -43 }, { -16, -43, 98 } });
            var eigen = matrix.Eigen();

            var product = 1.0;
            foreach (var value in eigen.RealEigenvalues)
            {
                product *= value;
            }

            Assert.IsTrue(Math.Abs(product - 36) < 1e-6);
        }

        [TestMethod]
        public void TestDiagonalMatrixIsTrivial()
        {
            var matrix = new Matrix(new double[,] { { 3, 0 }, { 0, 7 } });
            var eigen = matrix.Eigen();

            Assert.IsTrue(Math.Abs(eigen.RealEigenvalues[0] - 3) < Tolerance);
            Assert.IsTrue(Math.Abs(eigen.RealEigenvalues[1] - 7) < Tolerance);
        }

        [TestMethod]
        public void TestNonSquareMatrixThrows()
        {
            var matrix = new Matrix(new double[,] { { 1, 2, 3 }, { 4, 5, 6 } });

            Assert.ThrowsException<Exception>(() => matrix.Eigen());
        }

        [TestMethod]
        public void TestLargeSymmetricTridiagonal()
        {
            // Discrete 1D Laplacian: eigenvalues 2 - 2cos(k*pi/(n+1)) are known analytically.
            const int n = 20;
            var values = new double[n, n];
            for (var i = 0; i < n; i++)
            {
                values[i, i] = 2;
                if (i > 0) values[i, i - 1] = -1;
                if (i < n - 1) values[i, i + 1] = -1;
            }

            var eigen = new EigenDecomposition(new Matrix(values));

            for (var k = 0; k < n; k++)
            {
                var expected = 2 - 2 * Math.Cos((k + 1) * Math.PI / (n + 1));
                Assert.IsTrue(Math.Abs(eigen.RealEigenvalues[k] - expected) < 1e-8,
                    $"Eigenvalue {k}: expected {expected}, actual {eigen.RealEigenvalues[k]}");
            }
        }
    }
}
