using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Objects;
using System;

namespace NumericsTests
{
    [TestClass]
    public class MatrixTest
    {
        [TestMethod]
        public void TestMatrixTranspose()
        {
            var matrix = new Matrix(new double[,] { { 1, 3, 7 }, { 5, 2, 9 } });
            var transposematrix = matrix.Transpose();

            Assert.IsTrue(transposematrix.values[0, 0] == 1);
            Assert.IsTrue(transposematrix.values[0, 1] == 5);

            Assert.IsTrue(transposematrix.values[1, 0] == 3);
            Assert.IsTrue(transposematrix.values[1, 1] == 2);

            Assert.IsTrue(transposematrix.values[2, 0] == 7);
            Assert.IsTrue(transposematrix.values[2, 1] == 9);

        }

        [TestMethod]
        public void TestMatrixIdentity()
        {
            var matrix = new Matrix(new double[,] { { 1, 3, 7 }, { 5, 2, 9 } });

            Assert.IsTrue(matrix.identity[0, 0] == 1);
            Assert.IsTrue(matrix.identity[0, 1] == 0);
            Assert.IsTrue(matrix.identity[0, 2] == 0);
            Assert.IsTrue(matrix.identity[1, 0] == 0);
            Assert.IsTrue(matrix.identity[1, 1] == 1);
            Assert.IsTrue(matrix.identity[1, 2] == 0);

        }

        [TestMethod]
        public void TestMatrixDeterminant()
        {
            var matrix = new Matrix(new double[,] { { 1, 5, 2 }, { 0, 3, 7 }, { 2, -1, 4 } });
            var det = matrix.Determinant();

            Assert.IsTrue(det == 77);

        }

        [TestMethod]
        public void TestMatrixAdjugate()
        {
            var matrix = new Matrix(new double[,] { { 1, 5, 2 }, { 0, 3, 7 }, { 2, -1, 4 } });
            var adj = matrix.Adjugate();

            Assert.IsTrue(adj.values[0, 0] == 19);
            Assert.IsTrue(adj.values[0, 1] == -22);
            Assert.IsTrue(adj.values[0, 2] == 29);

            Assert.IsTrue(adj.values[1, 0] == 14);
            Assert.IsTrue(adj.values[1, 1] == 0);
            Assert.IsTrue(adj.values[1, 2] == -7);


            Assert.IsTrue(adj.values[2, 0] == -6);
            Assert.IsTrue(adj.values[2, 1] == 11);
            Assert.IsTrue(adj.values[2, 2] == 3);

        }

        [TestMethod]
        public void TestMatrixInverse()
        {
            var matrix = new Matrix(new double[,] { { 1, 5, 2 }, { 0, 3, 7 }, { 2, -1, 4 } });
            var inv = matrix.Inverse();

            Assert.IsTrue(Math.Round(inv.values[0, 0], 2) == Math.Round((double)19 / 77, 2));
            Assert.IsTrue(Math.Round(inv.values[0, 1], 2) == Math.Round((double)-2 / 7, 2));
            Assert.IsTrue(Math.Round(inv.values[0, 2], 2) == Math.Round((double)29 / 77, 2));

            Assert.IsTrue(Math.Round(inv.values[1, 0], 2) == Math.Round((double)2 / 11, 2));
            Assert.IsTrue(Math.Round(inv.values[1, 1], 2) == Math.Round((double)0, 2));
            Assert.IsTrue(Math.Round(inv.values[1, 2], 2) == Math.Round((double)-1 / 11, 2));


            Assert.IsTrue(Math.Round(inv.values[2, 0], 2) == Math.Round((double)-6 / 77, 2));
            Assert.IsTrue(Math.Round(inv.values[2, 1], 2) == Math.Round((double)1 / 7, 2));
            Assert.IsTrue(Math.Round(inv.values[2, 2], 2) == Math.Round((double)3 / 77, 2));

        }

        [TestMethod]
        public void TestMatrixAddition()
        {
            var a = new Matrix(new double[,] { { 5, 7, 2 }, { -2, 9, 4 } });
            var b = new Matrix(new double[,] { { 1, 3, 7 }, { 5, 2, 9 } });
            var result = a + b;

            Assert.IsTrue(result.values[0, 0] == 6);
            Assert.IsTrue(result.values[0, 1] == 10);
            Assert.IsTrue(result.values[0, 2] == 9);

            Assert.IsTrue(result.values[1, 0] == 3);
            Assert.IsTrue(result.values[1, 1] == 11);
            Assert.IsTrue(result.values[1, 2] == 13);

        }

        [TestMethod]
        public void TestMatrixSubtraction()
        {
            var a = new Matrix(new double[,] { { 5, 7, 2 }, { -2, 9, 4 } });
            var b = new Matrix(new double[,] { { 1, 3, 7 }, { 5, 2, 9 } });
            var result = a - b;

            Assert.IsTrue(result.values[0, 0] == 4);
            Assert.IsTrue(result.values[0, 1] == 4);
            Assert.IsTrue(result.values[0, 2] == -5);

            Assert.IsTrue(result.values[1, 0] == -7);
            Assert.IsTrue(result.values[1, 1] == 7);
            Assert.IsTrue(result.values[1, 2] == -5);

        }

        [TestMethod]
        public void TestMatrixMultiplier()
        {
            var a = new Matrix(new double[,] { { 4, 0 }, { -1, 5 } });
            var b = 3;
            var result = b * a;

            Assert.IsTrue(result.values[0, 0] == 12);
            Assert.IsTrue(result.values[0, 1] == 0);


            Assert.IsTrue(result.values[1, 0] == -3);
            Assert.IsTrue(result.values[1, 1] == 15);

        }

        [TestMethod]
        public void TestMatrixMultiplication()
        {
            var a = new Matrix(new double[,] { { 2, -1 }, { 3, 0 }, { 0, 4 } });
            var b = new Matrix(new double[,] { { 5, 1, -3 }, { -1, 0, 1 } });
            var result = a * b;

            Assert.IsTrue(result.values[0, 0] == 11);
            Assert.IsTrue(result.values[0, 1] == 2);
            Assert.IsTrue(result.values[0, 2] == -7);

            Assert.IsTrue(result.values[1, 0] == 15);
            Assert.IsTrue(result.values[1, 1] == 3);
            Assert.IsTrue(result.values[1, 2] == -9);

            Assert.IsTrue(result.values[2, 0] == -4);
            Assert.IsTrue(result.values[2, 1] == 0);
            Assert.IsTrue(result.values[2, 2] == 4);

        }
        [TestMethod]
        public void TestMatrixVectorMultiplication()
        {
            var a = new Matrix(new double[,] { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } });
            var b = new Vector(2, 1, 3);
            var result = a * b;

            Assert.IsTrue(result.x == 13);
            Assert.IsTrue(result.y == 31);
            Assert.IsTrue(result.z == 49);


        }
    }
}
