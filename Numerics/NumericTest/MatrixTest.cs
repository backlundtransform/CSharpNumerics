
using CSharpNumerics.Numerics.Objects;
using Xunit.Sdk;

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


        [TestMethod]
        public void TestPascalMatrix()
        {
            var pascalMatrix = new Matrix(new double[6, 6]).Pascal();
            Assert.IsTrue(pascalMatrix.values[1, 1] == 1);
            Assert.IsTrue(pascalMatrix.values[2, 1] == 2);
            Assert.IsTrue(pascalMatrix.values[3, 1] == 3);
            Assert.IsTrue(pascalMatrix.values[3, 2] == 3);
            Assert.IsTrue(pascalMatrix.values[4, 1] == 4);
            Assert.IsTrue(pascalMatrix.values[4, 2] == 6);
            Assert.IsTrue(pascalMatrix.values[5, 1] == 5);
            Assert.IsTrue(pascalMatrix.values[5, 2] == 10);

        }


        [TestMethod]
        public void TestTensorRank()
        {
            var tensor1D = new Tensor(6);         // shape: [6]
            var tensor2D = new Tensor(2, 3);      // shape: [2,3]
            var tensor3D = new Tensor(2, 3, 3);   // shape: [2,3,3]
            var tensor4D = new Tensor(2, 3, 3, 3);// shape: [2,3,3,3]

            Assert.AreEqual(1, tensor1D.Dimension);
            Assert.AreEqual(2, tensor2D.Dimension);
            Assert.AreEqual(3, tensor3D.Dimension);
            Assert.AreEqual(4, tensor4D.Dimension);
        }

        [TestMethod]
        public void TestTensorShape()
        {
            var tensor1D = new Tensor(6);
            var tensor3D = new Tensor(2, 3, 3);

            Assert.AreEqual(6, tensor1D.Shape[0]);
            Assert.AreEqual(2, tensor3D.Shape[0]);
            Assert.AreEqual(3, tensor3D.Shape[1]);
            Assert.AreEqual(3, tensor3D.Shape[2]);
        }

        [TestMethod]
        public void TestTensorAddition()
        {
            var tensor = new Tensor(2, 3);
            tensor.Fill(1.0);
            var result = tensor + tensor;

            foreach (var val in result.Values)
                Assert.AreEqual(2.0, val);
        }

        [TestMethod]
        public void TestTensorSubstraction()
        {
            var tensor = new Tensor(2, 3);
            tensor.Fill(1.0);
            var result = tensor - tensor;

            foreach (var val in result.Values)
                Assert.AreEqual(0.0, val);
        }

        [TestMethod]
        public void TestTensorMultiplication()
        {
            var tensor = new Tensor(2, 3);
            for (int i = 0; i < tensor.Values.Length; i++)
                tensor.Values[i] = i + 1; // 1,2,3,4,5,6

            var result = tensor * tensor;
            for (int i = 0; i < tensor.Values.Length; i++)
                Assert.AreEqual((i + 1) * (i + 1), result.Values[i]);
        }

        [TestMethod]
        public void TestTensorDivision()
        {
            var tensor = new Tensor(2, 3);
            for (int i = 0; i < tensor.Values.Length; i++)
                tensor.Values[i] = i + 1;

            var result = tensor / tensor;
            foreach (var val in result.Values)
                Assert.AreEqual(1.0, val);
        }


        [TestMethod]
        public void TestTensorAdditionFail()
        {
            var a = new Tensor(2, 3);
            var b = new Tensor(3, 2); 
            Assert.ThrowsException<ArgumentException>(() => { var c = a + b; });
        }


        [TestMethod]
        public void TestTensorDot()
        {
            var tensorX = new Tensor(2);
            tensorX.Values[0] = 1;
            tensorX.Values[1] = 2;

            var tensorY = new Tensor(2);
            tensorY.Values[0] = 3;
            tensorY.Values[1] = 4;

            double result = tensorX.Dot(tensorY);
            Assert.AreEqual(1 * 3 + 2 * 4, result); 
        }


        [TestMethod]
        public void TestMatixSlice()
        {
            var mat = new Matrix(new double[,] {
    {1,2,3},
    {4,5,6},
    {7,8,9}
});

        
            VectorN row1 = mat.RowSlice(1); // [4,5,6]

            Assert.IsTrue(row1.Values[0] == 4);
            Assert.IsTrue(row1.Values[1] == 5);
            Assert.IsTrue(row1.Values[2] == 6);


            // Hämta kolumn 2
            VectorN col2 = mat.ColumnSlice(2); // [3,6,9]
            Assert.IsTrue(col2.Values[0] == 3);
            Assert.IsTrue(col2.Values[1] == 6);
            Assert.IsTrue(col2.Values[2] == 9);


            // Hämta submatrix 0-2 rader, 1-3 kolumner
            Matrix sub = mat.Slice(0, 2, 1, 3);
            Assert.IsTrue(sub.values[0, 0] == 2);
            Assert.IsTrue(sub.values[0, 1] == 3);
            Assert.IsTrue(sub.values[1, 0] == 5);
            Assert.IsTrue(sub.values[1, 1] == 6);

        }
    }
}
