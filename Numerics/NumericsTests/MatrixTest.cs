using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Objects;
using System;
using System.Collections.Generic;

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
            Assert.IsTrue(pascalMatrix.values[2, 1] ==2);
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
        
            var tensor1D = new Tensor(new double[] { 1, 3, 7 , 5, 2, 9  });

            var tensor2D = new Tensor(new double[,] { { 1, 3, 7 }, { 5, 2, 9 } });

            var tensor3D = new Tensor(new double[,,] { { { 1, 3, 7 }, { 1, 3, 7 }, { 1, 3, 7 } }, { { 1, 3, 7 }, { 1, 3, 7 }, { 1, 3, 7 } } });

            var tensor4D = new Tensor(new double[,,,] { { { { 1, 3, 7 }, { 1, 3, 7 }, { 1, 3, 7 } }, { { 1, 3, 7 }, { 1, 3, 7 }, { 1, 3, 7 } }, { { 1, 3, 7 }, { 1, 3, 7 }, { 1, 3, 7 } } }, { { { 1, 3, 7 }, { 1, 3, 7 }, { 1, 3, 7 } }, { { 1, 3, 7 }, { 1, 3, 7 }, { 1, 3, 7 } }, { { 1, 3, 7 }, { 1, 3, 7 }, { 1, 3, 7 } } } });

           Assert.IsTrue(tensor1D.dimension ==1);
            Assert.IsTrue(tensor2D.dimension == 2);
            Assert.IsTrue(tensor3D.dimension == 3);
            Assert.IsTrue(tensor4D.dimension == 4);
        }


        [TestMethod]
        public void TestTensorShape()
        {


           var tensor1D = new Tensor(new double[] { 1, 3, 7, 5, 2, 9 });
           var tensor3D = new Tensor(new double[,,] { { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } }, { { 10, 11, 12 }, { 13, 14, 15 }, { 16, 17, 18 } } });

            Assert.IsTrue(tensor1D.shape[0] == 6);

            Assert.IsTrue(tensor3D.shape[0] == 2);

            Assert.IsTrue(tensor3D.shape[1] == 3);

            Assert.IsTrue(tensor3D.shape[2] == 3);


        }

        [TestMethod]
        public void TestTensorAddition()
        {


            var tensor3D = new Tensor(new double[,,] { { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } },
                { { 11, 12, 13 }, { 14, 15, 16 }, { 17, 18, 19 } },  { { 21, 22, 23 }, { 24, 25, 26}, { 27, 28, 29 } }  });

            var result = tensor3D + tensor3D;

            CollectionAssert.AreEqual(result.values,new double[,,] { { { 2, 4, 6 }, { 8, 10, 12}, { 14, 16, 18 } },
                { { 22, 24, 26}, { 28, 30, 32 }, { 34, 36, 38 } },  { { 42, 44, 46 }, { 48, 50, 52}, { 54, 56, 58 } } });

        }

        [TestMethod]
        public void TestTensorSubstraction()
        {


            var tensor3D = new Tensor(new double[,,] { { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } },
                { { 11, 12, 13 }, { 14, 15, 16 }, { 17, 18, 19 } },  { { 21, 22, 23 }, { 24, 25, 26}, { 27, 28, 29 } }  });

            var result = tensor3D - tensor3D;

            CollectionAssert.AreEqual(result.values, new double[,,] { { { 0, 0, 0 }, { 0, 0, 0}, { 0, 0, 0 } },
                { { 0,0, 0}, { 0, 0,0 }, { 0, 0, 0 } },  { { 0, 0, 0}, { 0, 0,0}, { 0, 0, 0 } } });

        }


        [TestMethod]
        public void TestTensorDivision()
        {

            var tensor3D = new Tensor(new double[,,] { { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } },
                { { 11, 12, 13 }, { 14, 15, 16 }, { 17, 18, 19 } },  { { 21, 22, 23 }, { 24, 25, 26}, { 27, 28, 29 } }  });

            var result = tensor3D / tensor3D;

            CollectionAssert.AreEqual(result.values, new double[,,] { { { 1, 1, 1 }, { 1, 1, 1}, { 1, 1, 1 } },
                { { 1,1, 1}, { 1, 1,1 }, { 1, 1, 1 } },  { { 1, 1, 1}, { 1, 1, 1}, { 1, 1, 1} } });

        }


        [TestMethod]
        public void TestTensorMultiplication()
        {


            var tensor3D = new Tensor(new double[,,] { { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } },
                { { 11, 12, 13 }, { 14, 15, 16 }, { 17, 18, 19 } },  { { 21, 22, 23 }, { 24, 25, 26}, { 27, 28, 29 } }  });

            var result = tensor3D *tensor3D;

            CollectionAssert.AreEqual(result.values, new double[,,] { { { 1*1, 2*2, 3*3 }, { 4*4, 5*5, 6*6 }, { 7*7, 8*8, 9*9 } },
                { { 11*11, 12*12, 13*13 }, { 14*14, 15*15, 16*16 }, { 17*17, 18*18, 19*19 } },  { { 21*21, 22*22, 23*23 }, { 24*24, 25*25, 26*26}, { 27*27, 28*28, 29*29 } }  });

        }


        [TestMethod]
        public void TestTensorAdditionFail()
        { 

            try
            {
                var tensor3D = new Tensor(new string[,] { { "WDWD" } });

                var result = tensor3D + tensor3D;
                Assert.Fail("no exception thrown");
            }
            catch (Exception ex)
            {
                Assert.IsTrue(ex is ArgumentException);
            }

        }


        [TestMethod]
        public void TestTensorDot()
        {
  
           var tensor1Dx = new Tensor(new double[] { 1,2});
           var tensor1Dy = new Tensor(new double[] { 3, 4 });
           var result = tensor1Dx.TensorDot(tensor1Dy);
          
            Assert.AreEqual(result.values.GetValue(0).ToString(),  new double[2] { 3, 4}.ToString());
            Assert.AreEqual(result.values.GetValue(1).ToString(), new double[2] { 6, 8 }.ToString());


        }


    }
}
