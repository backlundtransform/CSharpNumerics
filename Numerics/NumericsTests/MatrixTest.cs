using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Objects;


namespace NumericsTests
{
    [TestClass]
    public class MatrixTest
    {
        [TestMethod]
        public void TestMatrixTranspose()
        {
            var matrix= new Matrix(new double[,] { { 1, 3, 7 },  { 5, 2, 9} });
            var transposematrix = matrix.Transpose();

            Assert.IsTrue(transposematrix.Values[0,0] == 1);
            Assert.IsTrue(transposematrix.Values[0,1] == 5);

            Assert.IsTrue(transposematrix.Values[1,0] == 3);
            Assert.IsTrue(transposematrix.Values[1,1] == 2);

            Assert.IsTrue(transposematrix.Values[2,0] == 7);
            Assert.IsTrue(transposematrix.Values[2,1] == 9);

        }

        [TestMethod]
        public void TestMatrixIdentity()
        {
            var matrix = new Matrix(new double[,] { { 1, 3, 7 }, { 5, 2, 9 } });
       
            Assert.IsTrue(matrix.Identity[0, 0] == 1);
            Assert.IsTrue(matrix.Identity[0, 1] == 0);
            Assert.IsTrue(matrix.Identity[0, 2] == 0);

            Assert.IsTrue(matrix.Identity[1, 0] == 0);
            Assert.IsTrue(matrix.Identity[1, 1] == 1);
            Assert.IsTrue(matrix.Identity[1, 2] == 0);

        }

    }
}
