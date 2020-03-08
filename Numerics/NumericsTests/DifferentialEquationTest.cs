using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;

namespace NumericsTests
{
    [TestClass]
    public class DifferentialEquationTest
    {
        [TestMethod]
        public void TestRungeKutta()
        {
            Func<(double y, double t), double> func = ((double t, double y) v) => Math.Tan(v.y) + 1;

            var result = func.RungeKutta(1, 1.1, 0.025, 1);
            Assert.IsTrue(Math.Round(result, 2) == 1.34);
        }

        [TestMethod]
        public void TestRungeKuttaMatrix()
        {
            Func<(double y, double t), double> func = ((double t, double y) v) => Math.Tan(v.y) + 1;

            var result = func.RungeKutta(1, 1.1, 0.025, 1, new Matrix(new double[,] { { 0, 0 }, { 2.0/3.0, 0 }}),new double[] { 1.0/4.0, 3.0/4.0 }, new double[] { 0.0, 2.0/3.0 });
            Assert.IsTrue(Math.Round(result, 2) == 1.34);
        }

        [TestMethod]
        public void TestTrapets()
        {
            Func<(double y, double t), double> func = ((double t, double y) v) => Math.Tan(v.y) + 1;

            var result = func.TrapezoidalRule(1, 1.1, 0.00025, 1);
            Assert.IsTrue(Math.Round(result,2) == 1.34);
        }

        [TestMethod]
        public void TestLinearSystem()
        {
            var matrix = new Matrix(new double[,] { { 1, -2, 3 }, { -1, 1, -2 }, { 2, -1, -1 } });

            var vector = new Vector(7, -5, 4);

            var result = matrix.LinearSystemSolver(vector);
            Assert.IsTrue(result.x == 2);
            Assert.IsTrue(result.y == -1);
            Assert.IsTrue(result.z== 1);
        }

        [TestMethod]
        public void TestLinearSystemRange()
        {
            var matrix = new Matrix(new double[,] { { 1, -2, 3 }, { -1, 1, -2 }, { 2, -1, -1 } });

            var vector = new List<double>() { { 7 },{ -5 }, { 4 } };

            var result = matrix.LinearSystemSolver(vector);
            Assert.IsTrue(result[0] == 2);
            Assert.IsTrue(result[1] == -1);
            Assert.IsTrue(result[2] == 1);
        }


        [TestMethod]
        public void TestEigenvalueLarge()
        {
            var matrix = new Matrix(new double[,] { {3, -1 }, { 4, -2 }});
            var result = matrix.LargestEigenValue();
            Assert.IsTrue(result ==2);
           
        }

        [TestMethod]
        public void TestEigenvalueSmall()
        {
            var matrix = new Matrix(new double[,] { { 3, -1 }, { 4, -2 } });
            var result = matrix.SmallestEigenValue();
            Assert.IsTrue(result == -1);

        }


        [TestMethod]
        public void TestEigenvalues()
        {
            var matrix = new Matrix(new double[,] { { 3, -1 }, { 4, -2 } });

            var result = matrix.EigenValues();
            Assert.IsTrue(result.First() == -1);
            Assert.IsTrue(result.Last() == 2);

        }

    }
}
