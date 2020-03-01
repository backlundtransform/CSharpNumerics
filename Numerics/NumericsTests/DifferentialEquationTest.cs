using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Objects;
using System;

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
        public void TestRungeKutta2()
        {
            Func<(double y, double t), double> func = ((double t, double y) v) => Math.Tan(v.y) + 1;

            var result = func.RungeKutta(1, 1.1, 0.025, 1, new Matrix(new double[,] { { 0, 0 }, { 2.0/3.0, 0 }}),new double[] { 1.0/4.0, 3.0/4.0 }, new double[] { 0.0, 2.0/3.0 });
            Assert.IsTrue(Math.Round(result, 2) == 1.34);
        }

    }
}
