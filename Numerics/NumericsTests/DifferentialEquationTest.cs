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
            Func<(double y, double t), double> func = ((double t, double y) v) => Math.Tan(v.y) +1;
         
            var result = func.RungeKutta(1,1.1,0.025,1);
            Assert.IsTrue(Math.Round(result, 3) ==1.335);
        }
    }
}
