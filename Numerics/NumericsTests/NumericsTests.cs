using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;


namespace NumericsTests
{
    
    [TestClass]
    public class NumericsTests
    {
        private const double g = 9.8;
        [TestMethod]
        public void TestDerivateExponentiation()
        {
            Func<double, double> func = (double variable) => g * Math.Pow(variable, 2) / 2;
            var t = 5;
            var result = func.Derivate(t);
            Assert.IsTrue(Math.Round(result,2) == Math.Round(g * t, 2));
        }
        [TestMethod]
        public void TestDerivateExponential()
        {
            Func<double, double> func = (double variable) => Math.Pow(Math.E, variable);
            var t = 5;
            var result = func.Derivate(t);
            Assert.IsTrue(Math.Round(result, 2) == Math.Round(Math.Pow(Math.E, t),2));
        }

        [TestMethod]
        public void TestDerivateSinus()
        {
            Func<double, double> func = (double variable) => Math.Sin(variable);
            var t = Math.PI;
            var result = func.Derivate(t);
            Assert.IsTrue(Math.Round(result) == Math.Round(Math.Cos(t)));
        }
    }
}
