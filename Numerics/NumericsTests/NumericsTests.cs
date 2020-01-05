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
            Assert.IsTrue(Math.Round(result) == Math.Round(g * t));
        }
        [TestMethod]
        public void TestIntegrateExponentiation()
        {
            Func<double, double> func = (double variable) => g * variable;
            var lowerlimit = 0;
            var upperlimit = 5;
            var result = func.Integrate(lowerlimit, upperlimit);
            Assert.IsTrue(Math.Round(result) == Math.Round(g * Math.Pow(upperlimit, 2) / 2) - Math.Round(g * Math.Pow(lowerlimit, 2) / 2));
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
        public void TestIntegrateExponential()
        {
            Func<double, double> func = (double variable) => Math.Pow(Math.E, variable);
            var lowerlimit = 0;
            var upperlimit = 5;
            var result = func.Integrate(lowerlimit, upperlimit);
            Assert.IsTrue(Math.Round(result,2) == Math.Round(Math.Pow(Math.E, upperlimit) - Math.Pow(Math.E, lowerlimit),2));
        }

        [TestMethod]
        public void TestDerivateSinus()
        {
            Func<double, double> func = (double variable) => Math.Sin(variable);
            var t = Math.PI;
            var result = func.Derivate(t);
            Assert.IsTrue(Math.Round(result) == Math.Round(Math.Cos(t)));
        }

        [TestMethod]
        public void TestIntegrateCosinus()
        {
            Func<double, double> func = (double variable) => Math.Round(Math.Cos(variable));
            var lowerlimit = 0;
            var upperlimit = Math.PI/5;
            var result = func.Integrate(lowerlimit, upperlimit);
            Assert.IsTrue(Math.Round(result) == Math.Round(Math.Sin(upperlimit)- Math.Sin(lowerlimit)));
        }


    }
}
