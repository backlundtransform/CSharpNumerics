using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Objects;

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
        public void TestSecondDerivateExponentiation()
        {
            Func<double, double> func = (double variable) => g * Math.Pow(variable, 2) / 2;
            var t = 5;
            var result = func.Derivate(t,2);
            Assert.IsTrue(Math.Round(result) == Math.Round(g));
        }

        [TestMethod]
        public void TestThirdDerivateExponentiation()
        {
            Func<double, double> func = (double variable) => g * Math.Pow(variable, 3) / 6;
            var t = 5;
            var result = func.Derivate(t,3);
            Assert.IsTrue(Math.Round(result) == Math.Round(g));
        }

        [TestMethod]
        public void TestFourthDerivateExponentiation()
        {
            Func<double, double> func = (double variable) => g * Math.Pow(variable, 4) / 24;
            var t = 5;
            var result = func.Derivate(t, 4);
            Assert.IsTrue(Math.Round(result) == Math.Round(g));
        }

        [TestMethod]
        public void TestFifthDerivateExponentiation()
        {
            Func<double, double> func = (double variable) => g * Math.Pow(variable, 5) / 120;
            var t = 5;
            var result = func.Derivate(t, 5);
            Assert.IsTrue(Math.Round(result) == Math.Round(g));
        }


        [TestMethod]
        public void TestDerivatePartial()
        {
            Func<double[], double> func = (double[] variables) =>Math.Pow(variables[0], 2) + variables[1] * variables[0] + Math.Pow(variables[1], 2);
            var result = func.Derivate(new double[] {1,1}, 0);
            Assert.IsTrue(Math.Round(result) == 3);
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
        public void TestDerivateLinear()
        {
            Func<double, double> func = (double variable) => g * variable;
            var t = 5;
            var result = func.Derivate(t);
            Assert.IsTrue(Math.Round(result, 2) == Math.Round(g, 2));
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


        [TestMethod]
        public void TestDoubleIntegrate()
        {
            Func<(double, double), double> func = ((double x, double y) v) => (Math.Pow(v.x, 3) + Math.Pow(v.y, 2));
    
            var result = func.Integrate((1,4),(1,4));
            Assert.IsTrue(Math.Truncate(result) == 254);

            Func<(double, double), double> func2 = ((double x, double y) v) => 2*v.x*v.y + Math.Pow(v.y, 2);

            result = func2.Integrate((1, 4), (1, 4));
            Assert.IsTrue(Math.Truncate(result) == 175);
        }


        [TestMethod]
        public void TestTripleIntegral()
        {
            Func<Vector, double> func = (Vector v) => (Math.Pow(v.x, 3) + Math.Pow(v.y, 2))+v.z;

            var result = func.Integrate(new Vector(-2, -2,-2), new Vector(2, 2, 2));
           
            Assert.IsTrue(Math.Truncate(result) ==85);
        }

        [TestMethod]
        public void TestFactorial()
        {
           Assert.IsTrue(5.Factorial() == 120);
        }

    }
}
