using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Objects;
using System;

namespace NumericsTests
{
    [TestClass]
    public class ComplexNumberTests
    {

        [TestMethod]
        public void TestImaginaryPower()
        {
            var i = new ComplexNumber(3, 2);
            i.Pow(2);
            Assert.IsTrue(Math.Round(i.RealPart) == 5);
            Assert.IsTrue(Math.Round(i.ImaginaryPart) == 12);
        }

        [TestMethod]
        public void TestEulerIdentity()
        {
            var i = new ComplexNumber(0, Math.PI);
            i.Exponential();
            Assert.IsTrue(Math.Round(i.RealPart) == -1);

        }

        [TestMethod]
        public void TestFromPolar()
        {
            var i =  ComplexNumber.FromPolarCoordinates(2, Math.PI);
     
            Assert.IsTrue(Math.Round(i.RealPart) == -2);
        }

        [TestMethod]
        public void TestComplexToString()
        {
            var i = new ComplexNumber(3, 2);

            Assert.IsTrue(i.ToString() == "3+2*i");

        }

        [TestMethod]
        public void TestComplexAddition()
        {
            var a = new ComplexNumber(3, 2);
            var b = new ComplexNumber(5, 3);
            var result = a + b;
            Assert.IsTrue(result.RealPart == 8);
            Assert.IsTrue(result.ImaginaryPart == 5);
        }

        [TestMethod]
        public void TestComplexSubtraction()
        {
            var a = new ComplexNumber(3, 2);
            var b = new ComplexNumber(5, 3);
            var result = a - b;
            Assert.IsTrue(result.RealPart == -2);
            Assert.IsTrue(result.ImaginaryPart == -1);

        }

        [TestMethod]
        public void TestComplexMultiplication()
        {
            var a = new ComplexNumber(3, 2);
            var b = new ComplexNumber(5, 3);
            var result = a * b;
            Assert.IsTrue(result.RealPart == 9);
            Assert.IsTrue(result.ImaginaryPart == 19);

        }
        [TestMethod]
        public void TestComplexDivision()
        {
            var a = new ComplexNumber(3, 2);
            var b = new ComplexNumber(5, 3);
            var result = a / b;
            Assert.IsTrue(Math.Round(result.RealPart, 2) == 0.62);
            Assert.IsTrue(Math.Round(result.ImaginaryPart, 2) == 0.03);

        }
    }
}
