using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Objects;
using System;


namespace NumericsTests
{
   [TestClass]
   public class ComplexFunctionTest
    {
        [TestMethod]
        public void TestIsAnalytical()
        {

            double fx((double x, double y) p) => Math.Pow(Math.E ,p.x) *Math.Cos(p.y);
            double fy((double x, double y) p) => Math.Pow(Math.E, p.x) * Math.Sin(p.y);

            var w = new ComplexFunction(fx, fy);

            var rnd = new Random();

            Assert.IsTrue(w.IsAnalytical((rnd.Next(10), rnd.Next(10))));


            ComplexNumber fz(ComplexNumber z) => new ComplexNumber(Math.Pow(Math.E, z.realPart) * Math.Cos(z.imaginaryPart), Math.Pow(Math.E, z.realPart) * Math.Sin(z.imaginaryPart));

           w = new ComplexFunction(fz);


            Assert.IsTrue(w.IsAnalytical((rnd.Next(10), rnd.Next(10))));
        }

        [TestMethod]
        public void TestDerivate()
        {

            double fx((double x, double y) p) => Math.Pow(Math.E, p.x) * Math.Cos(p.y);

            double fy((double x, double y) p) => Math.Pow(Math.E, p.x) * Math.Sin(p.y);

            var fz = new ComplexFunction(fx, fy);

            var complexnumber = new ComplexNumber(fz.u((2,3)), fz.v((2, 3)));

            var complexderivate = fz.Derivate(new ComplexNumber(2, 3));

            Assert.IsTrue(Math.Round(complexnumber.realPart,2) == Math.Round(complexderivate.realPart, 2));
            Assert.IsTrue(Math.Round(complexnumber.imaginaryPart, 2) == Math.Round(complexderivate.imaginaryPart, 2));

        }


        [TestMethod]
        public void TestJacobian()
        {

            double fx((double x, double y) p) => Math.Pow(Math.E, p.x) * Math.Cos(p.y);
            double fy((double x, double y) p) => Math.Pow(Math.E, p.x) * Math.Sin(p.y);

            var fz = new ComplexFunction(fx, fy);

            var test = Math.Pow(fz.Derivate(new ComplexNumber(1, 3)).GetMagnitude(), 2);
            var jacobian = fz.Jacobian((1, 3)).Determinant();
            Assert.IsTrue(Math.Round(jacobian, 2) == Math.Round(test, 2));


        }

    }
}
