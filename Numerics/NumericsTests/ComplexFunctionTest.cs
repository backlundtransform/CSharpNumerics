using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;

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

            ComplexNumber fz(ComplexNumber z) => new ComplexNumber(Math.Pow(Math.E, z.realPart) * Math.Cos(z.imaginaryPart), Math.Pow(Math.E, z.realPart) * Math.Sin(z.imaginaryPart));

            var w = new ComplexFunction(fz);

            var complexnumber = new ComplexNumber(w.u((2,3)), w.v((2, 3)));

            var complexderivate = w.Derivate(new ComplexNumber(2, 3));

            Assert.IsTrue(Math.Round(complexnumber.realPart,2) == Math.Round(complexderivate.realPart, 2));
            Assert.IsTrue(Math.Round(complexnumber.imaginaryPart, 2) == Math.Round(complexderivate.imaginaryPart, 2));

        }


        [TestMethod]
        public void TestJacobian()
        {

            ComplexNumber fz(ComplexNumber z) => new ComplexNumber(Math.Pow(Math.E, z.realPart) * Math.Cos(z.imaginaryPart), Math.Pow(Math.E, z.realPart) * Math.Sin(z.imaginaryPart));

            var w = new ComplexFunction(fz);

            var test = Math.Pow(w.Derivate(new ComplexNumber(1, 3)).GetMagnitude(), 2);
            var jacobian = w.Jacobian((1, 3)).Determinant();
            Assert.IsTrue(Math.Round(jacobian, 2) == Math.Round(test, 2));

        }


        [TestMethod]
        public void TestFourierTransform()
        {

            var input = new List<ComplexNumber>() { new ComplexNumber(1, 0), new ComplexNumber(1, 0), new ComplexNumber(1, 0), new ComplexNumber(1, 0), new ComplexNumber(0, 0), new ComplexNumber(0, 0), new ComplexNumber(0, 0), new ComplexNumber(0, 0) };

            var transform = input.Fouriertransform().ToList();

            Assert.IsTrue(transform[0].realPart == 4);
            Assert.IsTrue(transform[0].imaginaryPart == 0);

            Assert.IsTrue(transform[1].realPart == 1);
            Assert.IsTrue(Math.Round(transform[1].imaginaryPart, 2) == -2.41);


            Assert.IsTrue(transform[2].realPart == 0);
            Assert.IsTrue(transform[2].imaginaryPart == 0);


            Assert.IsTrue(transform[3].realPart == 1);
            Assert.IsTrue(Math.Round(transform[3].imaginaryPart, 2) == -0.41);


            Assert.IsTrue(transform[4].realPart == 0);
            Assert.IsTrue(transform[4].imaginaryPart == 0);


            Assert.IsTrue(Math.Round(transform[5].realPart) == 1);
            Assert.IsTrue(Math.Round(transform[5].imaginaryPart, 2) == 0.41);

            Assert.IsTrue(transform[6].realPart == 0);
            Assert.IsTrue(transform[6].imaginaryPart == 0);

            Assert.IsTrue(Math.Round(transform[7].realPart) == 1);
            Assert.IsTrue(Math.Round(transform[7].imaginaryPart, 2) == 2.41);

            input = transform.InverseFouriertransform();


            Assert.IsTrue(Math.Round(input[0].realPart) == 1);
            Assert.IsTrue(Math.Round(input[0].imaginaryPart) == 0);

            Assert.IsTrue(Math.Round(input[1].realPart) == 1);
            Assert.IsTrue(Math.Round(input[1].imaginaryPart) == 0);


            Assert.IsTrue(Math.Round(input[2].realPart) == 1);
            Assert.IsTrue(Math.Round(input[2].imaginaryPart) == 0);


            Assert.IsTrue(Math.Round(input[3].realPart) == 1);
            Assert.IsTrue(Math.Round(input[3].imaginaryPart) == 0);


            Assert.IsTrue(Math.Round(input[4].realPart) == 0);
            Assert.IsTrue(Math.Round(input[4].imaginaryPart) == 0);


            Assert.IsTrue(Math.Round(input[5].realPart) == 0);
            Assert.IsTrue(Math.Round(input[5].imaginaryPart) == 0);

            Assert.IsTrue(Math.Round(input[6].realPart) == 0);
            Assert.IsTrue(Math.Round(input[6].imaginaryPart) == 0);

            Assert.IsTrue(Math.Round(input[7].realPart) == 0);
            Assert.IsTrue(Math.Round(input[7].imaginaryPart) == 0);


        }

    }
}
