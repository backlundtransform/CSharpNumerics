
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;

namespace NumericsTests
{
    [TestClass]
    public class TransFormTest
    {

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
