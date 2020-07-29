using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Objects;
using System;
using System.Drawing;



namespace NumericsTests
{
    [TestClass]
    public class ComplexNumberTests
    {

        [TestMethod]
        public void TestImaginaryPower()
        {
            var i = new ComplexNumber(3, 2);
            i = i.Pow(2);
            Assert.IsTrue(Math.Round(i.realPart) == 5);
            Assert.IsTrue(Math.Round(i.imaginaryPart) == 12);
        }

        [TestMethod]
        public void TestEulerIdentity()
        {
            var i = new ComplexNumber(0, Math.PI);
            i = i.Exponential();
            Assert.IsTrue(Math.Round(i.realPart) == -1);

        }

        [TestMethod]
        public void TestFromPolar()
        {
            var i = ComplexNumber.FromPolarCoordinates(2, Math.PI);

            Assert.IsTrue(Math.Round(i.realPart) == -2);
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
            Assert.IsTrue(result.realPart == 8);
            Assert.IsTrue(result.imaginaryPart == 5);
        }

        [TestMethod]
        public void TestComplexSubtraction()
        {
            var a = new ComplexNumber(3, 2);
            var b = new ComplexNumber(5, 3);
            var result = a - b;
            Assert.IsTrue(result.realPart == -2);
            Assert.IsTrue(result.imaginaryPart == -1);

        }

        [TestMethod]
        public void TestComplexMultiplication()
        {
            var a = new ComplexNumber(3, 2);
            var b = new ComplexNumber(5, 3);
            var result = a * b;
            Assert.IsTrue(result.realPart == 9);
            Assert.IsTrue(result.imaginaryPart == 19);

        }
        [TestMethod]
        public void TestComplexDivision()
        {
            var a = new ComplexNumber(3, 2);
            var b = new ComplexNumber(5, 3);
            var result = a / b;
            Assert.IsTrue(Math.Round(result.realPart, 2) == 0.62);
            Assert.IsTrue(Math.Round(result.imaginaryPart, 2) == 0.03);

        }


        [TestMethod]
        public void TestMandelbrot()
        {
            var maxValueExtent = 2.0;
            var  bitmap = new Bitmap(600, 600);
            var scale = 2 * maxValueExtent / Math.Min(bitmap.Width, bitmap.Height);
            for (int i = 0; i < bitmap.Height; i++)
            {
                var y = (bitmap.Height / 2 - i) * scale;
                for (var j = 0; j < bitmap.Width; j++)
                {
                    var x = (j - bitmap.Width / 2) * scale;

                   var maxIterations = 1000;
                   var maxNorm = maxValueExtent * maxValueExtent;

                    int iteration = 0;
                    ComplexNumber z = new ComplexNumber(0, 0);
                    while (z.GetMagnitude() < maxNorm && iteration < maxIterations)
                    {
                        z = z * z + new ComplexNumber(x, y);
                        iteration++;
                    }

                    var color = iteration < maxIterations ? (double)iteration / maxIterations : 0;
                    bitmap.SetPixel(j, i, Color.FromArgb(0, 0,
                (int)(256 * Math.Pow(color, 0.2))));
                }
            }

            bitmap.Save("mandelbrot.bmp");

        }

    }
}
