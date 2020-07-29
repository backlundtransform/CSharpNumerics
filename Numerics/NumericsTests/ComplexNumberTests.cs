using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Objects;
using System;
using System.Drawing;
using System.Drawing.Imaging;


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

            var  bitmap = new Bitmap(300, 300, PixelFormat.Format24bppRgb);

            var scale = 0.1;
            var max = 255;
            var norm = 4;


            for (var i = 0; i< bitmap.Height; i++)
            {
              var y = (bitmap.Height / 2 - i) * scale;
                for (var j = 0; j < bitmap.Width; j++)
                {
                    var k = 0;
                    var x = (bitmap.Height / 2 - i)*scale;
                    var z = new ComplexNumber(0, 0);
                    while (z.GetMagnitude() < norm && k < max)
                    {

                        z = z.Pow(2) + new ComplexNumber(x, y);

                       k++;

                    }
                  
                    bitmap.SetPixel(j, i, k < max ? Color.FromArgb(0, 0,k) : Color.FromArgb(0, 0, 0));

                }
            }
            bitmap.Save("mandelbrot.bmp");

        }
    }
}
