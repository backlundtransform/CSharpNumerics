using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Objects;
using System;

namespace NumericsTests
{
    [TestClass]
   public class VectorFieldTest
    {
        [TestMethod]
        public void TestGradient()
        {
            Func<Vector, double> func = (Vector p) => Math.Pow(p.x, 2) * Math.Pow(p.y, 3);
            var v=func.Gradient((1, -2, 0));
            Assert.IsTrue(Math.Round(v.x) == -16);
            Assert.IsTrue(Math.Round(v.y) == 12);
            Assert.IsTrue(Math.Round(v.z) == 0);
        }

        [TestMethod]
        public void TestDivergence()
        {
            double fx(Vector p) => Math.Sin(p.x * p.y);
            double fy(Vector p) => Math.Cos(p.x * p.y);
            double fz(Vector p) => Math.Pow(Math.E, p.z);
            var field = new VectorField(fx, fy, fz);
            var div = field.Divergence((1, 2, 2));
            Assert.IsTrue(Math.Round(div,2) == 5.65);
          
        }

        [TestMethod]
        public void TestCurl()
        {
            double fx(Vector p) => 4*p.z;
            double fy(Vector p) => p.y *Math.Pow(p.x,3);
            double fz(Vector p) => p.z * Math.Pow(p.y,2);
            var field = new VectorField(fx, fy, fz);
            var v = field.Curl((1, 4, 2));
            Assert.IsTrue(Math.Round(v.x) == 16);
            Assert.IsTrue(Math.Round(v.y) == 4);
            Assert.IsTrue(Math.Round(v.z) == 12);
        }

    }
}
