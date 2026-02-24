using Xunit.Sdk; 
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics;


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

            Assert.IsTrue(w.IsAnalytical((3, 3)));


            ComplexNumber fz(ComplexNumber z) => new ComplexNumber(Math.Pow(Math.E, z.realPart) * Math.Cos(z.imaginaryPart), Math.Pow(Math.E, z.realPart) * Math.Sin(z.imaginaryPart));

            w = new ComplexFunction(fz);

            Assert.IsTrue(w.IsAnalytical((5,6)));
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

    }
}
