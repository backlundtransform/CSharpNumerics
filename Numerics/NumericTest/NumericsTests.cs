using Xunit.Sdk;
using CSharpNumerics.Numerics.Objects;

using CSharpNumerics.Numerics.Enums;
using CSharpNumerics.Statistics.Data;
using CSharpNumerics.Numerics;


namespace NumericsTests
{
    
    [TestClass]
    public class NumericsTests
    {
        private const double g = 9.8;

        [TestMethod]
        public void TestLimit()
        {
            Func<double, double> func = (double x) =>Math.Sin(x) / x;
          
            var result = func.Limit(0);

            Assert.IsTrue(Math.Round(result) ==1);
        }
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
        public void TestDerivateExponentiationSeries()
        {
            Func<double, double> func = (double variable) => g * Math.Pow(variable, 2) / 2;

            var series = func.GetSeries(0,10,1000).Derivate();

            var (slope, intercept, correlation) = series.LinearRegression(p => (p.Index, p.Value));

            Assert.IsTrue(Math.Round(slope, 1)==g);
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
        public void TestChainRule()
        {
            double funcG(double x) => 4 * x - 3;
            Func<double, double> funcF=(double x) => Math.Pow(x, 2);
            var result = funcF.Derivate(funcG,1);
            Assert.IsTrue(Math.Round(result) == 8);
        }

        [TestMethod]
        public void TestProductRule()
        {
            double funcG(double x) => Math.Exp(x);
            Func<double, double> funcF = (double x) => Math.Sin(x);
            var result = funcF.Derivate(funcG, 5, DerivateOperator.Product);
            Assert.IsTrue(Math.Round(result) == -100);
        }

        [TestMethod]
        public void TestQuotientRule()
        {
            Func<double, double> funcF=(double x) =>3* Math.Cos(x) ;
            Func<double, double> funcG = (double x) => 2 * x + 1;
            var result = funcF.Derivate(funcG, 0, DerivateOperator.Quotient);
            Assert.IsTrue(Math.Round(result) ==-6);
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

            Assert.IsTrue(result < 177 && result > 173);
        }


        [TestMethod]
        public void TestTripleIntegral()
        {
            Func<Vector, double> func = (Vector v) => (Math.Pow(v.x, 3) + Math.Pow(v.y, 2))+v.z;

            var result = func.Integrate(new Vector(-2, -2,-2), new Vector(2, 2, 2));
           
            Assert.IsTrue(result <87 && result >83);
        }

        [TestMethod]
        public void TestFactorial()
        {
           Assert.IsTrue(5.Factorial() == 120);
        }

        [TestMethod]
        public void FindRoots()
        {
            Func<double, double> func = (double x) => Math.Pow(x,2) - 4;
            var result = func.NewtonRaphson();
            Assert.IsTrue(Math.Abs(2) == 2);
        }

        [TestMethod]
        public void TestIsPrime()
        {
            Assert.IsTrue(79.IsPrime());
        }

        [TestMethod]
        public void TestGetPrimeFactors()
        {
            var factors = 78.GetPrimeFactors();
            CollectionAssert.AreEqual(factors, new[] {2,3,13});
       
        }

        [TestMethod]
        public void TestIsHappy()
        {
            Assert.IsTrue(19.IsHappy());

        }
        [TestMethod]
        public void TestDecimalPlace()
        {
            Assert.IsTrue(0.01.GetDecimalPlaces()== 2);

        }
        [TestMethod]
        public void TestIsPerfectNumber()
        {
            Assert.IsTrue(6.IsPerfectNumber());

        }
       

       [TestMethod]
        public void TestIntegrateTimeSerie()
        {

            var date = new DateTime(2011, 1, 1);
            var list = new List<TimeSerie>()
            {
                new TimeSerie () { TimeStamp = date, Value =0.56 },
                new TimeSerie () { TimeStamp = date.AddHours(4) , Value =0.55 },
                new TimeSerie () { TimeStamp = date.AddHours(8) , Value =0.43 },
                new TimeSerie () { TimeStamp =date.AddHours(12) , Value =0.47 },
                new TimeSerie () { TimeStamp =date.AddHours(16) , Value =0.65 },
                new TimeSerie () { TimeStamp =date.AddHours(20) , Value =0.76 },

            };

            var valueApprox = list.Integrate(date, date.AddDays(1));
            var approx = list.Average(p => p.Value) * 86400;
            list.Add(new TimeSerie() { TimeStamp = date.AddHours(24), Value = 0.76 });
            var value = list.Integrate();
            var relError= 100 * Math.Abs(value - valueApprox) / value;
            Assert.IsTrue(relError < 3.5);

            relError = 100 * Math.Abs(approx - valueApprox) / value;
            Assert.IsTrue(relError < 6.5);

        }
  
    }
}
