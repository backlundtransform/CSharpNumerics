using Xunit.Sdk;

using CSharpNumerics.Statistics.Data;
using CSharpNumerics.Statistics;


namespace NumericsTests
{
    [TestClass]
   public class StatisticsTest
    {
        [TestMethod]
        public void TestMedian()
        {
           var timeserie = new List<TimeSerie>() { new TimeSerie() { Value=5 }, new TimeSerie() { Value = 1} , new TimeSerie() { Value = 3 }, new TimeSerie() { Value = 2 }, new TimeSerie() { Value = 4 } };
           var median= timeserie.Median(p => p.Value);
            Assert.IsTrue(median == 3);
            timeserie.AddRange(new List<TimeSerie>() {  new TimeSerie() { Value =8 }, new TimeSerie() { Value = 9 }, new TimeSerie() { Value = 6 } });
            median = timeserie.Median(p => p.Value);
            Assert.IsTrue(median == 4.5);
        }


        [TestMethod]
        public void TestStandardDeviation()
        {
            var timeserie = new List<TimeSerie>() { new TimeSerie() { Value = 9 }, new TimeSerie() { Value = 2 }, new TimeSerie() { Value = 5 }, new TimeSerie() { Value =4 }, new TimeSerie() { Value = 12 },
                           new TimeSerie() { Value = 7 }, new TimeSerie() { Value = 8 }, new TimeSerie() { Value = 11 }, new TimeSerie() { Value =9 }, new TimeSerie() { Value = 3 }
                          ,new TimeSerie() { Value = 7 }, new TimeSerie() { Value = 4 }, new TimeSerie() { Value = 12 }, new TimeSerie() { Value =5 }, new TimeSerie() { Value = 4 }
                          ,new TimeSerie() { Value = 10 }, new TimeSerie() { Value = 9 }, new TimeSerie() { Value = 6 }, new TimeSerie() { Value =9 }, new TimeSerie() { Value = 4 } };

            var standardDeviation = timeserie.StandardDeviation(p => p.Value);
            Assert.IsTrue(Math.Round(standardDeviation,3) == 2.983);
    
        }


        [TestMethod]
        public void TestVariance()
        {
            var timeserie = new List<TimeSerie>() { new TimeSerie() { Value = 206 }, new TimeSerie() { Value = 76 }, new TimeSerie() { Value =-224}, new TimeSerie() { Value = 36 },
                           new TimeSerie() { Value = -94}};

            var variance = timeserie.Variance(p => p.Value);
            Assert.IsTrue(variance == 21704);

        }


        [TestMethod]
        public void TestCovariance()
        {
            var serie = new List<Serie>() { new Serie() { Index =2, Value=10 }, new Serie() { Index  = 3, Value = 14 }, new Serie() {Index  =2.7, Value = 12 }, new Serie() { Index = 3.2, Value = 15}, new Serie() { Index = 4.1, Value = 20 } };

            var covariance = serie.Covariance(p => (p.Index,p.Value));
            Assert.IsTrue(Math.Round(covariance,2) == 2.85);

        }



        [TestMethod]
        public void TestCumulativeSum()
        {
            var timeserie = new List<TimeSerie>() { new TimeSerie() { Value = 1.0 }, new TimeSerie() { Value = 2.0 }, new TimeSerie() { Value = 3.0 }, new TimeSerie() { Value =4.0 }, new TimeSerie() { Value = 5.0 }};

            var cumSum= timeserie.CumulativeSum(p => p.Value);
            Assert.IsTrue(Math.Round(cumSum.Sum(),1) == 35.0);

        }

    

        [TestMethod]
        public void R2_Uncorrelated()
        {
            var data = new[] { (1.0, 5.0), (2.0, 1.0), (3.0, 4.0), (4.0, 6.0) };
            double r2 = data.CoefficientOfDetermination(p => (p.Item1, p.Item2));

            Assert.IsTrue(r2 < 0.2); 
        }

        [TestMethod]
        public void Covariance_PerfectPositive()
        {
            var data = new[] { (1.0, 2.0), (2.0, 4.0), (3.0, 6.0) };

            double cov = data.Covariance(p => (p.Item1, p.Item2));

            Assert.AreEqual(2.0, cov, 1e-10);
        }

        [TestMethod]
        public void Covariance_ConstantX()
        {
            var data = new[] { (2.0, 1.0), (2.0, 3.0), (2.0, 5.0) };

            double cov = data.Covariance(p => (p.Item1, p.Item2));

            Assert.AreEqual(0.0, cov, 1e-10);
        }

        [TestMethod]
        public void StandardDeviation_Simple()
        {
            var values = new double[] { 2, 4, 4, 4, 5, 5, 7, 9 };
            double sd = values.StandardDeviation();

            Assert.AreEqual(2.138089935, sd, 1e-6);
        }


        [TestMethod]
        public void R2_PerfectNegative()
        {
            var data = new[] { (1.0, 10.0), (2.0, 8.0), (3.0, 6.0) };
            double r2 = data.CoefficientOfDetermination(p => (p.Item1, p.Item2));

            Assert.AreEqual(1.0, r2, 1e-10);
        }


        [TestMethod]
        public void R2_PerfectPositive()
        {
            var data = new[] { (1.0, 2.0), (2.0, 4.0), (3.0, 6.0) };
            double r2 = data.CoefficientOfDetermination(p => (p.Item1, p.Item2));

            Assert.AreEqual(1.0, r2, 1e-10);
        }


        [TestMethod]
        public void TestConfidence()
        {

            var timeserie = new List<TimeSerie>() { new TimeSerie() { Value = 9 }, new TimeSerie() { Value = 2 }, new TimeSerie() { Value = 5 }, new TimeSerie() { Value =4 }, new TimeSerie() { Value = 12 },
                           new TimeSerie() { Value = 7 }, new TimeSerie() { Value = 8 }, new TimeSerie() { Value = 11 }, new TimeSerie() { Value =9 }, new TimeSerie() { Value = 3 }
                          ,new TimeSerie() { Value = 7 }, new TimeSerie() { Value = 4 }, new TimeSerie() { Value = 12 }, new TimeSerie() { Value =5 }, new TimeSerie() { Value = 4 }
                          ,new TimeSerie() { Value = 10 }, new TimeSerie() { Value = 9 }, new TimeSerie() { Value = 6 }, new TimeSerie() { Value =9 }, new TimeSerie() { Value = 4 } };

            var (lower,upper) = timeserie.ConfidenceIntervals(p => p.Value, 0.95);

            Assert.IsTrue(Math.Round(lower, 1) ==5.7);
            Assert.IsTrue(Math.Round(upper, 1) == 8.3);

        }

      
   }
}
