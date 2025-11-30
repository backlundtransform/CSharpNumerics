using Xunit.Sdk;
using Numerics.Methods;
using Numerics.Models;


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
        public void NormalDistribution()
        {
            var normalDistribution = Statistics.NormalDistribution(15, 100).GetSeries(0, 200, 100);
           
            Assert.IsTrue(Math.Round(normalDistribution.First(p=>p.Value==normalDistribution.Max(c=>c.Value)).Index) == 100);

         ;

        }


        [TestMethod]
        public void TestCumulativeSum()
        {
            var timeserie = new List<TimeSerie>() { new TimeSerie() { Value = 1.0 }, new TimeSerie() { Value = 2.0 }, new TimeSerie() { Value = 3.0 }, new TimeSerie() { Value =4.0 }, new TimeSerie() { Value = 5.0 }};

            var cumSum= timeserie.CumulativeSum(p => p.Value);
            Assert.IsTrue(Math.Round(cumSum.Sum(),1) == 35.0);

        }

        [TestMethod]
        public void TestLinearInterpolationTimeSerie()
        {
            var timeserie = new List<TimeSerie>() { new TimeSerie() { TimeStamp=new DateTime(2020,01,01),  Value = 1.0 }, new TimeSerie() { TimeStamp = new DateTime(2020, 01, 30), Value = 2.0 } };

            var value = timeserie.LinearInterpolationTimeSerie(new DateTime(2020, 01, 15));
            Assert.IsTrue(Math.Round(value, 1) == 1.5);

        }

        [TestMethod]
        public void TestlinearInterpolationSerie()
        {
            var serie = new List<Serie>() { new Serie() { Index= 0, Value = 1.0 }, new Serie() { Index = 2, Value = 2.0 }, new Serie() { Index = 3, Value = 2.5 } };

            var value = serie.LinearInterpolation(p=>(p.Index, p.Value),1);
            Assert.IsTrue(Math.Round(value, 1) == 1.5);

        }

        [TestMethod]
        public void R2_Uncorrelated()
        {
            var data = new[] { (1.0, 5.0), (2.0, 1.0), (3.0, 4.0), (4.0, 6.0) };
            double r2 = data.CoefficientOfDetermination(p => (p.Item1, p.Item2));

            Assert.IsTrue(r2 < 0.2); // slumpdata → låg förklaringsgrad
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
        public void TestLinearRegression()
        {
            var serie = new List<Serie>() {
              new Serie() { Index = 3.0, Value = 0.62},
              new Serie() { Index = 3.4, Value = 0.93 },
              new Serie() { Index = 3.8, Value = 1.08 },
              new Serie() { Index = 4.2, Value = 1.19},
              new Serie() { Index = 4.6, Value = 1.45 },
              new Serie() { Index = 5.0, Value = 1.54 },
              new Serie() { Index = 5.4, Value = 1.62 },
              new Serie() { Index = 5.8, Value = 1.92 },
              new Serie() { Index = 6.2, Value = 1.96 },
              new Serie() { Index = 6.6, Value = 2.10 },
              new Serie() { Index = 7.0, Value = 2.35},
              new Serie() { Index = 7.4, Value = 2.49 },
              new Serie() { Index = 7.8, Value = 2.58 } };


            var (slope, intercept, correlation) = serie.LinearRegression(p=>(p.Index, p.Value));
            Assert.IsTrue(Math.Round(slope, 3) == 0.395);
            Assert.IsTrue(Math.Round(intercept, 3) == -0.455);
            Assert.IsTrue(Math.Round(correlation, 3) == 0.995);
        }


        [TestMethod]
        public void LogisticRegression()
        {
            var series = new List<Serie>() {
              new Serie() { Index = 1, Value = 0.5},
              new Serie() { Index = 0, Value = 0.75 },
              new Serie() { Index = 0, Value = 1.0},
              new Serie() { Index = 0, Value = 1.25},
              new Serie() { Index = 0, Value = 1.50 },
              new Serie() { Index =1, Value = 1.75 },
              new Serie() { Index = 0, Value = 1.75 },
              new Serie() { Index = 0, Value = 2 },
              new Serie() { Index = 1, Value = 2.25 },
              new Serie() { Index = 0, Value = 2.5 },
              new Serie() { Index = 1, Value = 2.75 },
              new Serie() { Index = 0, Value = 3},
              new Serie() { Index = 1, Value = 3.25 },
              new Serie() { Index = 0, Value = 3.50 },
              new Serie() { Index = 1, Value = 4 },
              new Serie() { Index = 1, Value = 4.25 },
              new Serie() { Index =1, Value = 4.50 },
              new Serie() { Index = 1, Value = 4.75 },
              new Serie() { Index = 1, Value = 5 },
              new Serie() { Index = 1, Value = 5.5 }
            };
           
            var result = series.LogisticRegression(p => (p.Index, p.Value), 1.5046, -4.0777);

       
            Assert.IsTrue(Math.Round(result.First(p=>p.Index==2).Value, 2) == 0.26);
  
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

        [TestMethod]
        public void TestNearestNeighbors()
        {

            var timeserie = new List<(double x, double y, int classification)>() { (7, 7, 0), (7, 4, 0), (3, 4, 1), (1, 4, 1) };

            var classification = timeserie.KnearestNeighbors(p=> (p.x, p.y, p.classification),(3,7),3);

            Assert.IsTrue(classification == 1);

        }
   }
}
