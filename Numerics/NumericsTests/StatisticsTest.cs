using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Models;
using System;
using System.Collections.Generic;
using System.Linq;

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
    }
}
