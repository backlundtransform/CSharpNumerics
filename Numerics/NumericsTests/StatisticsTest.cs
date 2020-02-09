using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Models;
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

    }
}
