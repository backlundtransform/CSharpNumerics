using CSharpNumerics.Numerics.Enums;
using CSharpNumerics.Numerics.Models;

namespace NumericTest
{
    [TestClass]
    public class InterpolationTest
    {
        [TestMethod]
        public void TestLinearInterpolationTimeSerie()
        {
            var timeserie = new List<TimeSerie>() { new TimeSerie() { TimeStamp = new DateTime(2020, 01, 01), Value = 1.0 }, new TimeSerie() { TimeStamp = new DateTime(2020, 01, 30), Value = 2.0 } };

            var value = timeserie.LinearInterpolationTimeSerie(new DateTime(2020, 01, 15));
            Assert.IsTrue(Math.Round(value, 1) == 1.5);

        }

        [TestMethod]
        public void TestlinearInterpolationSerie()
        {
            var serie = new List<Serie>() { new Serie() { Index = 0, Value = 1.0 }, new Serie() { Index = 2, Value = 2.0 }, new Serie() { Index = 3, Value = 2.5 } };

            var value = serie.LinearInterpolation(p => (p.Index, p.Value), 1);
            Assert.IsTrue(Math.Round(value, 1) == 1.5);

            value = serie.Interpolate(p => (p.Index, p.Value), 1, InterpolationType.Linear);
            Assert.IsTrue(Math.Round(value, 1) == 1.5);
        }

        [TestMethod]
        public void TestLogarithmicInterpolationSerie()
        {
            var serie = new List<Serie>()
            {
                 new() { Index = 1, Value = 10.0 },
                 new() { Index = 10, Value = 100.0 }
            };


            double index = Math.Sqrt(10);

            var value = serie.LogarithmicInterpolation(p => (p.Index, p.Value), index);

            double expected = Math.Sqrt(10 * 100); 

            Assert.IsTrue(Math.Abs(value - expected) < 1e-6,
                $"Expected {expected}, got {value}");

            value = serie.Interpolate(p => (p.Index, p.Value), index, InterpolationType.Logarithmic);

            Assert.IsTrue(Math.Abs(value - expected) < 1e-6,
             $"Expected {expected}, got {value}");

        }

        [TestMethod]
        public void TestLinLogInterpolation()
        {

            var serie = new List<Serie>()
    {
        new Serie { Index = 0, Value = 10 },
        new Serie { Index = 10, Value = 100 }
    };


            double index = 5;


            double value = serie.LinLogInterpolation(p => (p.Index, p.Value), index);

            double ly1 = Math.Log(10);
            double ly2 = Math.Log(100);
            double t = (index - 0.0) / (10.0 - 0.0);
            double expected = Math.Exp(ly1 + t * (ly2 - ly1));

            Assert.IsTrue(Math.Abs(value - expected) < 1e-9);

            value = serie.Interpolate(p => (p.Index, p.Value), index, InterpolationType.LinLog);
            Assert.IsTrue(Math.Abs(value - expected) < 1e-9);
        }
        [TestMethod]
        public void TestLogLinInterpolation()
        {

            var serie = new List<Serie>()
    {
        new Serie { Index = 1, Value = 10 },
        new Serie { Index = 10, Value = 20 }
    };


            double index = Math.Sqrt(10);


            double value = serie.LogLinInterpolation(p => (p.Index, p.Value), index);


            double t = (Math.Log(index) - Math.Log(1)) / (Math.Log(10) - Math.Log(1));
            double expected = 10 + t * (20 - 10);

            Assert.IsTrue(Math.Abs(value - expected) < 1e-9);


            value = serie.Interpolate(p => (p.Index, p.Value), index, InterpolationType.LogLin);

            Assert.IsTrue(Math.Abs(value - expected) < 1e-9);
        }
    }
}
