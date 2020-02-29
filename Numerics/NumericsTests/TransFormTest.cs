using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;

namespace NumericsTests
{
    [TestClass]
    public class TransFormTest
    {
        [TestMethod]
        public void TestFastFourierTransform()
        {
            var input = new List<ComplexNumber>() { new ComplexNumber(1, 0), new ComplexNumber(1, 0), new ComplexNumber(1, 0), new ComplexNumber(1, 0), new ComplexNumber(0, 0), new ComplexNumber(0, 0), new ComplexNumber(0, 0), new ComplexNumber(0, 0) };

            var transform = input.FastFourierTransform().ToList();

            Assert.IsTrue(transform[0].realPart == 4);
            Assert.IsTrue(transform[0].imaginaryPart == 0);

            Assert.IsTrue(transform[1].realPart == 1);
            Assert.IsTrue(Math.Round(transform[1].imaginaryPart, 2) == -2.41);

            Assert.IsTrue(transform[2].realPart == 0);
            Assert.IsTrue(transform[2].imaginaryPart == 0);

            Assert.IsTrue(transform[3].realPart == 1);
            Assert.IsTrue(Math.Round(transform[3].imaginaryPart, 2) == -0.41);

            Assert.IsTrue(transform[4].realPart == 0);
            Assert.IsTrue(transform[4].imaginaryPart == 0);

            Assert.IsTrue(Math.Round(transform[5].realPart) == 1);
            Assert.IsTrue(Math.Round(transform[5].imaginaryPart, 2) == 0.41);

            Assert.IsTrue(transform[6].realPart == 0);
            Assert.IsTrue(transform[6].imaginaryPart == 0);

            Assert.IsTrue(Math.Round(transform[7].realPart) == 1);
            Assert.IsTrue(Math.Round(transform[7].imaginaryPart, 2) == 2.41);

            input = transform.InverseFastFourierTransform().ToList();

            Assert.IsTrue(Math.Round(input[0].realPart) == 1);
            Assert.IsTrue(Math.Round(input[0].imaginaryPart) == 0);

            Assert.IsTrue(Math.Round(input[1].realPart) == 1);
            Assert.IsTrue(Math.Round(input[1].imaginaryPart) == 0);

            Assert.IsTrue(Math.Round(input[2].realPart) == 1);
            Assert.IsTrue(Math.Round(input[2].imaginaryPart) == 0);

            Assert.IsTrue(Math.Round(input[3].realPart) == 1);
            Assert.IsTrue(Math.Round(input[3].imaginaryPart) == 0);

            Assert.IsTrue(Math.Round(input[4].realPart) == 0);
            Assert.IsTrue(Math.Round(input[4].imaginaryPart) == 0);

            Assert.IsTrue(Math.Round(input[5].realPart) == 0);
            Assert.IsTrue(Math.Round(input[5].imaginaryPart) == 0);

            Assert.IsTrue(Math.Round(input[6].realPart) == 0);
            Assert.IsTrue(Math.Round(input[6].imaginaryPart) == 0);

            Assert.IsTrue(Math.Round(input[7].realPart) == 0);
            Assert.IsTrue(Math.Round(input[7].imaginaryPart) == 0);
        }

        [TestMethod]
        public void TestDiscreteFourierTransform()
        {
            var input = new List<ComplexNumber>() { new ComplexNumber(1, 0), new ComplexNumber(1, 0), new ComplexNumber(1, 0), new ComplexNumber(1, 0), new ComplexNumber(0, 0), new ComplexNumber(0, 0), new ComplexNumber(0, 0), new ComplexNumber(0, 0) };

            var transform = input.DiscreteFourierTransform(-1).ToList();

            Assert.IsTrue(transform[0].realPart == 4);
            Assert.IsTrue(transform[0].imaginaryPart == 0);

            Assert.IsTrue(transform[1].realPart == 1);
            Assert.IsTrue(Math.Round(transform[1].imaginaryPart, 2) == -2.41);

            Assert.IsTrue(Math.Round(transform[2].realPart) == 0);
            Assert.IsTrue(Math.Round(transform[2].imaginaryPart) == 0);

            Assert.IsTrue(transform[3].realPart == 1);
            Assert.IsTrue(Math.Round(transform[3].imaginaryPart, 2) == -0.41);

            Assert.IsTrue(Math.Round(transform[4].realPart) == 0);
            Assert.IsTrue(Math.Round(transform[4].imaginaryPart) == 0);

            Assert.IsTrue(Math.Round(transform[5].realPart) == 1);
            Assert.IsTrue(Math.Round(transform[5].imaginaryPart, 2) == 0.41);

            Assert.IsTrue(Math.Round(transform[6].realPart) == 0);
            Assert.IsTrue(Math.Round(transform[6].imaginaryPart) == 0);

            Assert.IsTrue(Math.Round(transform[7].realPart) == 1);
            Assert.IsTrue(Math.Round(transform[7].imaginaryPart, 2) == 2.41);
        }
        [TestMethod]
        public void GaussianPulse()
        {
            Func<double, double> func = (double t) => 1 / (4 * Math.Sqrt(2 * Math.PI * 0.01)) * (Math.Exp(-t * t / (2 * 0.01)));
            var timeseries = func.GetSeries(-0.5, 0.5, 100);
            Assert.IsTrue(timeseries.Count() == 100);
            timeseries.Save(@"\timeserie.csv");

            var frequency = func.FastFourierTransform(-0.5, 0.5, 100).ToFrequencyResolution(100);
            Assert.IsTrue(frequency.Count() == 100);
            frequency.Save(@"\frequency.csv");
        }

        [TestMethod]
        public void TestLaplaceTransForm()
        {
            Func<double, double> func = (double t) => 1.0 / Math.Exp(2.0 * t);
            var result = func.LaplaceTransform(2);
            Assert.IsTrue(Math.Round(result, 2) == Math.Round(1.0 / 4.0, 2));
        }

        [TestMethod]
        public void TestInvertLaplaceTransForm()
        {
            Func<double, double> func = (double s) => 1.0 / (s + 2);
            var result = func.InverseLaplaceTransform(3);

            Assert.IsTrue(Math.Round(result, 3) == Math.Round(1.0 / Math.Exp(6), 3));
        }

        [TestMethod]
        public void TestLowPassFilter()
        {
         var rnd = new Random();

            Func<double, double> func = (double t) => Math.Sin(t);
            var series = func.GetSeries(-10, 10, 100);
            var noiseSignal = series.Select(p => p.Value + rnd.GenerateNoise(4)).ToList();
            var output = series.Select(p => p.Value).ToList();
            var result = noiseSignal.LowPassFilter(output, 0.25).ToList();

            Assert.IsTrue(Math.Round(result[32], 3) != Math.Round(noiseSignal[32], 3));

            Assert.IsTrue(Math.Round(result[32], 3) == Math.Round(output[32], 3));
        }
 
    }
}