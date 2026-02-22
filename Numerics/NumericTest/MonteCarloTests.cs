using CSharpNumerics.Statistics.Random;
using CSharpNumerics.Statistics.Distributions;
using CSharpNumerics.Statistics.MonteCarlo;

namespace NumericsTests
{
    [TestClass]
    public class RandomGeneratorTests
    {
        [TestMethod]
        public void NextGaussian_MeanApproximatesZero()
        {
            var rng = new RandomGenerator(42);
            double sum = 0;
            int n = 100_000;
            for (int i = 0; i < n; i++)
                sum += rng.NextGaussian();

            double mean = sum / n;
            Assert.IsTrue(Math.Abs(mean) < 0.02, $"Mean was {mean}, expected near 0");
        }

        [TestMethod]
        public void NextGaussian_WithParameters()
        {
            var rng = new RandomGenerator(42);
            double mu = 10.0;
            double sigma = 2.0;
            double sum = 0;
            int n = 100_000;
            for (int i = 0; i < n; i++)
                sum += rng.NextGaussian(mu, sigma);

            double mean = sum / n;
            Assert.IsTrue(Math.Abs(mean - mu) < 0.05, $"Mean was {mean}, expected near {mu}");
        }

        [TestMethod]
        public void NextUniform_RespectsRange()
        {
            var rng = new RandomGenerator(7);
            for (int i = 0; i < 10_000; i++)
            {
                double val = rng.NextUniform(2.0, 5.0);
                Assert.IsTrue(val >= 2.0 && val < 5.0, $"Value {val} out of range [2, 5)");
            }
        }

        [TestMethod]
        public void NextExponential_MeanApproximatesInverseLambda()
        {
            var rng = new RandomGenerator(99);
            double lambda = 3.0;
            double sum = 0;
            int n = 100_000;
            for (int i = 0; i < n; i++)
                sum += rng.NextExponential(lambda);

            double mean = sum / n;
            Assert.IsTrue(Math.Abs(mean - 1.0 / lambda) < 0.02, $"Mean was {mean}, expected near {1.0/lambda}");
        }

        [TestMethod]
        public void NextPoisson_MeanApproximatesLambda()
        {
            var rng = new RandomGenerator(55);
            double lambda = 4.0;
            double sum = 0;
            int n = 100_000;
            for (int i = 0; i < n; i++)
                sum += rng.NextPoisson(lambda);

            double mean = sum / n;
            Assert.IsTrue(Math.Abs(mean - lambda) < 0.1, $"Mean was {mean}, expected near {lambda}");
        }

        [TestMethod]
        public void NextBernoulli_MeanApproximatesP()
        {
            var rng = new RandomGenerator(12);
            double p = 0.3;
            double sum = 0;
            int n = 100_000;
            for (int i = 0; i < n; i++)
                sum += rng.NextBernoulli(p);

            double mean = sum / n;
            Assert.IsTrue(Math.Abs(mean - p) < 0.02, $"Mean was {mean}, expected near {p}");
        }

        [TestMethod]
        public void NextBinomial_MeanApproximatesNP()
        {
            var rng = new RandomGenerator(21);
            int trials = 20;
            double p = 0.4;
            double sum = 0;
            int n = 100_000;
            for (int i = 0; i < n; i++)
                sum += rng.NextBinomial(trials, p);

            double mean = sum / n;
            double expected = trials * p;
            Assert.IsTrue(Math.Abs(mean - expected) < 0.1, $"Mean was {mean}, expected near {expected}");
        }

        [TestMethod]
        public void Shuffle_PreservesElements()
        {
            var rng = new RandomGenerator(33);
            var arr = new[] { 1, 2, 3, 4, 5 };
            rng.Shuffle(arr);

            Array.Sort(arr);
            CollectionAssert.AreEqual(new[] { 1, 2, 3, 4, 5 }, arr);
        }

        [TestMethod]
        public void Sample_ReturnsCorrectCount()
        {
            var rng = new RandomGenerator(77);
            var source = new[] { 10.0, 20.0, 30.0, 40.0, 50.0 };
            var sampled = rng.Sample(source, 3);
            Assert.AreEqual(3, sampled.Length);
        }

        [TestMethod]
        public void GaussianSamples_ReturnsCorrectCount()
        {
            var rng = new RandomGenerator(88);
            var samples = rng.GaussianSamples(500);
            Assert.AreEqual(500, samples.Length);
        }

        [TestMethod]
        public void Reproducibility_SameSeedSameResults()
        {
            var rng1 = new RandomGenerator(42);
            var rng2 = new RandomGenerator(42);

            for (int i = 0; i < 100; i++)
                Assert.AreEqual(rng1.NextGaussian(), rng2.NextGaussian(), 1e-15);
        }
    }

    [TestClass]
    public class DistributionTests
    {
        [TestMethod]
        public void NormalDistribution_PdfPeakAtMean()
        {
            var dist = new NormalDistribution(5.0, 2.0);
            double pdfAtMean = dist.Pdf(5.0);
            double pdfAway = dist.Pdf(10.0);

            Assert.IsTrue(pdfAtMean > pdfAway, "PDF should peak at the mean");
        }

        [TestMethod]
        public void NormalDistribution_CdfAtMeanIsHalf()
        {
            var dist = new NormalDistribution(0.0, 1.0);
            double cdf = dist.Cdf(0.0);
            Assert.IsTrue(Math.Abs(cdf - 0.5) < 0.001, $"CDF(0) = {cdf}, expected 0.5");
        }

        [TestMethod]
        public void NormalDistribution_CdfAtTwoSigma()
        {
            var dist = new NormalDistribution(0.0, 1.0);
            double cdf = dist.Cdf(2.0);
            // Expected ≈ 0.9772
            Assert.IsTrue(Math.Abs(cdf - 0.9772) < 0.001, $"CDF(2) = {cdf}, expected ~0.9772");
        }

        [TestMethod]
        public void NormalDistribution_InverseCdfRoundTrips()
        {
            var dist = new NormalDistribution(3.0, 1.5);
            double x = 4.5;
            double p = dist.Cdf(x);
            double xBack = dist.InverseCdf(p);
            Assert.IsTrue(Math.Abs(x - xBack) < 0.01, $"InverseCdf round-trip: {x} -> {p} -> {xBack}");
        }

        [TestMethod]
        public void NormalDistribution_SampleMean()
        {
            var dist = new NormalDistribution(10.0, 3.0);
            var rng = new RandomGenerator(42);
            var samples = dist.Samples(rng, 100_000);
            double mean = samples.Average();
            Assert.IsTrue(Math.Abs(mean - 10.0) < 0.1, $"Sample mean = {mean}, expected ~10.0");
        }

        [TestMethod]
        public void UniformDistribution_PdfConstant()
        {
            var dist = new UniformDistribution(2.0, 5.0);
            double expected = 1.0 / 3.0;

            Assert.IsTrue(Math.Abs(dist.Pdf(3.0) - expected) < 1e-10);
            Assert.IsTrue(Math.Abs(dist.Pdf(4.5) - expected) < 1e-10);
            Assert.AreEqual(0.0, dist.Pdf(1.0));
            Assert.AreEqual(0.0, dist.Pdf(6.0));
        }

        [TestMethod]
        public void UniformDistribution_CdfLinear()
        {
            var dist = new UniformDistribution(0.0, 10.0);
            Assert.IsTrue(Math.Abs(dist.Cdf(5.0) - 0.5) < 1e-10);
            Assert.IsTrue(Math.Abs(dist.Cdf(0.0) - 0.0) < 1e-10);
            Assert.IsTrue(Math.Abs(dist.Cdf(10.0) - 1.0) < 1e-10);
        }

        [TestMethod]
        public void ExponentialDistribution_PdfAndCdf()
        {
            var dist = new ExponentialDistribution(2.0);
            // PDF(0) = λ = 2
            Assert.IsTrue(Math.Abs(dist.Pdf(0.0) - 2.0) < 1e-10);
            // CDF(0) = 0
            Assert.IsTrue(Math.Abs(dist.Cdf(0.0)) < 1e-10);
            // CDF(∞) → 1, check at large value
            Assert.IsTrue(dist.Cdf(10.0) > 0.999);
        }

        [TestMethod]
        public void ExponentialDistribution_MeanAndVariance()
        {
            var dist = new ExponentialDistribution(0.5);
            Assert.IsTrue(Math.Abs(dist.Mean - 2.0) < 1e-10);
            Assert.IsTrue(Math.Abs(dist.Variance - 4.0) < 1e-10);
        }

        [TestMethod]
        public void PoissonDistribution_Pmf()
        {
            var dist = new PoissonDistribution(3.0);
            // P(X=0) = e^(-3) ≈ 0.04979
            Assert.IsTrue(Math.Abs(dist.Pdf(0) - Math.Exp(-3)) < 1e-6);
            // P(X=3) = e^(-3)*3^3/6 ≈ 0.22404
            double expected = Math.Exp(-3) * 27.0 / 6.0;
            Assert.IsTrue(Math.Abs(dist.Pdf(3) - expected) < 1e-6);
        }

        [TestMethod]
        public void PoissonDistribution_CdfSumsToOne()
        {
            var dist = new PoissonDistribution(5.0);
            // CDF at large value should be very close to 1
            Assert.IsTrue(dist.Cdf(30) > 0.9999);
        }

        [TestMethod]
        public void BernoulliDistribution_Pmf()
        {
            var dist = new BernoulliDistribution(0.7);
            Assert.IsTrue(Math.Abs(dist.Pdf(1.0) - 0.7) < 1e-10);
            Assert.IsTrue(Math.Abs(dist.Pdf(0.0) - 0.3) < 1e-10);
            Assert.AreEqual(0.0, dist.Pdf(0.5));
        }

        [TestMethod]
        public void BinomialDistribution_MeanAndVariance()
        {
            var dist = new BinomialDistribution(10, 0.3);
            Assert.IsTrue(Math.Abs(dist.Mean - 3.0) < 1e-10);
            Assert.IsTrue(Math.Abs(dist.Variance - 2.1) < 1e-10);
        }

        [TestMethod]
        public void BinomialDistribution_PmfSumsToOne()
        {
            var dist = new BinomialDistribution(8, 0.5);
            double sum = 0;
            for (int k = 0; k <= 8; k++)
                sum += dist.Pdf(k);
            Assert.IsTrue(Math.Abs(sum - 1.0) < 1e-10, $"PMF sum = {sum}");
        }

        [TestMethod]
        public void BinomialDistribution_SampleMean()
        {
            var dist = new BinomialDistribution(20, 0.4);
            var rng = new RandomGenerator(42);
            var samples = dist.Samples(rng, 100_000);
            double mean = samples.Average();
            Assert.IsTrue(Math.Abs(mean - 8.0) < 0.1, $"Sample mean = {mean}, expected ~8.0");
        }
    }

    [TestClass]
    public class MonteCarloTests
    {
        [TestMethod]
        public void EstimatePi_WithInlineFunction()
        {
            // Classic Monte Carlo: hit-or-miss π estimation
            var sim = new MonteCarloSimulator(42);
            var result = sim.Run(rng =>
            {
                double x = rng.NextUniform(0, 1);
                double y = rng.NextUniform(0, 1);
                return x * x + y * y <= 1.0 ? 1.0 : 0.0;
            }, iterations: 500_000);

            double piEstimate = result.Mean * 4.0;
            Assert.IsTrue(Math.Abs(piEstimate - Math.PI) < 0.02,
                $"π estimate = {piEstimate}, expected ~{Math.PI}");
        }

        [TestMethod]
        public void IntegrateUniform_SinX()
        {
            // ∫₀^π sin(x) dx = 2
            var sim = new MonteCarloSimulator(99);
            var result = sim.IntegrateUniform(Math.Sin, 0, Math.PI, 500_000);

            Assert.IsTrue(Math.Abs(result.Mean - 2.0) < 0.02,
                $"Integral estimate = {result.Mean}, expected ~2.0");
        }

        [TestMethod]
        public void Integrate_WithDistribution()
        {
            // E[X²] where X ~ N(0,1) should be ≈ 1 (variance of standard normal)
            var sim = new MonteCarloSimulator(55);
            var normal = new NormalDistribution(0, 1);
            var result = sim.Integrate(x => x * x, normal, 200_000);

            Assert.IsTrue(Math.Abs(result.Mean - 1.0) < 0.02,
                $"E[X²] = {result.Mean}, expected ~1.0");
        }

        [TestMethod]
        public void Run_WithModel_DiceSum()
        {
            // Sum of two dice: E[X] = 7
            var sim = new MonteCarloSimulator(42);
            var model = new TwoDiceModel();
            var result = sim.Run(model, 200_000);

            Assert.IsTrue(Math.Abs(result.Mean - 7.0) < 0.05,
                $"Mean dice sum = {result.Mean}, expected ~7.0");
        }

        [TestMethod]
        public void Result_ConfidenceInterval()
        {
            var sim = new MonteCarloSimulator(42);
            var result = sim.Run(rng => rng.NextGaussian(100, 15), 100_000);

            var (lower, upper) = result.ConfidenceInterval(0.95);
            Assert.IsTrue(lower < 100 && upper > 100,
                $"95% CI [{lower}, {upper}] should contain 100");
            Assert.IsTrue(upper - lower < 70,
                "95% CI should be reasonable width for N(100,15)");
        }

        [TestMethod]
        public void Result_Probability()
        {
            var sim = new MonteCarloSimulator(42);
            var result = sim.Run(rng => rng.NextGaussian(0, 1), 200_000);

            // P(X > 0) should be ≈ 0.5
            double pPositive = result.Probability(x => x > 0);
            Assert.IsTrue(Math.Abs(pPositive - 0.5) < 0.01,
                $"P(X>0) = {pPositive}, expected ~0.5");
        }

        [TestMethod]
        public void Result_Histogram_BinsCorrectly()
        {
            var sim = new MonteCarloSimulator(42);
            var result = sim.Run(rng => rng.NextUniform(0, 1), 10_000);

            var histogram = result.Histogram(10);
            Assert.AreEqual(10, histogram.Length);

            int total = 0;
            foreach (var (_, count) in histogram)
                total += count;
            Assert.AreEqual(10_000, total, "Histogram should contain all samples");
        }

        [TestMethod]
        public void Result_Percentile()
        {
            var sim = new MonteCarloSimulator(42);
            var result = sim.Run(rng => rng.NextUniform(0, 100), 100_000);

            double p50 = result.Percentile(50);
            Assert.IsTrue(Math.Abs(p50 - 50.0) < 1.0,
                $"50th percentile = {p50}, expected ~50");

            double p90 = result.Percentile(90);
            Assert.IsTrue(Math.Abs(p90 - 90.0) < 1.0,
                $"90th percentile = {p90}, expected ~90");
        }

        [TestMethod]
        public void RunMultivariate_TwoDimensions()
        {
            var sim = new MonteCarloSimulator(42);
            var results = sim.RunMultivariate(rng => new[]
            {
                rng.NextGaussian(0, 1),
                rng.NextGaussian(10, 2)
            }, 100_000);

            Assert.AreEqual(2, results.Length);
            Assert.IsTrue(Math.Abs(results[0].Mean) < 0.05, $"Dim0 mean = {results[0].Mean}");
            Assert.IsTrue(Math.Abs(results[1].Mean - 10.0) < 0.05, $"Dim1 mean = {results[1].Mean}");
        }

        [TestMethod]
        public void ConvergenceCurve_ApproachesTrue()
        {
            var sim = new MonteCarloSimulator(42);
            var curve = sim.ConvergenceCurve(rng => rng.NextGaussian(5.0, 1.0), 50_000);

            Assert.AreEqual(50_000, curve.Length);
            // Last value (after many samples) should be close to true mean
            Assert.IsTrue(Math.Abs(curve[curve.Length - 1] - 5.0) < 0.05,
                $"Final running mean = {curve[curve.Length - 1]}, expected ~5.0");
        }

        [TestMethod]
        public void StandardError_DecreasesWithMoreSamples()
        {
            var sim = new MonteCarloSimulator(42);
            var smallResult = sim.Run(rng => rng.NextGaussian(), 1_000);
            var largeResult = sim.Run(rng => rng.NextGaussian(), 100_000);

            // Note: different seeds, but SEM ∝ 1/√n so large should have smaller SE
            // We check the property holds structurally
            Assert.IsTrue(largeResult.StandardError < smallResult.StandardError * 0.5,
                $"Large SE={largeResult.StandardError}, Small SE={smallResult.StandardError}");
        }

        /// <summary>
        /// A simple model: sum of two fair dice.
        /// </summary>
        private class TwoDiceModel : IMonteCarloModel
        {
            public double Evaluate(RandomGenerator rng)
            {
                return rng.NextInt(1, 7) + rng.NextInt(1, 7);
            }
        }
    }
}
