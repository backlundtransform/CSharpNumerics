using CSharpNumerics.Statistics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;

namespace NumericsTests
{
    [TestClass]
    public class HypothesisTestsExtensionsTests
    {
        #region One-sample t-test

        [TestMethod]
        public void TTest_OneSample_TwoSided_KnownResult()
        {
            // Sample with known mean ≈ 5, test H₀: μ = 5
            var sample = new double[] { 4.8, 5.1, 5.3, 4.9, 5.0, 5.2, 4.7, 5.1 };
            var result = sample.TTest(mu: 5.0, alternative: Alternative.TwoSided);

            Assert.IsFalse(result.RejectNull, "Should not reject H₀ when sample mean ≈ μ");
            Assert.IsTrue(result.PValue > 0.05);
            Assert.AreEqual(Alternative.TwoSided, result.Alternative);
            Assert.AreEqual(7, result.DegreesOfFreedom, 1e-10);
        }

        [TestMethod]
        public void TTest_OneSample_RejectsWhenMeanFarFromMu()
        {
            var sample = new double[] { 10.1, 10.3, 9.8, 10.5, 10.0, 10.2, 9.9, 10.4 };
            var result = sample.TTest(mu: 5.0, alternative: Alternative.TwoSided);

            Assert.IsTrue(result.RejectNull, "Should reject H₀ when sample mean far from μ");
            Assert.IsTrue(result.PValue < 0.001);
        }

        [TestMethod]
        public void TTest_OneSample_Greater()
        {
            var sample = new double[] { 11.0, 12.0, 10.5, 11.5, 12.5 };
            var result = sample.TTest(mu: 10.0, alternative: Alternative.Greater);

            Assert.IsTrue(result.RejectNull);
            Assert.AreEqual(Alternative.Greater, result.Alternative);
        }

        [TestMethod]
        public void TTest_OneSample_Less()
        {
            var sample = new double[] { 3.0, 2.5, 3.2, 2.8, 3.1 };
            var result = sample.TTest(mu: 5.0, alternative: Alternative.Less);

            Assert.IsTrue(result.RejectNull);
            Assert.AreEqual(Alternative.Less, result.Alternative);
        }

        [TestMethod]
        public void TTest_OneSample_ConfidenceInterval_BracketsMean()
        {
            var sample = new double[] { 4.8, 5.1, 5.3, 4.9, 5.0, 5.2, 4.7, 5.1 };
            var result = sample.TTest(mu: 5.0);
            double mean = sample.Average();

            Assert.IsTrue(result.ConfidenceIntervalLower <= mean);
            Assert.IsTrue(result.ConfidenceIntervalUpper >= mean);
        }

        [TestMethod]
        public void TTest_OneSample_EffectSize_CohensD()
        {
            var sample = new double[] { 10.1, 10.3, 9.8, 10.5, 10.0 };
            var result = sample.TTest(mu: 5.0);

            // Large effect expected
            Assert.IsTrue(Math.Abs(result.EffectSize) > 1.0);
        }

        [TestMethod]
        public void TTest_OneSample_ToString_IsReadable()
        {
            var sample = new double[] { 5.0, 5.1, 4.9 };
            var result = sample.TTest(mu: 5.0);
            string s = result.ToString();
            Assert.IsTrue(s.Contains("Statistic="));
            Assert.IsTrue(s.Contains("p="));
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void TTest_OneSample_TooFewObservations_Throws()
        {
            new double[] { 5.0 }.TTest(mu: 5.0);
        }

        #endregion

        #region Two-sample t-test

        [TestMethod]
        public void TTest_TwoSample_EqualMeans_NotRejected()
        {
            var s1 = new double[] { 5.1, 4.9, 5.0, 5.2, 4.8 };
            var s2 = new double[] { 5.0, 5.1, 4.9, 5.0, 5.2 };
            var result = s1.TTest(s2, Alternative.TwoSided);

            Assert.IsFalse(result.RejectNull);
        }

        [TestMethod]
        public void TTest_TwoSample_DifferentMeans_Rejected()
        {
            var s1 = new double[] { 10.1, 10.3, 9.8, 10.5, 10.0, 10.2 };
            var s2 = new double[] { 5.1, 4.9, 5.3, 5.0, 5.2, 4.8 };
            var result = s1.TTest(s2);

            Assert.IsTrue(result.RejectNull);
            Assert.IsTrue(result.PValue < 0.001);
        }

        [TestMethod]
        public void TTest_TwoSample_Greater_OneSided()
        {
            var s1 = new double[] { 10.0, 11.0, 10.5, 10.8 };
            var s2 = new double[] { 5.0, 5.5, 5.2, 5.1 };
            var result = s1.TTest(s2, Alternative.Greater);

            var resultTwoSided = s1.TTest(s2, Alternative.TwoSided);

            Assert.IsTrue(result.RejectNull);
            Assert.IsTrue(result.PValue < resultTwoSided.PValue);
        }

        [TestMethod]
        public void TTest_TwoSample_WelchDf_NotInteger()
        {
            var s1 = new double[] { 10.0, 11.0, 10.5, 10.8, 9.5 };
            var s2 = new double[] { 5.0, 5.5 };
            var result = s1.TTest(s2);

            // Welch df is not necessarily integer
            Assert.IsTrue(result.DegreesOfFreedom > 0);
        }

        #endregion

        #region Paired t-test

        [TestMethod]
        public void PairedTTest_SignificantDifference()
        {
            var before = new double[] { 60, 62, 58, 65, 61 };
            var after = new double[] { 70, 73, 68, 75, 72 };
            var result = before.PairedTTest(after);

            Assert.IsTrue(result.RejectNull);
        }

        [TestMethod]
        public void PairedTTest_NoSignificantDifference()
        {
            var before = new double[] { 60, 62, 58, 65, 61 };
            var after = new double[] { 61, 61, 59, 64, 62 };
            var result = before.PairedTTest(after);

            Assert.IsFalse(result.RejectNull);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void PairedTTest_UnequalLengths_Throws()
        {
            new double[] { 1, 2, 3 }.PairedTTest(new double[] { 1, 2 });
        }

        #endregion

        #region Z-test

        [TestMethod]
        public void ZTest_NotRejected_WhenMeanCloseToMu()
        {
            var sample = new double[] { 100.1, 99.8, 100.3, 99.9, 100.0 };
            var result = sample.ZTest(mu: 100.0, populationStdDev: 0.5);

            Assert.IsFalse(result.RejectNull);
        }

        [TestMethod]
        public void ZTest_Rejected_WhenMeanFarFromMu()
        {
            var sample = new double[] { 105.0, 104.5, 105.2, 104.8, 105.1 };
            var result = sample.ZTest(mu: 100.0, populationStdDev: 0.5);

            Assert.IsTrue(result.RejectNull);
            Assert.IsTrue(result.PValue < 0.001);
        }

        [TestMethod]
        public void ZTest_OneSided_Greater()
        {
            var sample = new double[] { 105.0, 104.5, 105.2, 104.8, 105.1 };
            var resultTwoSided = sample.ZTest(mu: 100.0, populationStdDev: 0.5, alternative: Alternative.TwoSided);
            var resultGreater = sample.ZTest(mu: 100.0, populationStdDev: 0.5, alternative: Alternative.Greater);

            // One-sided p-value should be half of two-sided
            Assert.AreEqual(resultTwoSided.PValue / 2, resultGreater.PValue, 1e-10);
        }

        #endregion

        #region Chi-squared test

        [TestMethod]
        public void ChiSquaredTest_GoodFit_NotRejected()
        {
            // Fair die: 100 rolls, roughly uniform
            var observed = new double[] { 18, 16, 17, 15, 17, 17 };
            var expected = new double[] { 16.67, 16.67, 16.67, 16.67, 16.67, 16.67 };
            var result = observed.ChiSquaredTest(expected);

            Assert.IsFalse(result.RejectNull);
        }

        [TestMethod]
        public void ChiSquaredTest_BadFit_Rejected()
        {
            var observed = new double[] { 50, 10, 10, 10, 10, 10 };
            var expected = new double[] { 16.67, 16.67, 16.67, 16.67, 16.67, 16.67 };
            var result = observed.ChiSquaredTest(expected);

            Assert.IsTrue(result.RejectNull);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void ChiSquaredTest_UnequalLengths_Throws()
        {
            new double[] { 10, 20 }.ChiSquaredTest(new double[] { 15 });
        }

        #endregion

        #region F-test

        [TestMethod]
        public void FTest_EqualVariances_NotRejected()
        {
            var s1 = new double[] { 10.1, 10.3, 9.8, 10.5, 10.0, 10.2 };
            var s2 = new double[] { 10.2, 10.0, 10.4, 9.9, 10.1, 10.3 };
            var result = s1.FTest(s2);

            Assert.IsFalse(result.RejectNull);
        }

        [TestMethod]
        public void FTest_UnequalVariances_Rejected()
        {
            var s1 = new double[] { 10.0, 20.0, 5.0, 25.0, 15.0, 0.0 }; // high variance
            var s2 = new double[] { 10.0, 10.1, 10.0, 10.1, 10.0, 10.1 }; // low variance
            var result = s1.FTest(s2, Alternative.TwoSided);

            Assert.IsTrue(result.RejectNull);
        }

        #endregion

        #region Mann-Whitney U test

        [TestMethod]
        public void MannWhitneyU_SameDistribution_NotRejected()
        {
            var s1 = new double[] { 1.1, 2.3, 3.0, 4.2, 5.1, 6.0, 7.3, 8.1, 9.2, 10.0,
                                    11.1, 12.0, 13.2, 14.1, 15.0, 16.3, 17.1, 18.0, 19.2, 20.1 };
            var s2 = new double[] { 1.5, 2.7, 3.3, 4.5, 5.4, 6.3, 7.6, 8.4, 9.5, 10.3,
                                    11.4, 12.3, 13.5, 14.4, 15.3, 16.6, 17.4, 18.3, 19.5, 20.4 };
            var result = s1.MannWhitneyUTest(s2);

            Assert.IsFalse(result.RejectNull);
        }

        [TestMethod]
        public void MannWhitneyU_DifferentDistributions_Rejected()
        {
            var s1 = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                    11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };
            var s2 = new double[] { 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                                    60, 61, 62, 63, 64, 65, 66, 67, 68, 69 };
            var result = s1.MannWhitneyUTest(s2);

            Assert.IsTrue(result.RejectNull);
        }

        #endregion

        #region Wilcoxon signed-rank test

        [TestMethod]
        public void WilcoxonSignedRank_SignificantDifference()
        {
            // All positive differences → reject H₀: median = 0
            var diffs = new double[] { 3.0, 4.5, 2.1, 5.3, 3.8, 4.0, 2.9, 5.1, 3.5, 4.2 };
            var result = diffs.WilcoxonSignedRankTest();

            Assert.IsTrue(result.RejectNull);
        }

        [TestMethod]
        public void WilcoxonSignedRank_NoSignificantDifference()
        {
            // Balanced around 0
            var diffs = new double[] { 0.1, -0.2, 0.15, -0.1, 0.05, -0.15, 0.2, -0.05, 0.1, -0.1 };
            var result = diffs.WilcoxonSignedRankTest();

            Assert.IsFalse(result.RejectNull);
        }

        #endregion

        #region One-way ANOVA

        [TestMethod]
        public void Anova_EqualMeans_NotRejected()
        {
            var groups = new List<IEnumerable<double>>
            {
                new double[] { 5.1, 4.9, 5.0, 5.2, 4.8, 5.1 },
                new double[] { 5.0, 5.1, 4.9, 5.0, 5.2, 5.0 },
                new double[] { 4.9, 5.0, 5.1, 4.8, 5.2, 5.0 }
            };
            var result = groups.Anova();

            Assert.IsFalse(result.RejectNull);
        }

        [TestMethod]
        public void Anova_DifferentMeans_Rejected()
        {
            var groups = new List<IEnumerable<double>>
            {
                new double[] { 2.0, 2.1, 1.9, 2.2, 1.8 },
                new double[] { 5.0, 5.1, 4.9, 5.2, 4.8 },
                new double[] { 8.0, 8.1, 7.9, 8.2, 7.8 }
            };
            var result = groups.Anova();

            Assert.IsTrue(result.RejectNull);
            Assert.IsTrue(result.PValue < 0.001);
            Assert.IsTrue(result.EffectSize > 0.9); // very large η²
        }

        #endregion
    }
}
