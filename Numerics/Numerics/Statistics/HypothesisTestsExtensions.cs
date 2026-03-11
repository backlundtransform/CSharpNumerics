using CSharpNumerics.Statistics.Distributions;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Statistics
{
    /// <summary>
    /// Specifies the direction of the alternative hypothesis.
    /// </summary>
    public enum Alternative
    {
        /// <summary>H₁: parameter ≠ value (two-sided test).</summary>
        TwoSided,
        /// <summary>H₁: parameter &lt; value (left-tailed test).</summary>
        Less,
        /// <summary>H₁: parameter &gt; value (right-tailed test).</summary>
        Greater
    }

    /// <summary>
    /// Result of a hypothesis test, including confidence interval and effect size.
    /// </summary>
    public class HypothesisTestResult
    {
        /// <summary>The test statistic (t, z, χ², F, U, W, …).</summary>
        public double TestStatistic { get; }

        /// <summary>The p-value.</summary>
        public double PValue { get; }

        /// <summary>Whether the null hypothesis is rejected at α = 1 − ConfidenceLevel.</summary>
        public bool RejectNull { get; }

        /// <summary>Confidence level used (e.g. 0.95).</summary>
        public double ConfidenceLevel { get; }

        /// <summary>Degrees of freedom (where applicable).</summary>
        public double DegreesOfFreedom { get; }

        /// <summary>Alternative hypothesis direction.</summary>
        public Alternative Alternative { get; }

        /// <summary>Lower bound of the confidence interval for the estimated parameter.</summary>
        public double ConfidenceIntervalLower { get; }

        /// <summary>Upper bound of the confidence interval for the estimated parameter.</summary>
        public double ConfidenceIntervalUpper { get; }

        /// <summary>Effect size (Cohen's d for t-tests, φ for chi-squared, η² for ANOVA).</summary>
        public double EffectSize { get; }

        public HypothesisTestResult(
            double testStatistic, double pValue, double confidenceLevel,
            Alternative alternative, double degreesOfFreedom = 0,
            double ciLower = double.NaN, double ciUpper = double.NaN,
            double effectSize = double.NaN)
        {
            TestStatistic = testStatistic;
            PValue = pValue;
            ConfidenceLevel = confidenceLevel;
            Alternative = alternative;
            DegreesOfFreedom = degreesOfFreedom;
            ConfidenceIntervalLower = ciLower;
            ConfidenceIntervalUpper = ciUpper;
            EffectSize = effectSize;
            RejectNull = pValue < (1.0 - confidenceLevel);
        }

        public override string ToString() =>
            $"Statistic={TestStatistic:F4}, p={PValue:F6}, df={DegreesOfFreedom:F1}, " +
            $"CI=[{ConfidenceIntervalLower:F4}, {ConfidenceIntervalUpper:F4}], " +
            $"d={EffectSize:F4}, Reject={RejectNull} (α={1.0 - ConfidenceLevel:F2}, {Alternative})";
    }

    /// <summary>
    /// Fluent hypothesis-test extensions on <see cref="IEnumerable{double}"/>.
    /// All methods support one-sided and two-sided alternatives.
    /// </summary>
    public static class HypothesisTestsExtensions
    {
        #region One-sample t-test

        /// <summary>
        /// One-sample Student's t-test.
        /// H₀: μ = <paramref name="mu"/>.
        /// </summary>
        /// <param name="sample">Sample values.</param>
        /// <param name="mu">Hypothesized population mean.</param>
        /// <param name="alternative">Direction of the alternative hypothesis.</param>
        /// <param name="confidenceLevel">Confidence level (default 0.95).</param>
        public static HypothesisTestResult TTest(this IEnumerable<double> sample,
            double mu, Alternative alternative = Alternative.TwoSided,
            double confidenceLevel = 0.95)
        {
            var values = sample.ToList();
            int n = values.Count;
            if (n < 2) throw new ArgumentException("Sample must contain at least 2 observations.");

            double mean = values.Average();
            double s = Math.Sqrt(values.Sum(x => (x - mean) * (x - mean)) / (n - 1));
            double se = s / Math.Sqrt(n);
            double t = (mean - mu) / se;
            double df = n - 1;

            var tDist = new StudentTDistribution(df);
            double pValue = ComputePValue(t, tDist, alternative);

            // Confidence interval for the mean
            double tCrit = tDist.InverseCdf(1.0 - (1.0 - confidenceLevel) / 2.0);
            double ciLower = mean - tCrit * se;
            double ciUpper = mean + tCrit * se;

            // Cohen's d
            double d = (mean - mu) / s;

            return new HypothesisTestResult(t, pValue, confidenceLevel, alternative, df, ciLower, ciUpper, d);
        }

        #endregion

        #region Two-sample t-test (Welch)

        /// <summary>
        /// Two-sample Welch's t-test (unequal variances).
        /// H₀: μ₁ = μ₂.
        /// </summary>
        /// <param name="sample1">First sample.</param>
        /// <param name="sample2">Second sample.</param>
        /// <param name="alternative">Direction of the alternative hypothesis.</param>
        /// <param name="confidenceLevel">Confidence level (default 0.95).</param>
        public static HypothesisTestResult TTest(this IEnumerable<double> sample1,
            IEnumerable<double> sample2, Alternative alternative = Alternative.TwoSided,
            double confidenceLevel = 0.95)
        {
            var v1 = sample1.ToList();
            var v2 = sample2.ToList();
            if (v1.Count < 2 || v2.Count < 2)
                throw new ArgumentException("Both samples must contain at least 2 observations.");

            int n1 = v1.Count, n2 = v2.Count;
            double mean1 = v1.Average(), mean2 = v2.Average();
            double s1sq = v1.Sum(x => (x - mean1) * (x - mean1)) / (n1 - 1);
            double s2sq = v2.Sum(x => (x - mean2) * (x - mean2)) / (n2 - 1);

            double se = Math.Sqrt(s1sq / n1 + s2sq / n2);
            double t = (mean1 - mean2) / se;

            // Welch-Satterthwaite degrees of freedom
            double num = Math.Pow(s1sq / n1 + s2sq / n2, 2);
            double den = Math.Pow(s1sq / n1, 2) / (n1 - 1) + Math.Pow(s2sq / n2, 2) / (n2 - 1);
            double df = num / den;

            var tDist = new StudentTDistribution(df);
            double pValue = ComputePValue(t, tDist, alternative);

            double tCrit = tDist.InverseCdf(1.0 - (1.0 - confidenceLevel) / 2.0);
            double diff = mean1 - mean2;
            double ciLower = diff - tCrit * se;
            double ciUpper = diff + tCrit * se;

            // Cohen's d (pooled)
            double sPooled = Math.Sqrt(((n1 - 1) * s1sq + (n2 - 1) * s2sq) / (n1 + n2 - 2));
            double d = sPooled > 0 ? diff / sPooled : 0;

            return new HypothesisTestResult(t, pValue, confidenceLevel, alternative, df, ciLower, ciUpper, d);
        }

        #endregion

        #region Paired t-test

        /// <summary>
        /// Paired t-test on matched observations.
        /// H₀: μ_d = 0.
        /// </summary>
        /// <param name="before">Measurements before treatment.</param>
        /// <param name="after">Measurements after treatment.</param>
        /// <param name="alternative">Direction of the alternative hypothesis.</param>
        /// <param name="confidenceLevel">Confidence level (default 0.95).</param>
        public static HypothesisTestResult PairedTTest(this IEnumerable<double> before,
            IEnumerable<double> after, Alternative alternative = Alternative.TwoSided,
            double confidenceLevel = 0.95)
        {
            var b = before.ToList();
            var a = after.ToList();
            if (b.Count != a.Count)
                throw new ArgumentException("Before and after samples must have equal length.");

            var diffs = new double[b.Count];
            for (int i = 0; i < b.Count; i++)
                diffs[i] = a[i] - b[i];

            return diffs.AsEnumerable().TTest(mu: 0.0, alternative, confidenceLevel);
        }

        #endregion

        #region Z-test

        /// <summary>
        /// One-sample z-test when the population standard deviation is known.
        /// H₀: μ = <paramref name="mu"/>.
        /// </summary>
        /// <param name="sample">Sample values.</param>
        /// <param name="mu">Hypothesized population mean.</param>
        /// <param name="populationStdDev">Known population standard deviation.</param>
        /// <param name="alternative">Direction of the alternative hypothesis.</param>
        /// <param name="confidenceLevel">Confidence level (default 0.95).</param>
        public static HypothesisTestResult ZTest(this IEnumerable<double> sample,
            double mu, double populationStdDev, Alternative alternative = Alternative.TwoSided,
            double confidenceLevel = 0.95)
        {
            var values = sample.ToList();
            int n = values.Count;
            if (n < 1) throw new ArgumentException("Sample must contain at least 1 observation.");

            double mean = values.Average();
            double se = populationStdDev / Math.Sqrt(n);
            double z = (mean - mu) / se;

            var normal = new NormalDistribution();
            double pValue = ComputePValueNormal(z, normal, alternative);

            double zCrit = normal.InverseCdf(1.0 - (1.0 - confidenceLevel) / 2.0);
            double ciLower = mean - zCrit * se;
            double ciUpper = mean + zCrit * se;

            double d = (mean - mu) / populationStdDev;

            return new HypothesisTestResult(z, pValue, confidenceLevel, alternative, ciLower: ciLower, ciUpper: ciUpper, effectSize: d);
        }

        #endregion

        #region Chi-squared goodness-of-fit

        /// <summary>
        /// Chi-squared goodness-of-fit test.
        /// H₀: observed frequencies match expected frequencies.
        /// </summary>
        /// <param name="observed">Observed frequencies.</param>
        /// <param name="expected">Expected frequencies.</param>
        /// <param name="confidenceLevel">Confidence level (default 0.95).</param>
        public static HypothesisTestResult ChiSquaredTest(this IEnumerable<double> observed,
            IEnumerable<double> expected, double confidenceLevel = 0.95)
        {
            var obs = observed.ToList();
            var exp = expected.ToList();
            if (obs.Count != exp.Count)
                throw new ArgumentException("Observed and expected arrays must have equal length.");

            int k = obs.Count;
            double chi2 = 0;
            for (int i = 0; i < k; i++)
            {
                double diff = obs[i] - exp[i];
                chi2 += diff * diff / exp[i];
            }

            double df = k - 1;
            var chi2Dist = new ChiSquaredDistribution(df);
            double pValue = 1.0 - chi2Dist.Cdf(chi2);

            // Cramér's V (for goodness-of-fit, reduces to sqrt(χ²/n))
            double n = obs.Sum();
            double phi = Math.Sqrt(chi2 / n);

            return new HypothesisTestResult(chi2, pValue, confidenceLevel, Alternative.TwoSided, df, effectSize: phi);
        }

        #endregion

        #region F-test for variance equality

        /// <summary>
        /// F-test for equality of two population variances.
        /// H₀: σ₁² = σ₂².
        /// </summary>
        /// <param name="sample1">First sample.</param>
        /// <param name="sample2">Second sample.</param>
        /// <param name="alternative">Direction of the alternative hypothesis.</param>
        /// <param name="confidenceLevel">Confidence level (default 0.95).</param>
        public static HypothesisTestResult FTest(this IEnumerable<double> sample1,
            IEnumerable<double> sample2, Alternative alternative = Alternative.TwoSided,
            double confidenceLevel = 0.95)
        {
            var v1 = sample1.ToList();
            var v2 = sample2.ToList();
            if (v1.Count < 2 || v2.Count < 2)
                throw new ArgumentException("Both samples must contain at least 2 observations.");

            double mean1 = v1.Average(), mean2 = v2.Average();
            double s1sq = v1.Sum(x => (x - mean1) * (x - mean1)) / (v1.Count - 1);
            double s2sq = v2.Sum(x => (x - mean2) * (x - mean2)) / (v2.Count - 1);

            double f = s1sq / s2sq;
            double df1 = v1.Count - 1;
            double df2 = v2.Count - 1;

            var fDist = new FDistribution(df1, df2);
            double pValue;
            switch (alternative)
            {
                case Alternative.Greater:
                    pValue = 1.0 - fDist.Cdf(f);
                    break;
                case Alternative.Less:
                    pValue = fDist.Cdf(f);
                    break;
                default: // TwoSided
                    double pRight = 1.0 - fDist.Cdf(f);
                    double pLeft = fDist.Cdf(f);
                    pValue = 2.0 * Math.Min(pRight, pLeft);
                    break;
            }

            return new HypothesisTestResult(f, pValue, confidenceLevel, alternative, df1,
                effectSize: s1sq / s2sq);
        }

        #endregion

        #region Mann-Whitney U test

        /// <summary>
        /// Mann-Whitney U test (non-parametric alternative to the two-sample t-test).
        /// H₀: the two populations have the same distribution.
        /// Uses normal approximation for n₁, n₂ ≥ 20; otherwise exact for small samples.
        /// </summary>
        /// <param name="sample1">First sample.</param>
        /// <param name="sample2">Second sample.</param>
        /// <param name="alternative">Direction of the alternative hypothesis.</param>
        /// <param name="confidenceLevel">Confidence level (default 0.95).</param>
        public static HypothesisTestResult MannWhitneyUTest(this IEnumerable<double> sample1,
            IEnumerable<double> sample2, Alternative alternative = Alternative.TwoSided,
            double confidenceLevel = 0.95)
        {
            var a = sample1.ToList();
            var b = sample2.ToList();
            int n1 = a.Count, n2 = b.Count;
            if (n1 < 1 || n2 < 1)
                throw new ArgumentException("Both samples must be non-empty.");

            // Combine and rank
            var combined = new List<(double value, int group)>();
            for (int i = 0; i < n1; i++) combined.Add((a[i], 0));
            for (int i = 0; i < n2; i++) combined.Add((b[i], 1));
            combined.Sort((x, y) => x.value.CompareTo(y.value));

            // Assign average ranks
            var ranks = new double[combined.Count];
            int idx = 0;
            while (idx < combined.Count)
            {
                int start = idx;
                while (idx < combined.Count && combined[idx].value == combined[start].value)
                    idx++;
                double avgRank = (start + 1 + idx) / 2.0; // 1-based
                for (int j = start; j < idx; j++)
                    ranks[j] = avgRank;
            }

            double r1 = 0;
            for (int i = 0; i < combined.Count; i++)
                if (combined[i].group == 0)
                    r1 += ranks[i];

            double u1 = r1 - (double)n1 * (n1 + 1) / 2.0;
            double u2 = (double)n1 * n2 - u1;
            double u = Math.Min(u1, u2);

            // Normal approximation
            double mu = (double)n1 * n2 / 2.0;
            double sigma = Math.Sqrt((double)n1 * n2 * (n1 + n2 + 1) / 12.0);
            double z = (u1 - mu) / sigma;

            var normal = new NormalDistribution();
            double pValue = ComputePValueNormal(z, normal, alternative);

            // Effect size: rank-biserial correlation r = 1 - 2U / (n1*n2)
            double r = 1.0 - 2.0 * u / ((double)n1 * n2);

            return new HypothesisTestResult(u1, pValue, confidenceLevel, alternative, effectSize: r);
        }

        #endregion

        #region Wilcoxon signed-rank test

        /// <summary>
        /// Wilcoxon signed-rank test (non-parametric alternative to the one-sample or paired t-test).
        /// H₀: the median of the differences equals zero.
        /// </summary>
        /// <param name="sample">Differences (or paired after − before).</param>
        /// <param name="alternative">Direction of the alternative hypothesis.</param>
        /// <param name="confidenceLevel">Confidence level (default 0.95).</param>
        public static HypothesisTestResult WilcoxonSignedRankTest(this IEnumerable<double> sample,
            Alternative alternative = Alternative.TwoSided,
            double confidenceLevel = 0.95)
        {
            var diffs = sample.Where(d => d != 0).ToList(); // exclude zero diffs
            int n = diffs.Count;
            if (n < 1) throw new ArgumentException("Sample must contain at least 1 non-zero observation.");

            // Sort by absolute value, assign ranks
            var abs = diffs.Select((d, i) => (absVal: Math.Abs(d), sign: Math.Sign(d), idx: i))
                           .OrderBy(x => x.absVal).ToList();

            var ranks = new double[n];
            int k = 0;
            while (k < n)
            {
                int start = k;
                while (k < n && abs[k].absVal == abs[start].absVal)
                    k++;
                double avgRank = (start + 1 + k) / 2.0;
                for (int j = start; j < k; j++)
                    ranks[j] = avgRank;
            }

            double wPlus = 0, wMinus = 0;
            for (int i = 0; i < n; i++)
            {
                if (abs[i].sign > 0) wPlus += ranks[i];
                else wMinus += ranks[i];
            }

            double w = wPlus;

            // Normal approximation
            double mu = (double)n * (n + 1) / 4.0;
            double sigma = Math.Sqrt((double)n * (n + 1) * (2 * n + 1) / 24.0);
            double z = (w - mu) / sigma;

            var normal = new NormalDistribution();
            double pValue = ComputePValueNormal(z, normal, alternative);

            // Effect size: r = z / √n
            double r = z / Math.Sqrt(n);

            return new HypothesisTestResult(w, pValue, confidenceLevel, alternative, effectSize: r);
        }

        #endregion

        #region One-way ANOVA

        /// <summary>
        /// One-way ANOVA: tests whether multiple group means are equal.
        /// H₀: μ₁ = μ₂ = … = μₖ.
        /// </summary>
        /// <param name="groups">Each inner list is one group's observations.</param>
        /// <param name="confidenceLevel">Confidence level (default 0.95).</param>
        public static HypothesisTestResult Anova(this IEnumerable<IEnumerable<double>> groups,
            double confidenceLevel = 0.95)
        {
            var g = groups.Select(x => x.ToList()).ToList();
            int k = g.Count;
            int N = g.Sum(x => x.Count);
            double grandMean = g.SelectMany(x => x).Average();

            double ssBetween = 0;
            double ssWithin = 0;
            foreach (var group in g)
            {
                double groupMean = group.Average();
                ssBetween += group.Count * (groupMean - grandMean) * (groupMean - grandMean);
                ssWithin += group.Sum(x => (x - groupMean) * (x - groupMean));
            }

            double dfBetween = k - 1;
            double dfWithin = N - k;
            double msBetween = ssBetween / dfBetween;
            double msWithin = ssWithin / dfWithin;
            double f = msBetween / msWithin;

            var fDist = new FDistribution(dfBetween, dfWithin);
            double pValue = 1.0 - fDist.Cdf(f);

            // η² = SS_between / SS_total
            double eta2 = ssBetween / (ssBetween + ssWithin);

            return new HypothesisTestResult(f, pValue, confidenceLevel, Alternative.TwoSided, dfBetween,
                effectSize: eta2);
        }

        #endregion

        #region Helpers

        private static double ComputePValue(double t, StudentTDistribution dist, Alternative alt)
        {
            switch (alt)
            {
                case Alternative.Less:
                    return dist.Cdf(t);
                case Alternative.Greater:
                    return 1.0 - dist.Cdf(t);
                default: // TwoSided
                    return 2.0 * (1.0 - dist.Cdf(Math.Abs(t)));
            }
        }

        private static double ComputePValueNormal(double z, NormalDistribution dist, Alternative alt)
        {
            switch (alt)
            {
                case Alternative.Less:
                    return dist.Cdf(z);
                case Alternative.Greater:
                    return 1.0 - dist.Cdf(z);
                default: // TwoSided
                    return 2.0 * (1.0 - dist.Cdf(Math.Abs(z)));
            }
        }

        #endregion
    }
}
