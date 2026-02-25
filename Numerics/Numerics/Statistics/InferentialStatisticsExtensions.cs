using CSharpNumerics.Statistics.Distributions;
using System;
using System.Collections.Generic;
using System.Linq;



namespace CSharpNumerics.Statistics;


/// <summary>
/// Represents the result of a statistical hypothesis test.
/// </summary>
public class StatisticalTestResult
{
    /// <summary>The test statistic (t, z, χ², F, etc.).</summary>
    public double TestStatistic { get; }

    /// <summary>The p-value of the test.</summary>
    public double PValue { get; }

    /// <summary>Whether the null hypothesis is rejected at the given confidence level.</summary>
    public bool RejectNull { get; }

    /// <summary>The confidence level used (e.g. 0.95).</summary>
    public double ConfidenceLevel { get; }

    /// <summary>Degrees of freedom (where applicable).</summary>
    public double DegreesOfFreedom { get; }

    public StatisticalTestResult(double testStatistic, double pValue, double confidenceLevel, double degreesOfFreedom = 0)
    {
        TestStatistic = testStatistic;
        PValue = pValue;
        ConfidenceLevel = confidenceLevel;
        DegreesOfFreedom = degreesOfFreedom;
        RejectNull = pValue < (1.0 - confidenceLevel);
    }

    public override string ToString() =>
        $"Statistic={TestStatistic:F4}, p={PValue:F6}, df={DegreesOfFreedom:F1}, Reject={RejectNull} (α={1.0 - ConfidenceLevel:F2})";
}

public static class InferentialStatisticsExtensions
{
    // ──────────────────────────────────────────────
    //  REGRESSION
    // ──────────────────────────────────────────────

    /// <summary>
    /// Fits a simple linear regression y = slope * x + intercept.
    /// Returns slope, intercept and Pearson correlation coefficient r.
    /// </summary>
    public static (double slope, double intercept, double correlation) LinearRegression<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y)> func)
    {
        var serie = enumerable.Select(func);
        var sxx = serie.Sum(p => Math.Pow(p.x, 2)) - 1.0 / serie.Count() * Math.Pow(serie.Sum(p => p.x), 2);
        var syy = serie.Sum(p => Math.Pow(p.y, 2)) - 1.0 / serie.Count() * Math.Pow(serie.Sum(p => p.y), 2);
        var sxy = serie.Sum(p => p.x * p.y) - 1.0 / serie.Count() * (serie.Sum(p => p.x) * serie.Sum(p => p.y));

        var slope = sxy / sxx;
        var intercept = serie.Average(p => p.y) - slope * serie.Average(p => p.x);
        var correlation = sxy / Math.Sqrt(sxx * syy);

        return (slope, intercept, correlation);
    }

    /// <summary>
    /// Fits an exponential model y = exp(intercept) * exp(slope * x) by applying linear regression to log(y).
    /// </summary>
    /// <typeparam name="T">The element type.</typeparam>
    /// <param name="enumerable">The input sequence.</param>
    /// <param name="func">Projection from element to (x,y).</param>
    /// <returns>A function f(x) representing the fitted exponential curve.</returns>
    public static Func<double, double> ExponentialRegression<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y)> func)
    {
        var (slope, intercept, correlation) = enumerable.Select(func).LinearRegression(p => (p.x, Math.Log(p.y)));
        return (double x) => Math.Exp(intercept) * Math.Exp(slope * x);
    }

    /// <summary>
    /// Fits a logarithmic model y = a + b * ln(x) by applying linear regression to (ln(x), y).
    /// </summary>
    /// <typeparam name="T">The element type.</typeparam>
    /// <param name="enumerable">The input sequence.</param>
    /// <param name="func">Projection from element to (x,y). x must be > 0.</param>
    /// <returns>A function f(x) representing the fitted logarithmic curve.</returns>
    public static Func<double, double> LogarithmicRegression<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y)> func)
    {
        var (slope, intercept, _) = enumerable.Select(func).LinearRegression(p => (Math.Log(p.x), p.y));
        return (double x) => intercept + slope * Math.Log(x);
    }

    /// <summary>
    /// Fits a power model y = a * x^b by applying linear regression to (ln(x), ln(y)).
    /// </summary>
    /// <typeparam name="T">The element type.</typeparam>
    /// <param name="enumerable">The input sequence.</param>
    /// <param name="func">Projection from element to (x,y). Both x and y must be > 0.</param>
    /// <returns>A function f(x) representing the fitted power curve.</returns>
    public static Func<double, double> PowerRegression<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y)> func)
    {
        var (slope, intercept, _) = enumerable.Select(func).LinearRegression(p => (Math.Log(p.x), Math.Log(p.y)));
        double a = Math.Exp(intercept);
        return (double x) => a * Math.Pow(x, slope);
    }

    /// <summary>
    /// Fits a polynomial of the specified degree using ordinary least squares (normal equations).
    /// Returns the polynomial as a function and the coefficients in ascending order [a0, a1, ..., an].
    /// </summary>
    /// <typeparam name="T">The element type.</typeparam>
    /// <param name="enumerable">The input sequence.</param>
    /// <param name="func">Projection from element to (x,y).</param>
    /// <param name="degree">Polynomial degree (1 = linear, 2 = quadratic, etc.).</param>
    /// <returns>Tuple of (fitted function, coefficients array).</returns>
    public static (Func<double, double> predict, double[] coefficients) PolynomialRegression<T>(
        this IEnumerable<T> enumerable, Func<T, (double x, double y)> func, int degree)
    {
        if (degree < 1)
            throw new ArgumentException("Degree must be at least 1.");

        var data = enumerable.Select(func).ToList();
        int n = data.Count;
        int m = degree + 1;

        // Build the Vandermonde-style normal equations: (XᵀX) β = Xᵀy
        double[] sums = new double[2 * degree + 1];
        for (int i = 0; i < sums.Length; i++)
            sums[i] = data.Sum(p => Math.Pow(p.x, i));

        double[] rhs = new double[m];
        for (int i = 0; i < m; i++)
            rhs[i] = data.Sum(p => Math.Pow(p.x, i) * p.y);

        // Solve via Gaussian elimination
        double[,] matrix = new double[m, m + 1];
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < m; j++)
                matrix[i, j] = sums[i + j];
            matrix[i, m] = rhs[i];
        }

        double[] coefficients = GaussianElimination(matrix, m);

        Func<double, double> predict = (double x) =>
        {
            double result = 0;
            for (int i = 0; i < coefficients.Length; i++)
                result += coefficients[i] * Math.Pow(x, i);
            return result;
        };

        return (predict, coefficients);
    }

    // ──────────────────────────────────────────────
    //  CORRELATION
    // ──────────────────────────────────────────────

    /// <summary>
    /// Computes the Pearson correlation coefficient r and its two-tailed p-value.
    /// </summary>
    /// <typeparam name="T">The element type.</typeparam>
    /// <param name="enumerable">The input sequence.</param>
    /// <param name="func">Projection from element to (x,y).</param>
    /// <returns>(r, pValue)</returns>
    public static (double r, double pValue) PearsonCorrelation<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y)> func)
    {
        var data = enumerable.Select(func).ToList();
        int n = data.Count;

        double meanX = data.Average(p => p.x);
        double meanY = data.Average(p => p.y);

        double sxy = data.Sum(p => (p.x - meanX) * (p.y - meanY));
        double sxx = data.Sum(p => Math.Pow(p.x - meanX, 2));
        double syy = data.Sum(p => Math.Pow(p.y - meanY, 2));

        double r = sxy / Math.Sqrt(sxx * syy);

        // t-statistic for testing H₀: ρ = 0
        double t = r * Math.Sqrt((n - 2) / (1.0 - r * r));
        var tDist = new StudentTDistribution(n - 2);
        double pValue = tDist.TwoTailedPValue(t);

        return (r, pValue);
    }

    /// <summary>
    /// Computes the Spearman rank correlation coefficient ρ and its two-tailed p-value.
    /// </summary>
    /// <typeparam name="T">The element type.</typeparam>
    /// <param name="enumerable">The input sequence.</param>
    /// <param name="func">Projection from element to (x,y).</param>
    /// <returns>(rho, pValue)</returns>
    public static (double rho, double pValue) SpearmanCorrelation<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y)> func)
    {
        var data = enumerable.Select(func).ToList();
        int n = data.Count;

        var rankX = ComputeRanks(data.Select(p => p.x).ToArray());
        var rankY = ComputeRanks(data.Select(p => p.y).ToArray());

        // Pearson on the ranks
        var ranked = rankX.Zip(rankY, (rx, ry) => (x: rx, y: ry)).ToList();
        double meanRx = ranked.Average(p => p.x);
        double meanRy = ranked.Average(p => p.y);

        double sxy = ranked.Sum(p => (p.x - meanRx) * (p.y - meanRy));
        double sxx = ranked.Sum(p => Math.Pow(p.x - meanRx, 2));
        double syy = ranked.Sum(p => Math.Pow(p.y - meanRy, 2));

        double rho = sxy / Math.Sqrt(sxx * syy);

        double t = rho * Math.Sqrt((n - 2) / (1.0 - rho * rho));
        var tDist = new StudentTDistribution(n - 2);
        double pValue = tDist.TwoTailedPValue(t);

        return (rho, pValue);
    }

    // ──────────────────────────────────────────────
    //  HYPOTHESIS TESTS
    // ──────────────────────────────────────────────

    /// <summary>
    /// One-sample t-test: tests whether the population mean equals a hypothesized value.
    /// H₀: μ = hypothesizedMean.
    /// </summary>
    /// <typeparam name="T">The element type.</typeparam>
    /// <param name="enumerable">The input sequence.</param>
    /// <param name="func">Projection from element to numeric value.</param>
    /// <param name="hypothesizedMean">The value to test against.</param>
    /// <param name="confidenceLevel">Confidence level (e.g. 0.95).</param>
    public static StatisticalTestResult OneSampleTTest<T>(this IEnumerable<T> enumerable, Func<T, double> func,
        double hypothesizedMean, double confidenceLevel = 0.95)
    {
        var values = enumerable.Select(func).ToList();
        int n = values.Count;
        double mean = values.Average();
        double s = Math.Sqrt(values.Sum(x => Math.Pow(x - mean, 2)) / (n - 1));
        double t = (mean - hypothesizedMean) / (s / Math.Sqrt(n));
        double df = n - 1;

        var tDist = new StudentTDistribution(df);
        double pValue = tDist.TwoTailedPValue(t);

        return new StatisticalTestResult(t, pValue, confidenceLevel, df);
    }

    /// <summary>
    /// Two-sample independent t-test (Welch's): tests whether two population means are equal.
    /// H₀: μ₁ = μ₂.
    /// </summary>
    /// <typeparam name="T">The element type.</typeparam>
    /// <param name="sample1">First sample.</param>
    /// <param name="sample2">Second sample.</param>
    /// <param name="func">Projection from element to numeric value.</param>
    /// <param name="confidenceLevel">Confidence level.</param>
    public static StatisticalTestResult TwoSampleTTest<T>(this IEnumerable<T> sample1, IEnumerable<T> sample2,
        Func<T, double> func, double confidenceLevel = 0.95)
    {
        var v1 = sample1.Select(func).ToList();
        var v2 = sample2.Select(func).ToList();

        int n1 = v1.Count, n2 = v2.Count;
        double mean1 = v1.Average(), mean2 = v2.Average();
        double s1 = v1.Sum(x => Math.Pow(x - mean1, 2)) / (n1 - 1);
        double s2 = v2.Sum(x => Math.Pow(x - mean2, 2)) / (n2 - 1);

        double t = (mean1 - mean2) / Math.Sqrt(s1 / n1 + s2 / n2);

        // Welch-Satterthwaite degrees of freedom
        double num = Math.Pow(s1 / n1 + s2 / n2, 2);
        double den = Math.Pow(s1 / n1, 2) / (n1 - 1) + Math.Pow(s2 / n2, 2) / (n2 - 1);
        double df = num / den;

        var tDist = new StudentTDistribution(df);
        double pValue = tDist.TwoTailedPValue(t);

        return new StatisticalTestResult(t, pValue, confidenceLevel, df);
    }

    /// <summary>
    /// Paired t-test: tests whether the mean difference of paired observations is zero.
    /// H₀: μ_d = 0.
    /// </summary>
    /// <typeparam name="T">The element type.</typeparam>
    /// <param name="enumerable">The input sequence of paired observations.</param>
    /// <param name="func">Projection from element to (before, after) pair.</param>
    /// <param name="confidenceLevel">Confidence level.</param>
    public static StatisticalTestResult PairedTTest<T>(this IEnumerable<T> enumerable, Func<T, (double before, double after)> func,
        double confidenceLevel = 0.95)
    {
        var diffs = enumerable.Select(func).Select(p => p.after - p.before).ToList();
        int n = diffs.Count;
        double meanD = diffs.Average();
        double sD = Math.Sqrt(diffs.Sum(d => Math.Pow(d - meanD, 2)) / (n - 1));

        double t = meanD / (sD / Math.Sqrt(n));
        double df = n - 1;

        var tDist = new StudentTDistribution(df);
        double pValue = tDist.TwoTailedPValue(t);

        return new StatisticalTestResult(t, pValue, confidenceLevel, df);
    }

    /// <summary>
    /// One-sample z-test: tests whether the population mean equals a hypothesized value
    /// when the population standard deviation is known.
    /// H₀: μ = hypothesizedMean.
    /// </summary>
    /// <typeparam name="T">The element type.</typeparam>
    /// <param name="enumerable">The input sequence.</param>
    /// <param name="func">Projection from element to numeric value.</param>
    /// <param name="hypothesizedMean">The value to test against.</param>
    /// <param name="populationStdDev">Known population standard deviation.</param>
    /// <param name="confidenceLevel">Confidence level.</param>
    public static StatisticalTestResult ZTest<T>(this IEnumerable<T> enumerable, Func<T, double> func,
        double hypothesizedMean, double populationStdDev, double confidenceLevel = 0.95)
    {
        var values = enumerable.Select(func).ToList();
        int n = values.Count;
        double mean = values.Average();

        double z = (mean - hypothesizedMean) / (populationStdDev / Math.Sqrt(n));

        var normal = new NormalDistribution();
        double pValue = 2.0 * (1.0 - normal.Cdf(Math.Abs(z)));

        return new StatisticalTestResult(z, pValue, confidenceLevel);
    }

    /// <summary>
    /// Chi-squared goodness-of-fit test: tests whether observed frequencies match expected frequencies.
    /// H₀: the data follows the expected distribution.
    /// </summary>
    /// <typeparam name="T">The element type.</typeparam>
    /// <param name="enumerable">The input sequence.</param>
    /// <param name="func">Projection from element to (observed, expected).</param>
    /// <param name="confidenceLevel">Confidence level.</param>
    public static StatisticalTestResult ChiSquaredTest<T>(this IEnumerable<T> enumerable,
        Func<T, (double observed, double expected)> func, double confidenceLevel = 0.95)
    {
        var data = enumerable.Select(func).ToList();
        double chiSq = data.Sum(p => Math.Pow(p.observed - p.expected, 2) / p.expected);
        double df = data.Count - 1;

        var chiDist = new ChiSquaredDistribution(df);
        double pValue = chiDist.UpperTailPValue(chiSq);

        return new StatisticalTestResult(chiSq, pValue, confidenceLevel, df);
    }

    /// <summary>
    /// One-way ANOVA: tests whether the means of multiple groups are equal.
    /// H₀: μ₁ = μ₂ = ... = μₖ.
    /// </summary>
    /// <typeparam name="T">The element type.</typeparam>
    /// <param name="enumerable">The input sequence.</param>
    /// <param name="func">Projection from element to (value, groupId).</param>
    /// <param name="confidenceLevel">Confidence level.</param>
    public static StatisticalTestResult OneWayAnova<T>(this IEnumerable<T> enumerable,
        Func<T, (double value, int group)> func, double confidenceLevel = 0.95)
    {
        var data = enumerable.Select(func).ToList();
        var groups = data.GroupBy(p => p.group).ToList();
        int k = groups.Count;
        int totalN = data.Count;
        double grandMean = data.Average(p => p.value);

        // Between-group sum of squares
        double ssBetween = groups.Sum(g =>
            g.Count() * Math.Pow(g.Average(p => p.value) - grandMean, 2));

        // Within-group sum of squares
        double ssWithin = groups.Sum(g =>
        {
            double groupMean = g.Average(p => p.value);
            return g.Sum(p => Math.Pow(p.value - groupMean, 2));
        });

        double dfBetween = k - 1;
        double dfWithin = totalN - k;

        double msBetween = ssBetween / dfBetween;
        double msWithin = ssWithin / dfWithin;
        double f = msBetween / msWithin;

        var fDist = new FDistribution(dfBetween, dfWithin);
        double pValue = fDist.UpperTailPValue(f);

        return new StatisticalTestResult(f, pValue, confidenceLevel, dfBetween);
    }

    // ──────────────────────────────────────────────
    //  PRIVATE HELPERS
    // ──────────────────────────────────────────────

    /// <summary>
    /// Solves an augmented matrix [A|b] via Gaussian elimination with partial pivoting.
    /// </summary>
    private static double[] GaussianElimination(double[,] matrix, int n)
    {
        for (int col = 0; col < n; col++)
        {
            // Partial pivot
            int maxRow = col;
            for (int row = col + 1; row < n; row++)
            {
                if (Math.Abs(matrix[row, col]) > Math.Abs(matrix[maxRow, col]))
                    maxRow = row;
            }
            for (int j = col; j <= n; j++)
            {
                double tmp = matrix[col, j];
                matrix[col, j] = matrix[maxRow, j];
                matrix[maxRow, j] = tmp;
            }

            // Eliminate below
            for (int row = col + 1; row < n; row++)
            {
                double factor = matrix[row, col] / matrix[col, col];
                for (int j = col; j <= n; j++)
                    matrix[row, j] -= factor * matrix[col, j];
            }
        }

        // Back substitute
        double[] result = new double[n];
        for (int i = n - 1; i >= 0; i--)
        {
            result[i] = matrix[i, n];
            for (int j = i + 1; j < n; j++)
                result[i] -= matrix[i, j] * result[j];
            result[i] /= matrix[i, i];
        }

        return result;
    }

    /// <summary>
    /// Computes fractional ranks for an array (average rank for ties).
    /// </summary>
    private static double[] ComputeRanks(double[] values)
    {
        int n = values.Length;
        var indexed = values.Select((v, i) => new { Value = v, Index = i })
                            .OrderBy(x => x.Value)
                            .ToList();

        double[] ranks = new double[n];
        int i = 0;
        while (i < n)
        {
            int j = i;
            while (j < n - 1 && indexed[j + 1].Value == indexed[i].Value)
                j++;

            double avgRank = (i + j) / 2.0 + 1.0;
            for (int k = i; k <= j; k++)
                ranks[indexed[k].Index] = avgRank;
            i = j + 1;
        }

        return ranks;
    }
}
