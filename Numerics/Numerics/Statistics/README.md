
## � Data

The `CSharpNumerics.Statistics.Data` namespace provides core data structures for working with indexed and time-indexed data.

### Serie & TimeSerie

Simple data-point classes for single-value series:

```csharp
// Numeric index
var point = new Serie { Index = 1.0, Value = 42.0 };

// DateTime index
var tp = new TimeSerie { TimeStamp = DateTime.Now, Value = 99.5 };
```

### Series

Columnar data structure for multi-column numeric data with optional grouping.

```csharp
// Create from a list of Serie
var series = Series.FromSerie(serieList);

// Load from CSV (auto-indexes rows, supports target and group columns)
var csv = Series.FromCsv("data.csv", sep: ',', targetColumn: "Price", groupColumn: "Region");

// Access
int[] index     = csv.Index;
double[][] data = csv.Data;   // data[col][row]
string[] cols   = csv.Cols;
int[] groups    = csv.Groups;

// Convert to Matrix (optionally exclude a column)
Matrix m = csv.ToMatrix(excludeCol: 0);
```

### TimeSeries

Columnar data structure indexed by `DateTime`, designed for time-based analysis.

```csharp
// Create from a list of TimeSerie
var ts = TimeSeries.FromTimeSerie(timeSerieList, columnName: "Temperature");

// Load from CSV (first column = DateTime, remaining = features)
var csv = TimeSeries.FromCsv("timeseries.csv");

// Access
DateTime[] time = csv.Time;
double[][] data = csv.Data;   // data[col][row]
string[] cols   = csv.Cols;
int rows        = csv.RowCount;
int columns     = csv.ColumnCount;

// Convert to Matrix
Matrix m = csv.ToMatrix();
```

### DataExtensions

Utility methods for generating and resampling series data.

```csharp
// Generate Serie from a function
Func<double, double> f = x => Math.Sin(x);
List<Serie> serie = f.GetSeries(minValue: 0, maxValue: 10, stepSize: 100);

// Generate TimeSerie from a function
Func<DateTime, double> g = t => t.Hour * 1.5;
List<TimeSerie> ts = g.GetTimeSeries(startTime, endTime, stepSize: 50);

// Resample a TimeSerie to equidistant steps (interpolation)
List<TimeSerie> resampled = timeSeries.GenerateTimeSerieWithEquivalentSteps(
    minutes: 15, startDate: start, endDate: end);

// Resample with grouping operator (Average, Max, Min, Sum, Median)
List<TimeSerie> grouped = timeSeries.GenerateTimeSerieWithEquivalentSteps(
    minutes: 60, startDate: start, endDate: end,
    groupOperator: GroupOperator.Average, multiplier: 1, shouldInterpolate: true);
```

---

## �📊 Statistics

### Descriptive Statistics

All methods follow the generic selector pattern `IEnumerable<T>.Method(Func<T, double>)`.

```csharp
var data = new[] { 3.0, 7.0, 7.0, 2.0, 5.0, 1.0, 9.0, 4.0 };

double median   = data.Median(x => x);           // 4.5
double variance = data.Variance(x => x);          // population variance (÷n)
double sVar     = data.SampleVariance(x => x);    // sample variance (÷(n-1))
double std      = data.StandardDeviation(x => x); // population std dev
double mode     = data.Mode(x => x);              // 7.0
double range    = data.Range(x => x);             // 8.0
double p90      = data.Percentile(x => x, 90);    // 90th percentile
double iqr      = data.InterquartileRange(x => x);// Q3 - Q1
double skew     = data.Skewness(x => x);          // sample skewness
double kurt     = data.Kurtosis(x => x);          // excess kurtosis (normal = 0)
```

Bivariate selectors use `Func<T, (double x, double y)>`:
```csharp
var pairs = new[] { (1.0, 5.0), (2.0, 1.0), (3.0, 4.0), (4.0, 6.0) };
double cov = pairs.Covariance(p => (p.Item1, p.Item2));
double r2  = pairs.CoefficientOfDetermination(p => (p.Item1, p.Item2));
```

Cumulative sum and confidence intervals:
```csharp
double[] cumSum = data.CumulativeSum(x => x).ToArray();
var (lo, hi) = data.ConfidenceIntervals(x => x, 0.95);
```

### Inferential Statistics — Regression

```csharp
var points = new[] { (1.0, 2.1), (2.0, 3.9), (3.0, 6.2), (4.0, 7.8) };

// Linear regression: y = slope * x + intercept
var (slope, intercept, r) = points.LinearRegression(p => (p.Item1, p.Item2));

// Exponential: y = a * e^(bx)
Func<double, double> expFit = points.ExponentialRegression(p => (p.Item1, p.Item2));

// Logarithmic: y = a + b * ln(x)
Func<double, double> logFit = points.LogarithmicRegression(p => (p.Item1, p.Item2));

// Power: y = a * x^b
Func<double, double> powFit = points.PowerRegression(p => (p.Item1, p.Item2));

// Polynomial (degree N): y = a₀ + a₁x + a₂x² + ...
var (predict, coefficients) = points.PolynomialRegression(p => (p.Item1, p.Item2), degree: 2);
double yHat = predict(2.5);
```

### Inferential Statistics — Correlation

```csharp
// Pearson correlation coefficient with p-value
var (r, pValue) = data.PearsonCorrelation(p => (p.X, p.Y));

// Spearman rank correlation with p-value
var (rho, pVal) = data.SpearmanCorrelation(p => (p.X, p.Y));
```

### Inferential Statistics — Hypothesis Tests

All tests return a `StatisticalTestResult` with `TestStatistic`, `PValue`, `RejectNull`, `ConfidenceLevel`, and `DegreesOfFreedom`.

```csharp
// One-sample t-test: is the mean equal to 50?
StatisticalTestResult t1 = scores.OneSampleTTest(s => s.Value, hypothesizedMean: 50);

// Two-sample t-test (Welch's): are two group means different?
StatisticalTestResult t2 = groupA.TwoSampleTTest(groupB, s => s.Value);

// Paired t-test: before/after comparison
StatisticalTestResult t3 = patients.PairedTTest(p => (p.Before, p.After));

// Z-test: known population standard deviation
StatisticalTestResult z = measurements.ZTest(m => m.Value, hypothesizedMean: 100, populationStdDev: 15);

// Chi-squared goodness-of-fit
var bins = new[] { (observed: 50.0, expected: 40.0), (30.0, 40.0), (20.0, 20.0) };
StatisticalTestResult chi = bins.ChiSquaredTest(b => (b.observed, b.expected));

// One-way ANOVA: compare means across groups
StatisticalTestResult anova = allSamples.OneWayAnova(s => (s.Value, s.GroupId));

// Inspect the result
Console.WriteLine(anova.TestStatistic); // F-value
Console.WriteLine(anova.PValue);        // p-value
Console.WriteLine(anova.RejectNull);    // true/false at 95% confidence
```

### Hypothesis Tests — Fluent API

Cleaner extensions directly on `IEnumerable<double>` with `Alternative` enum (TwoSided, Less, Greater),
confidence intervals, and effect sizes. All return `HypothesisTestResult`.

```csharp
var sample = new[] { 12.1, 11.3, 10.8, 11.9, 12.5, 11.0 };

// One-sample t-test: H₀: μ = 10
var result = sample.TTest(mu: 10.0, alternative: Alternative.TwoSided);
result.PValue                // p-value
result.RejectNull            // true
result.ConfidenceIntervalLower / result.ConfidenceIntervalUpper
result.EffectSize            // Cohen's d

// One-sided: is the mean greater than 10?
var greater = sample.TTest(mu: 10.0, alternative: Alternative.Greater);
```

```csharp
// Two-sample Welch's t-test
var groupA = new[] { 10.1, 10.3, 9.8, 10.5 };
var groupB = new[] { 5.1, 4.9, 5.3, 5.0 };
var t2 = groupA.TTest(groupB, Alternative.TwoSided);

// Paired t-test
var before = new[] { 60.0, 62.0, 58.0, 65.0 };
var after  = new[] { 70.0, 73.0, 68.0, 75.0 };
var paired = before.PairedTTest(after, Alternative.Greater);
```

```csharp
// Z-test (known σ)
var z = sample.ZTest(mu: 10.0, populationStdDev: 1.5, alternative: Alternative.TwoSided);

// Chi-squared goodness-of-fit
var observed = new[] { 50.0, 30.0, 20.0 };
var expected = new[] { 40.0, 30.0, 30.0 };
var chi = observed.ChiSquaredTest(expected);

// F-test for equality of variances
var fTest = groupA.FTest(groupB, Alternative.TwoSided);
```

```csharp
// Non-parametric tests
var u = groupA.MannWhitneyUTest(groupB, Alternative.TwoSided);   // Mann-Whitney U
var w = diffs.WilcoxonSignedRankTest(Alternative.TwoSided);        // Wilcoxon signed-rank

// One-way ANOVA
var groups = new List<IEnumerable<double>> { groupA, groupB, groupC };
var anova = groups.Anova(confidenceLevel: 0.95);
anova.EffectSize   // η²
```


## 📐 Fitting

The `CSharpNumerics.Statistics.Fitting` namespace provides a comprehensive curve-fitting toolkit: ordinary, weighted, nonlinear, and robust least squares, plus goodness-of-fit metrics, residual analysis, and parameter estimation.

### Least Squares (OLS)

Polynomial and multiple-regression fitting via the normal equations.

```csharp
using CSharpNumerics.Statistics.Fitting;
using CSharpNumerics.Numerics.Objects;

// Linear fit: y = a₀ + a₁x
var x = new VectorN(new[] { 1.0, 2, 3, 4, 5 });
var y = new VectorN(new[] { 2.1, 3.9, 6.2, 7.8, 10.1 });
FittingResult lr = LeastSquaresFitter.Fit(x, y, degree: 1);
// lr.Coefficients[0] → intercept, lr.Coefficients[1] → slope

// Polynomial fit: y = a₀ + a₁x + a₂x²
FittingResult poly = LeastSquaresFitter.Fit(x, y, degree: 2);

// Multiple regression with a custom design matrix (include intercept column)
double[,] X = { { 1, 2.0, 3.5 }, { 1, 4.0, 1.2 }, { 1, 6.0, 5.8 }, /* … */ };
var response = new VectorN(new[] { 10.0, 20, 30 });
FittingResult mr = LeastSquaresFitter.Fit(X, response);
```

### Weighted Least Squares (WLS)

Down-weight unreliable observations.

```csharp
var weights = new VectorN(new[] { 1.0, 1, 0.01, 1, 1 }); // low weight on outlier
FittingResult wls = WeightedLeastSquaresFitter.Fit(x, y, weights, degree: 1);

// With a design matrix
FittingResult wlsM = WeightedLeastSquaresFitter.Fit(designMatrix, y, weights);
```

### Nonlinear Least Squares (Levenberg-Marquardt)

Fit any user-defined model with numerical Jacobian.

```csharp
// Exponential model: y = a · exp(b · x)
Func<double, VectorN, double> model = (xi, p) => p[0] * Math.Exp(p[1] * xi);
var init = new VectorN(new[] { 1.0, 0.1 });

FittingResult nlr = NonlinearLeastSquaresFitter.Fit(model, x, y, init);
// nlr.Coefficients[0] → a, nlr.Coefficients[1] → b

// Logistic model
Func<double, VectorN, double> logistic = (xi, p) =>
    p[0] / (1.0 + Math.Exp(-p[1] * (xi - p[2])));
FittingResult logFit = NonlinearLeastSquaresFitter.Fit(logistic, x, y,
    initialParameters: new VectorN(new[] { 10.0, 1.0, 5.0 }),
    maxIterations: 500, tolerance: 1e-12);
```

### Robust Fitting (IRLS)

Resist outliers using Iteratively Reweighted Least Squares with configurable weight functions.

```csharp
// Huber weights (default)
FittingResult robust = RobustFitter.Fit(x, y, degree: 1);

// Tukey bisquare weights
FittingResult tukey = RobustFitter.Fit(x, y, degree: 1,
    weightFunction: RobustWeightFunction.TukeyBisquare);

// Available weight functions: Huber, TukeyBisquare, Cauchy, Andrews
// Custom tuning constant
FittingResult custom = RobustFitter.Fit(x, y, degree: 2,
    weightFunction: RobustWeightFunction.Huber, tuningConstant: 2.0);
```

### Goodness of Fit & Residual Analysis

```csharp
var observed  = new VectorN(new[] { 1.0, 2, 3, 4, 5 });
var predicted = new VectorN(new[] { 1.1, 2.0, 2.9, 4.1, 4.9 });

double r2    = GoodnessOfFit.RSquared(observed, predicted);
double adjR2 = GoodnessOfFit.AdjustedRSquared(r2, n: 5, parameterCount: 2);
double rmse  = GoodnessOfFit.RMSE(observed, predicted);
double mae   = GoodnessOfFit.MAE(observed, predicted);
double sse   = GoodnessOfFit.SSE(observed, predicted);
double sst   = GoodnessOfFit.SST(observed);
double mse   = GoodnessOfFit.MSE(observed, predicted);
double aic   = GoodnessOfFit.AIC(observed, predicted, parameterCount: 2);
double bic   = GoodnessOfFit.BIC(observed, predicted, parameterCount: 2);

// Full residual diagnostics
ResidualAnalysis ra = GoodnessOfFit.AnalyzeResiduals(observed, predicted);
ra.StandardizedResiduals   // residuals / σ
ra.DurbinWatsonStatistic   // ≈ 2 → no autocorrelation
ra.MeanResidual            // should be near 0
ra.Skewness                // 0 = symmetric
ra.Kurtosis                // 0 = normal tails
```

### Parameter Estimation

Confidence/prediction intervals, MLE, and method-of-moments estimators.

```csharp
FittingResult fit = LeastSquaresFitter.Fit(x, y);

// 95 % confidence intervals for coefficients (returns VectorN lower/upper)
var (lower, upper) = ParameterEstimation.ConfidenceIntervals(fit, confidenceLevel: 0.95);

// Prediction interval at a new point
double[,] designMatrix = { /* the X used for fitting */ };
var xNew = new VectorN(new[] { 1.0, 3.0 }); // [intercept, x=3]
var (lo, hi) = ParameterEstimation.PredictionInterval(fit, designMatrix, xNew);

// Maximum Likelihood Estimation (numerical optimisation)
Func<VectorN, double> negLogL = theta => { /* −ln L(θ|data) */ };
VectorN mle = ParameterEstimation.MaximumLikelihoodEstimate(negLogL, initialParameters);

// Method of Moments (all accept VectorN)
var (mean, var)      = ParameterEstimation.MethodOfMomentsNormal(data);
double rate          = ParameterEstimation.MethodOfMomentsExponential(data);
double lambda        = ParameterEstimation.MethodOfMomentsPoisson(data);
var (shape, scale)   = ParameterEstimation.MethodOfMomentsGamma(data);
```

### FittingResult

Every fitter returns a `FittingResult` with:

| Property | Description |
|---|---|
| `Coefficients` | `VectorN` — estimated parameters [a₀, a₁, …] |
| `StandardErrors` | `VectorN` — SE for each coefficient |
| `Residuals` | `VectorN` — observed − fitted |
| `FittedValues` | `VectorN` — model predictions |
| `RSquared` | R² |
| `AdjustedRSquared` | Penalised R² |
| `RMSE` | Root mean squared error |
| `MAE` | Mean absolute error |
| `DegreesOfFreedom` | n − p |
| `ObservationCount` | n |
| `ParameterCount` | p |


## 📡 Time Series Analysis

The `CSharpNumerics.Statistics.TimeSeriesAnalysis` namespace provides specialised algorithms for periodic signal detection, detrending, phase folding, and peak/event fitting in time-series data — including support for **unevenly sampled** observations.

### Lomb–Scargle Periodogram

Spectral analysis for unevenly sampled data. The classical Lomb (1976) / Scargle (1982) algorithm with the Townsend τ-offset for numerical stability.

```csharp
using CSharpNumerics.Statistics.TimeSeriesAnalysis;
using CSharpNumerics.Numerics.Objects;

// Observations at irregular times
var times  = new VectorN(new[] { 0.0, 1.1, 2.3, 3.0, 4.7, 5.2, /* … */ });
var values = new VectorN(new[] { 0.1, 0.9, 0.0, -0.8, 0.2, 0.95, /* … */ });

// Frequency-range scan
PeriodogramResult ls = LombScarglePeriodogram.Compute(
    times, values,
    minFrequency: 0.01, maxFrequency: 1.0, numFrequencies: 1000);

double bestFreq   = ls.BestFrequency;   // strongest spectral peak
double bestPeriod  = ls.BestPeriod;      // 1 / bestFreq
double bestPower   = ls.BestPower;

// Period-range scan (linearly spaced in period)
PeriodogramResult lsP = LombScarglePeriodogram.ComputeByPeriod(
    times, values, minPeriod: 1, maxPeriod: 50, numPeriods: 500);

// False-alarm probability
double fap = LombScarglePeriodogram.FalseAlarmProbability(
    bestPower, numObservations: times.Length, numFrequencies: 1000);
```

### Box-fitting Least Squares (BLS)

Transit detection via the Kovács, Zucker & Mazeh (2002) algorithm. Phase-folds the light curve at each trial period and searches for the optimal box-shaped dip.

```csharp
// Search for transits between 2 and 20 days
BLSResult bls = BoxFittingLeastSquares.Compute(
    times, flux,
    minPeriod: 2, maxPeriod: 20, numPeriods: 500,
    minDurationFraction: 0.01, maxDurationFraction: 0.15, numBins: 200);

double period   = bls.BestPeriod;
double epoch    = bls.BestEpoch;          // mid-transit time
double depth    = bls.BestDepth;          // flux decrement
double durFrac  = bls.BestDurationFraction;

// Or supply explicit trial periods
var trialPeriods = new VectorN(new[] { 3.0, 5.0, 7.0, 10.0 });
BLSResult bls2 = BoxFittingLeastSquares.Compute(times, flux, trialPeriods);
```

### Phase Folding

Fold a time series onto a single cycle of a trial period — essential for visualising periodic signals.

```csharp
// Raw fold (sorted by phase)
var (phases, vals) = PhaseFolding.Fold(times, values, period: 12.5, epoch: 0);

// Binned fold (average per bin)
var (binCenters, binMeans) = PhaseFolding.FoldAndBin(
    times, values, period: 12.5, numBins: 20);

// Binned fold with per-bin standard deviations
var (centers, means, stds) = PhaseFolding.FoldAndBinWithErrors(
    times, values, period: 12.5, numBins: 20);
```

### Trend Fitting / Detrending

Remove systematic trends before period analysis. All methods return both the detrended signal and the estimated trend.

```csharp
// Polynomial detrending (uses LeastSquaresFitter internally)
var (detrended, trend) = TimeSeriesDetrending.PolynomialDetrend(times, values, degree: 2);

// Moving-average (centred window, handles edges)
var (det, ma) = TimeSeriesDetrending.MovingAverageDetrend(values, windowSize: 15);

// Median filter (robust to outliers)
var (det2, med) = TimeSeriesDetrending.MedianFilterDetrend(values, windowSize: 15);

// Savitzky–Golay smoothing (local polynomial least squares)
var (det3, sg) = TimeSeriesDetrending.SavitzkyGolayDetrend(values, windowSize: 11, polyOrder: 3);

// Differencing (output length = n − order)
VectorN diff1 = TimeSeriesDetrending.Difference(values, order: 1);
VectorN diff2 = TimeSeriesDetrending.Difference(values, order: 2);
```

### Peak / Event Fitting

Detect and fit peaks (or valleys) in a time series with Gaussian or Lorentzian profiles via Levenberg-Marquardt.

```csharp
// Find peaks with prominence and distance filtering
List<int> peaks = PeakFitting.FindPeaks(values, prominence: 0.5, minDistance: 5);
List<int> dips  = PeakFitting.FindValleys(values, prominence: 0.3);

// Fit a Gaussian to a region around a known peak
PeakResult g = PeakFitting.FitGaussian(xRegion, yRegion,
    centerEstimate: 10.0, widthEstimate: 2.0);
// g.Center, g.Amplitude, g.Width, g.Baseline, g.RSquared

// Fit a Lorentzian
PeakResult lor = PeakFitting.FitLorentzian(xRegion, yRegion,
    centerEstimate: 10.0, widthEstimate: 2.0);

// Find and fit all peaks automatically
List<PeakResult> fitted = PeakFitting.FindAndFitPeaks(
    x, y, prominence: 0.5, shape: PeakShape.Gaussian, fitRadius: 15);
```

| Property | Description |
|----------|-------------|
| `Center` | Fitted peak centre (x-axis) |
| `Amplitude` | Fitted height above baseline |
| `Width` | σ (Gaussian) or γ (Lorentzian) |
| `Baseline` | Constant offset |
| `Shape` | `PeakShape.Gaussian` or `PeakShape.Lorentzian` |
| `RSquared` | Goodness-of-fit for the peak region |
| `PeakIndex` | Index in the original series |

### False Alarm Probability

Estimate the statistical significance of periodogram peaks using analytical (Baluev 2008) or bootstrap methods.

```csharp
using CSharpNumerics.Statistics.TimeSeriesAnalysis;

// Analytical FAP (fast) — Baluev 2008 approximation
double fap = FalseAlarmProbability.AnalyticalFAP(
    peakPower: 12.5, numFrequencies: 1000, numDataPoints: 500);
// fap ≈ 0.0001 → highly significant detection

// Bootstrap FAP (accurate, slower) — shuffles values and re-evaluates periodogram
double fapBoot = FalseAlarmProbability.BootstrapFAP(
    times, values, peakPower: 12.5,
    nBootstrap: 1000, minPeriod: 1.0, maxPeriod: 50.0,
    numFrequencies: 500, seed: 42);
```

### Sliding Window Statistics

Running (sliding window) descriptive statistics with centered windows. Edge windows are computed with fewer points. Even window sizes are incremented to the next odd value for symmetric centering.

```csharp
using CSharpNumerics.Statistics.TimeSeriesAnalysis;

var data = new[] { 1.0, 3.0, 5.0, 7.0, 2.0, 4.0, 6.0 };

double[] mean   = SlidingWindowStatistics.RunningMean(data, windowSize: 3);
double[] std    = SlidingWindowStatistics.RunningStd(data, windowSize: 3);
double[] mad    = SlidingWindowStatistics.RunningMAD(data, windowSize: 5);
double[] median = SlidingWindowStatistics.RunningMedian(data, windowSize: 3);

// All outputs have the same length as the input array
// MAD × 1.4826 ≈ standard deviation for normally distributed data
```

---

## 🛡️ Robust Statistics

The `CSharpNumerics.Statistics.Robust` namespace provides outlier-resistant statistical methods.

### Sigma Clipping

Iterative sigma-clipping: removes data points outside [mean − σ_low × σ, mean + σ_high × σ] and recomputes until convergence.

```csharp
using CSharpNumerics.Statistics.Robust;

var data = new[] { 1.0, 2.0, 2.5, 2.3, 100.0, 2.1, 1.9, -50.0, 2.4 };

// Full result with mask, mean, std, iteration count
ClipResult result = SigmaClipping.Clip(data, sigmaLow: 3.0, sigmaHigh: 3.0, maxIter: 10);
result.Mean           // mean of retained points
result.Std            // std of retained points
result.RetainedCount  // number of non-clipped points
result.ClippedCount   // number of removed outliers
result.Iterations     // convergence iterations
result.Mask           // bool[] — true = retained

// Symmetric shorthand
ClipResult sym = SigmaClipping.Clip(data, sigma: 2.5, maxIter: 5);

// Just get the cleaned array
double[] clean = SigmaClipping.Apply(data, sigmaLow: 3.0, sigmaHigh: 3.0);
// clean ≈ { 1.0, 2.0, 2.5, 2.3, 2.1, 1.9, 2.4 }
```

---


## 🎲 Random

The `CSharpNumerics.Statistics.Random` namespace provides `RandomGenerator` — a seedable random number engine with advanced sampling methods built on top of `System.Random`.

```csharp
using CSharpNumerics.Statistics.Random;

// Seedable for reproducible results
var rng = new RandomGenerator(seed: 42);

// Uniform
double u = rng.NextUniform(2.0, 5.0);

// Gaussian (Box-Muller transform)
double g = rng.NextGaussian(mean: 0, standardDeviation: 1);

// Exponential (inverse transform sampling)
double e = rng.NextExponential(lambda: 2.0);

// Poisson (Knuth / rejection)
int p = rng.NextPoisson(lambda: 4.0);

// Bernoulli & Binomial
int coin = rng.NextBernoulli(0.5);
int hits = rng.NextBinomial(n: 20, p: 0.3);

// Batch sampling
double[] gaussianSamples = rng.GaussianSamples(1000);
double[] uniformSamples   = rng.UniformSamples(1000, min: 0, max: 10);

// Shuffle & sample without replacement
var deck = new[] { 1, 2, 3, 4, 5 };
rng.Shuffle(deck);
var hand = rng.Sample(deck, k: 3);
```

---

## 📈 Distributions

The `CSharpNumerics.Statistics.Distributions` namespace provides probability distributions with a common `IDistribution` interface:

| Distribution | Class | Parameters |
|---|---|---|
| Uniform | `UniformDistribution` | a, b |
| Normal | `NormalDistribution` | μ, σ |
| Exponential | `ExponentialDistribution` | λ |
| Poisson | `PoissonDistribution` | λ |
| Bernoulli | `BernoulliDistribution` | p |
| Binomial | `BinomialDistribution` | n, p |
| Student's t | `StudentTDistribution` | ν (degrees of freedom) |
| Chi-squared | `ChiSquaredDistribution` | k (degrees of freedom) |
| F | `FDistribution` | d1, d2 (degrees of freedom) |

Every distribution exposes `Mean`, `Variance`, `StandardDeviation`, `Pdf(x)`, `Cdf(x)`, `Sample(rng)` and `Samples(rng, count)`.

```csharp
using CSharpNumerics.Statistics.Distributions;
using CSharpNumerics.Statistics.Random;

var normal = new NormalDistribution(mu: 100, sigma: 15);

double pdf = normal.Pdf(100);      // peak value
double cdf = normal.Cdf(115);      // P(X ≤ 115) ≈ 0.8413
double q   = normal.InverseCdf(0.95); // z such that P(X ≤ z) = 0.95

var rng     = new RandomGenerator(42);
double[] samples = normal.Samples(rng, 10_000);
```

```csharp
var poisson = new PoissonDistribution(lambda: 3.0);
double pmf  = poisson.Pdf(2);     // P(X = 2)
double cumul = poisson.Cdf(4);    // P(X ≤ 4)

var binomial = new BinomialDistribution(n: 10, p: 0.3);
// PMF sums to 1
double total = Enumerable.Range(0, 11).Sum(k => binomial.Pdf(k));
```

Student's t, chi-squared and F distributions are used by the hypothesis testing methods, but can also be used directly:

```csharp
var t = new StudentTDistribution(degreesOfFreedom: 29);
double pTwoTail = t.TwoTailedPValue(2.045);   // two-tailed p-value
double quantile  = t.InverseCdf(0.975);         // critical value

var chi2 = new ChiSquaredDistribution(degreesOfFreedom: 5);
double pUpper = chi2.UpperTailPValue(11.07);   // P(X ≥ 11.07)

var f = new FDistribution(d1: 3, d2: 20);
double pF = f.UpperTailPValue(3.10);           // P(F ≥ 3.10)
```

---

## 🎯 Monte Carlo

The `CSharpNumerics.Statistics.MonteCarlo` namespace provides a general-purpose Monte Carlo simulation engine.

### Quick start — estimate π

```csharp
using CSharpNumerics.Statistics.MonteCarlo;

var sim = new MonteCarloSimulator(seed: 42);
var result = sim.Run(rng =>
{
    double x = rng.NextUniform(0, 1);
    double y = rng.NextUniform(0, 1);
    return x * x + y * y <= 1.0 ? 1.0 : 0.0;
}, iterations: 1_000_000);

double piEstimate = result.Mean * 4.0; // ≈ 3.1416
```

### Monte Carlo integration

```csharp
// ∫₀^π sin(x) dx = 2
var result = sim.IntegrateUniform(Math.Sin, a: 0, b: Math.PI, iterations: 500_000);
// result.Mean ≈ 2.0

// E[X²] where X ~ N(0,1) → should be ≈ 1.0
var normal = new NormalDistribution(0, 1);
var result2 = sim.Integrate(x => x * x, normal, iterations: 200_000);
```

### Custom models via IMonteCarloModel

```csharp
public class TwoDiceModel : IMonteCarloModel
{
    public double Evaluate(RandomGenerator rng)
    {
        return rng.NextInt(1, 7) + rng.NextInt(1, 7);
    }
}

var result = sim.Run(new TwoDiceModel(), iterations: 200_000);
// result.Mean ≈ 7.0
```

### Result analysis

`MonteCarloResult` provides descriptive statistics, percentiles, histograms and probability queries:

```csharp
result.Mean                         // Sample mean
result.Variance                     // Sample variance (n-1)
result.StandardDeviation            // √Variance
result.StandardError                // σ / √n (precision of estimate)
result.Median                       // 50th percentile
result.Min / result.Max

result.Percentile(95);              // 95th percentile
result.ConfidenceInterval(0.95);    // (lower, upper) from empirical percentiles
result.Probability(x => x > 10);   // Fraction of samples satisfying predicate

var histogram = result.Histogram(bins: 20);  // (binCenter, count)[]
```

### Multivariate & convergence

```csharp
// Multiple outputs per trial
var results = sim.RunMultivariate(rng => new[]
{
    rng.NextGaussian(0, 1),
    rng.NextGaussian(10, 2)
}, iterations: 100_000);
// results[0].Mean ≈ 0, results[1].Mean ≈ 10

// Track running mean to verify convergence
double[] curve = sim.ConvergenceCurve(
    rng => rng.NextGaussian(5.0, 1.0), iterations: 50_000);
```

---

