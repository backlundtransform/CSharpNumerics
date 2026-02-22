
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

```csharp
var noise = new Random().GenerateNoise(4);
double median = ts.Median(p => p.Value);
double std = ts.StandardDeviation(p => p.Value);
```
Coefficient of determination:
```csharp
var data = new[] { (1.0, 5.0), (2.0, 1.0), (3.0, 4.0), (4.0, 6.0) };
double r2 = data.CoefficientOfDetermination(p => (p.Item1, p.Item2));
```


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

##  ✨ Interpolation

CSharpNumerics provides a unified interpolation API supporting linear and logarithmic scales:

* **Linear**
* **Log–Log** (log x, log y)
* **Lin–Log** (lin x, log y)
* **Log–Lin** (log x, lin y)

All methods are routed through one central function.


```csharp
public enum InterpolationType
{
    Linear,
    Logarithmic,   // log–log
    LogLin,
    LinLog
}
```

```csharp
double Interpolate<T>(
    this IEnumerable<T> source,
    Func<T, (double x, double y)> selector,
    double index,
    InterpolationType type);
```

Example:

```csharp
var data = new List<Serie>
{
    new Serie { Index = 1, Value = 10 },
    new Serie { Index = 10, Value = 100 }
};

double y = data.Interpolate(
    p => (p.Index, p.Value),
    3.5,
    InterpolationType.Linear
);
```
