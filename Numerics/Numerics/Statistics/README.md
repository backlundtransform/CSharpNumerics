
## 📊 Statistics

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

Regression:

```csharp
var (slope, intercept, corr) = serie.LinearRegression(p => (p.Index, p.Value));
var expFunc = serie.ExponentialRegression(p => (p.Index, p.Value));
```


K-nearest neighbors:

```csharp
var data = new List<(double x, double y, int c)> { (7,7,0), (7,4,0), (3,4,1), (1,4,1) };
int classification = data.KnearestNeighbors(p => (p.x, p.y, p.c), (3,7), 3);
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
