
CSharpNumerics includes a **lightweight, fully numerical machine learning framework** designed for **research, experimentation, and educational use**.
The focus is on **transparency**, **mathematical clarity**, and **pipeline-based model evaluation** — not black-box automation.
## ➡️ Pipeline Grid
All models are implemented directly on top of the library’s **Matrix** and **Vector** primitives.

Models can be combined with:

* **Scalers** (e.g. `StandardScaler`)
* **Feature selectors** (e.g. `SelectKBest`)
* **Dimensionality reducers** (e.g. `PCA`)
* **Cross-validation strategies**
* **Hyperparameter search grids**


```csharp
var pipelineGrid = new PipelineGrid()
    .AddModel<RandomForest>(g => g
        .Add("NumTrees", 50, 100, 200)
        .Add("MaxDepth", 5, 8, 10))
    .AddModel<Logistic>(g => g
        .Add("LearningRate", 0.05, 0.1)
        .Add("MaxIterations", 1000, 2000)
        .AddScaler<StandardScaler>(s => { })
        .AddSelector<SelectKBest>(s => s
           .Add("K", 1, 2)))
    .AddModel<DecisionTree>(g => g
        .Add("MaxDepth", 3, 5, 8))
    .AddModel<KNearestNeighbors>(g => g
        .Add("K", 3, 5, 7));
```

CSharpNumerics supports multiple cross-validation strategies for **time series** and **tabular data**:

---

## ⏩ Rolling Cross-Validation

Train on first folds, validate on the next fold, then roll forward. Works for **classification** and **regression**.

**Example visualization**

```
Train: [1 2 3] | Test: [4]
Train: [1 2 3 4] | Test: [5]
Train: [1 2 3 4 5] | Test: [6]
...
```

```csharp
var cv = new RollingCrossValidator(pipelineGrid);
var result = cv.Run(X, y); 
var bestModel = result.BestPipeline; 
var score = result.BestScore;
```

**Key points:**

* Always respects **temporal order**
* Prevents **data leakage**
* Works well for **time series forecasting**

---

## 🔁 K-Fold Cross-Validation

Split data into **K equally sized folds**. Each fold is used once as test while remaining folds form the training set. Works for **classification** and **regression** on tabular data.

**Visualization (K = 5)**

```
Data: [ 1 2 3 4 5 ]

Fold 1: Train [2 3 4 5] | Test [1]
Fold 2: Train [1 3 4 5] | Test [2]
Fold 3: Train [1 2 4 5] | Test [3]
Fold 4: Train [1 2 3 5] | Test [4]
Fold 5: Train [1 2 3 4] | Test [5]
```

```csharp
var cv = new KFoldCrossValidator(pipelineGrid, folds: 5);
var result = cv.Run(X, y);

var bestModel = result.BestPipeline;
var score = result.BestScore;
```

**Key points:**

* **Order of samples does not matter**
* No temporal assumptions
* All samples are evaluated **exactly once**

---

## 🧮 Stratified K-Fold Cross-Validation

Used for **classification** with **imbalanced classes**. Ensures that each fold preserves the **class proportions**.

**Example visualization (K = 5)**

```
Class distribution in dataset: 90% class 0, 10% class 1

Fold 1: Train -> 80% class0 / 20% class1 | Test -> 90% class0 / 10% class1
Fold 2: Train -> 80% class0 / 20% class1 | Test -> 90% class0 / 10% class1
...
```

```csharp
var cv = new StratifiedKFoldCrossValidator(pipelineGrid, folds: 5);
var result = cv.Run(X, y); // y contains class labels

var bestModel = result.BestPipeline;
var score = result.BestScore;
```

**Key points:**

* Maintains **class distribution** in every fold
* Works **only for classification**
* Ideal for **imbalanced datasets**

---

## 🔀 ShuffleSplit Cross-Validation

Randomly splits data into a **training set** and a **test set** multiple times. Works for **classification** and **regression**. Unlike K-Fold, not all samples are guaranteed to appear in a test set.

**Example visualization (3 splits, 20% test size)**

```
Split 1: Train [1 2 3 4] | Test [5]
Split 2: Train [1 3 4 5] | Test [2]
Split 3: Train [2 3 4 5] | Test [1]
...
```

```csharp
var cv = new ShuffleSplitCrossValidator(
    pipelineGrid,
    n_splits: 5,
    testSize: 0.2,
    trainSize: 0.8,
    randomState: 42);

var result = cv.Run(X, y);

var bestModel = result.BestPipeline;
var score = result.BestScore;
```

**Key points:**

* **Randomly shuffles** data before each split
* Can perform **multiple iterations** (`n_splits`)
* **Does not guarantee** all samples are tested exactly once
* Useful for **large datasets** where full K-Fold is costly
* Can be combined with **Pipelines**, **Series**, or **TimeSeries**

---

## 🎲 Monte Carlo Cross-Validation

Runs **many random train/test splits** (typically 100–1000) and collects the resulting score into a full **probability distribution**. Built on top of the library's `MonteCarloSimulator` engine.

Unlike ShuffleSplit which returns a single aggregate score, Monte Carlo CV returns a complete `MonteCarloResult` with **confidence intervals**, **histograms**, **standard error**, and a **convergence curve**.

**Example visualization (200 iterations, 20% test size)**

```
Iteration   1: Train [random 80%] | Test [random 20%] → score = 0.88
Iteration   2: Train [random 80%] | Test [random 20%] → score = 0.85
Iteration   3: Train [random 80%] | Test [random 20%] → score = 0.91
...
Iteration 200: Train [random 80%] | Test [random 20%] → score = 0.87

→ Mean = 0.87, StdDev = 0.03, 95% CI = [0.84, 0.90]
```

**Standard usage (drop-in `ICrossValidator`)**

```csharp
var cv = new MonteCarloCrossValidator(
    pipelineGrid,
    iterations: 200,
    testSize: 0.2,
    seed: 42);

var result = cv.Run(X, y);

var bestModel = result.BestPipeline;
var score = result.BestScore;
```

**Extended usage (full score distributions)**

```csharp
var cv = new MonteCarloCrossValidator(
    pipelineGrid,
    iterations: 200,
    testSize: 0.2,
    seed: 42);

var detailed = cv.RunDetailed(X, y);

// Confidence interval for the best pipeline
var ci = detailed.BestConfidenceInterval;       // e.g. (0.84, 0.90)
double stdDev = detailed.BestScoreStdDev;        // e.g. 0.03

// Convergence curve — verify that enough iterations were run
double[] convergence = detailed.ConvergenceCurve;

// Full MonteCarloResult per pipeline
foreach (var (pipeline, mcResult) in detailed.DetailedScores)
{
    Console.WriteLine($"{pipeline} → {mcResult.Mean:F3} ± {mcResult.StandardDeviation:F3}");
    Console.WriteLine($"  95% CI: [{mcResult.ConfidenceInterval(0.95).lower:F3}, {mcResult.ConfidenceInterval(0.95).upper:F3}]");
    Console.WriteLine($"  SE: {mcResult.StandardError:F4}");
    
    // Histogram of score distribution
    var histogram = mcResult.Histogram(10);
}
```

**Key points:**

* Quantifies **model evaluation uncertainty** — not just a point estimate
* Reports **confidence intervals** for scores (e.g. "accuracy = 0.87 ± 0.03")
* **Convergence curve** shows whether enough iterations were run
* **Histogram** visualizes the full score distribution
* **Standard error** decreases with more iterations ($SE = \sigma / \sqrt{n}$)
* All pipelines are evaluated on **identical random splits** (fair comparison)
* Implements `ICrossValidator` — **drop-in replacement** for other validators
* Reproducible via `seed` parameter

---
## 📅 Leave-One-Out Cross-Validation

Train on all rows except one, test on the held-out row, then iterate. Works for **tabular or grouped data**.

**Example visualization**

```
Data: [ 1 2 3 4 5 ]

Fold 1: Train [2 3 4 5] | Test [1]
Fold 2: Train [1 3 4 5] | Test [2]
Fold 3: Train [1 2 4 5] | Test [3]
Fold 4: Train [1 2 3 5] | Test [4]
Fold 5: Train [1 2 3 4] | Test [5]
```

```csharp
var cv = new LeaveOneOutCrossValidator(pipelineGrid);
var result = cv.Run(X, y);

var bestModel = result.BestPipeline;
var score = result.BestScore;
```

**Key points:**

* Extreme case of K-Fold where **K = n**
* Guarantees each sample is used as test **exactly once**
* Can be combined with **groups** if needed

---

## 📦 Grouped Cross-Validation

Used when **samples belong to groups** and all samples from the same group must stay together. Works for **classification** and **regression**.

**Example visualization 📊 Series**

```
Groups: [A] [B] [C] [D] [E]

Fold 1: Train -> B, C, D, E | Test -> A
Fold 2: Train -> A, C, D, E | Test -> B
Fold 3: Train -> A, B, D, E | Test -> C
...
```

```csharp
var cv = new LeaveOneOutCrossValidator(pipelineGrid);
var result = cv.Run(series, targetColumn: "Target", groupColumn: "Department"); 
```

**Key points:**

* Groups can be anything: **customer, company, department, gender**
* Ensures **all group members stay together**
* Often called **Leave-One-Group-Out**

**Example visualization ⏱️ TimeSeries**

Train on all groups except one, test on the held-out group, then iterate. Groups can be days, weeks, or custom intervals.
```
Groups:  [Day1] [Day2] [Day3] [Day4] [Day5]

Fold1: Train -> Day2-Day5 | Test -> Day1

Fold2: Train -> Day1,Day3-Day5 | Test -> Day2

Fold3: Train -> Day1-Day2,Day4-Day5 | Test -> Day3
...
```
```csharp
var ts = TimeSeries.FromCsv("data.csv");

var cv = new LeaveOneOutCrossValidator(pipelineGrid);
var result = cv.Run(ts, "Target", new DailyGrouping());
```
**Key points:**

* Order matters
* Leakage must be avoided
* Grouping often represents **time intervals**

| Validator                       | Uses grouping | Temporal awareness | Notes                                                                                 |
| ------------------------------- | ------------- | ------------------ | ------------------------------------------------------------------------------------- |
| `KFoldCrossValidator`           | ❌             | ❌                  | Classic tabular K-Fold; all samples used exactly once.                                |
| `LeaveOneOutCrossValidator`     | ✅ (optional)  | ❌                  | Extreme case of K-Fold; can act as Leave-One-Group-Out if groups are provided.        |
| `RollingCrossValidator`         | ✅ (implicit)  | ✅                  | Designed for time series; respects temporal order to prevent leakage.                 |
| `ShuffleSplitCrossValidator`    | ❌             | ❌                  | Random train/test splits; multiple iterations; not all rows guaranteed to be tested.  |
| `MonteCarloCrossValidator`      | ❌             | ❌                  | Like ShuffleSplit but returns full score distribution with CI, histogram & convergence. |
| `StratifiedKFoldCrossValidator` | ❌             | ❌                  | Maintains class proportions; only for classification; useful for imbalanced datasets. |


---

## 🧪 SupervisedExperiment (Fluent API)

A high-level fluent API for **supervised ML experiments** — the supervised counterpart to `ClusteringExperiment`. Combines a `PipelineGrid` with one or more **cross-validation strategies** in a single call, then ranks all pipeline × CV combinations automatically.

### ➡️ Quick Start

```csharp
var result = SupervisedExperiment
    .For(X, y)
    .WithGrid(new PipelineGrid()
        .AddModel<KNearestNeighbors>(g => g
            .Add("K", 3, 5, 7)
            .AddScaler<StandardScaler>(s => { }))
        .AddModel<DecisionTree>(g => g
            .Add("MaxDepth", 3, 5, 10)))
    .WithCrossValidator(CrossValidatorConfig.KFold(folds: 5))
    .Run();

Console.WriteLine(result.BestModelName);   // e.g. "DecisionTree"
Console.WriteLine(result.BestScore);       // e.g. 0.92
Pipeline best = result.BestPipeline;
```

### ➡️ Multiple Cross-Validators

Add several cross-validation strategies — the **first one added is the primary** used for final ranking (same convention as evaluators in `ClusteringExperiment`):

```csharp
var result = SupervisedExperiment
    .For(X, y)
    .WithGrid(new PipelineGrid()
        .AddModel<KNearestNeighbors>(g => g
            .Add("K", 3, 5, 7)
            .AddScaler<StandardScaler>(s => { }))
        .AddModel<DecisionTree>(g => g
            .Add("MaxDepth", 3, 5, 10))
        .AddModel<RandomForest>(g => g
            .Add("NumTrees", 50, 100)
            .Add("MaxDepth", 5, 8))
        .AddModel<LinearSVC>(g => g
            .Add("C", 0.1, 1.0, 10.0)
            .Add("LearningRate", 0.001, 0.01)
            .Add("Epochs", 500, 1000)))
    .WithCrossValidators(
        CrossValidatorConfig.KFold(folds: 5),
        CrossValidatorConfig.StratifiedKFold(folds: 10),
        CrossValidatorConfig.ShuffleSplit(nSplits: 10, testSize: 0.2))
    .Run();

// Ranked by primary CV (KFold)
foreach (var r in result.Rankings.Take(5))
{
    Console.WriteLine(
        $"{r.ModelName,-25} " +
        $"KFold={r.Scores["KFold"],6:F3}  " +
        $"Stratified={r.Scores["StratifiedKFold"],6:F3}  " +
        $"Shuffle={r.Scores["ShuffleSplit"],6:F3}");
}

// Best by a specific cross-validator
var bestStratified = result.BestBy("StratifiedKFold");
```

### ➡️ Regression Example

```csharp
var result = SupervisedExperiment
    .For(X, y)
    .WithGrid(new PipelineGrid()
        .AddModel<Linear>(g => { })
        .AddModel<Ridge>(g => g
            .Add("Alpha", 0.01, 0.1, 1.0))
        .AddModel<Lasso>(g => g
            .Add("Alpha", 0.01, 0.1))
        .AddModel<ElasticNet>(g => g
            .Add("Lambda", 0.01, 0.1)
            .Add("L1Ratio", 0.3, 0.5, 0.7)))
    .WithCrossValidators(
        CrossValidatorConfig.KFold(folds: 5),
        CrossValidatorConfig.MonteCarlo(iterations: 100, testSize: 0.2, seed: 42))
    .Run();

Console.WriteLine($"Best: {result.BestModelName} → R² = {result.BestR2:F4}");
```

### ➡️ Simple Pipeline Mode

Instead of a grid, you can also pass manually-constructed pipelines:

```csharp
var result = SupervisedExperiment
    .For(X, y)
    .WithPipelines(
        new Pipeline(new Linear(), new Dictionary<string, object>()),
        new Pipeline(new Ridge(), new Dictionary<string, object> { { "Alpha", 0.1 } }))
    .WithCrossValidator(CrossValidatorConfig.KFold(folds: 5))
    .Run();
```

### 📐 Available Cross-Validator Configs

| Factory method | CV strategy | Notes |
|---|---|---|
| `CrossValidatorConfig.KFold(folds)` | K-Fold | Standard tabular CV |
| `CrossValidatorConfig.StratifiedKFold(folds)` | Stratified K-Fold | Classification only, preserves class ratios |
| `CrossValidatorConfig.ShuffleSplit(nSplits, testSize, ...)` | Shuffle-Split | Random train/test splits |
| `CrossValidatorConfig.LeaveOneOut()` | Leave-One-Out | K = n, expensive |
| `CrossValidatorConfig.MonteCarlo(iterations, testSize, seed, ...)` | Monte Carlo CV | Full score distributions with CI |
| `CrossValidatorConfig.Custom(name, factory)` | Any `ICrossValidator` | Wrap your own implementation |

### 📐 Result Properties

| Property | Description |
|---|---|
| `Rankings` | All pipelines ranked by primary CV score (descending) |
| `Best` / `BestPipeline` / `BestModelName` / `BestScore` | Convenience accessors for the top-ranked pipeline |
| `BestBy(cvName)` | Best pipeline according to a specific cross-validator |
| `CVResults[cvName]` | Raw `CrossValidationResult` per CV (confusion matrix, R², actual/predicted) |
| `BestConfusionMatrix` | Confusion matrix from the best pipeline's primary CV (classification) |
| `BestR2` | Coefficient of determination from the best pipeline's primary CV (regression) |
| `BestActualValues` / `BestPredictedValues` | Actual vs predicted from the best pipeline |
| `ScoreSummary(cvName?)` | Descriptive statistics across all pipeline scores (see below) |
| `ScorePercentile(result, cvName?)` | Percentile rank (0–100) of a specific pipeline |
| `RankCorrelation(cv1, cv2)` | Spearman ρ between rankings from two CVs |
| `ScoreConsistency(result)` | StdDev of a pipeline's scores across CVs |
| `TotalDuration` | Wall-clock time for the entire experiment |

**Key points:**

* Reuses existing `PipelineGrid` — same `AddModel<T>()`, `AddScaler<T>()`, `AddSelector<T>()` syntax
* Reuses all existing `ICrossValidator` implementations — no new CV code needed
* First cross-validator added is the **primary** for ranking (matches ClusteringExperiment convention)
* `CrossValidatorConfig` is a **lazy factory** — pipelines are bound at `Run()` time
* Each pipeline is scored by **all** cross-validators for multi-perspective evaluation

### 📊 Score Distribution Statistics (Supervised & Clustering)

Both `SupervisedExperimentResult` and `ClusteringExperimentResult` expose **descriptive statistics** computed from the library's `DescriptiveStatisticsExtensions`. These help answer: *"how sensitive is my score to the pipeline configuration?"*

**`ScoreSummary()`** — overview of all pipeline scores for a given CV/evaluator:

```csharp
var result = SupervisedExperiment
    .For(X, y)
    .WithGrid(grid)
    .WithCrossValidator(CrossValidatorConfig.KFold(folds: 5))
    .Run();

var summary = result.ScoreSummary();          // primary CV
var summary2 = result.ScoreSummary("ShuffleSplit"); // specific CV

Console.WriteLine(summary);
// N=8, Mean=0.8520, Median=0.8650, StdDev=0.0430, IQR=0.0510,
// Range=[0.7800, 0.9200], Skew=-0.320, Kurt=-0.810
```

| Property | Meaning |
|---|---|
| `Mean` / `Median` | Central tendency of scores across all tested pipelines |
| `StandardDeviation` | How spread out the scores are — high = very model-sensitive problem |
| `InterquartileRange` | Middle 50 % spread (robust to outliers) |
| `Skewness` | Positive = few good pipelines; Negative = most pipelines are good |
| `Kurtosis` | High = score outliers (a few configs are much better/worse) |
| `ConfidenceInterval` | 95 % CI for the mean score |
| `Range` / `Min` / `Max` | Full score span |

**`ScorePercentile()`** — how a specific pipeline ranks relative to all others:

```csharp
double pct = result.ScorePercentile(result.Best);       // e.g. 100.0
double pct2 = result.ScorePercentile(result.Rankings[3]); // e.g. 62.5
```

**`RankCorrelation()`** *(supervised only)* — Spearman ρ between two CV rankings:

```csharp
var result = SupervisedExperiment.For(X, y)
    .WithGrid(grid)
    .WithCrossValidators(
        CrossValidatorConfig.KFold(folds: 5),
        CrossValidatorConfig.StratifiedKFold(folds: 10))
    .Run();

var (rho, pValue) = result.RankCorrelation("KFold", "StratifiedKFold");
// rho ≈ 0.94, p < 0.01 → CVs agree strongly on ranking
```

**`ScoreConsistency()`** *(supervised only)* — StdDev of a pipeline's score across CVs:

```csharp
double consistency = result.ScoreConsistency(result.Best);
// Low → model performs the same regardless of CV strategy
// High → performance is CV-dependent (potential overfitting signal)
```

The same `ScoreSummary()` and `ScorePercentile()` are available on `ClusteringExperimentResult`:

```csharp
var clustResult = ClusteringExperiment.For(X)
    .WithAlgorithm(new KMeans())
    .TryClusterCounts(2, 10)
    .WithEvaluators(new SilhouetteEvaluator(), new CalinskiHarabaszEvaluator())
    .Run();

var silSummary = clustResult.ScoreSummary("Silhouette");
var chSummary  = clustResult.ScoreSummary("CalinskiHarabasz");
double pct = clustResult.ScorePercentile(clustResult.Best);
```

---

## 📊 Classification Models

All classifiers implement `IClassificationModel` and operate directly on `Matrix` and `Vector` primitives.

**Logistic Regression**

Class: `Logistic`

Hyperparameters:

* `LearningRate`
* `MaxIterations`
* `FitIntercept` 
* `RegularizationStrength`
* `Tolerance` 


**Decision Tree (Classifier)**

Class: `DecisionTree`

Hyperparameters:

* `MaxDepth`
* `MinSamplesSplit`


**Random Forest**

Class: `RandomForest`

Hyperparameters:

* `NumTrees`
* `MaxDepth`
* `MinSamplesSplit`



**K-Nearest Neighbors**

Class: `KNearestNeighbors`

Hyperparameters:

* `K`

**Naive Bayes**

Class: `NaiveBayes`

Hyperparameters:
(No tunable hyperparameters)

**Support Vector Classifier (Linear)**

Class: `LinearSVC`

Hyperparameters:

* `C` (regularization strength)
* `LearningRate`
* `Epochs`

**Support Vector Classifier (Kernel)**

Class: `KernelSVC`

Hyperparameters:

* `C`
* `Kernel` (RBF, Polynomial)
* `LearningRate`
* `Epochs`
* `Gamma`
* `Degree` (for polynomial kernel)

**Multilayer Perceptron (Classifier)**

Class: `MLPClassifier`

Hyperparameters:

* `HiddenLayers` (e.g. `64`, `64,32`)
* `LearningRate`
* `Epochs`
* `Activation` (ReLU, Tanh, Sigmoid)

---

## 📈 Regression Models

All regressors implement `IRegressionModel`.


**Linear**

Class: `Linear`

Hyperparameters:

* `LearningRate`
* `FitIntercept`

**Ridge Regression (L2)**

Class: `Ridge`

Hyperparameters:

* `Alpha`
* `FitIntercept`


**Lasso Regression (L1)**

Class:`Lasso`

Hyperparameters:

* `Alpha`
* `MaxIterations`

**Elastic Net (L1 + L2)**

Class: `ElasticNet`

Hyperparameters:

* `Lambda`
* `L1Ratio`

**Support Vector Regression (Linear)**

Class: `LinearSVR`

Hyperparameters:

* `C`
* `Epsilon`
* `LearningRate`
* `Epochs`

**Support Vector Regression (Kernel)**

Class: `KernelSVR`

Hyperparameters:

* `C`
* `LearningRate`
* `Epochs`
* `Kernel`
* `Gamma`
* `Degree`


**Multilayer Perceptron (Regressor)**

Class: `MLPRegressor`

Hyperparameters:

* `HiddenLayers`
* `LearningRate`
* `Epochs`
* `BatchSize`
* `L2`
* `Activation`

---

## 🔬 Clustering (Unsupervised)

CSharpNumerics includes a full **unsupervised clustering framework** with the same philosophy as the supervised pipeline: easy API, pluggable algorithms, and transparent results.

All clustering models implement `IClusteringModel` and operate directly on `Matrix` primitives — no target vector needed.

---

### ➡️ ClusteringExperiment (Fluent API)

The simplest way to try clustering — just pick an algorithm, a range of K values, and an evaluator:

```csharp
var experiment = ClusteringExperiment
    .For(X)
    .WithAlgorithm(new KMeans())
    .TryClusterCounts(2, 10)
    .WithEvaluator(new SilhouetteEvaluator())
    .WithScaler(new StandardScaler())
    .Run();

Console.WriteLine(experiment.BestClusterCount);   // e.g. 4
Console.WriteLine(experiment.BestScore);           // e.g. 0.72
VectorN labels = experiment.BestLabels;
```

**Multi-algorithm comparison**

```csharp
var experiment = ClusteringExperiment
    .For(X)
    .WithAlgorithms(new KMeans(), new AgglomerativeClustering())
    .TryClusterCounts(2, 8)
    .WithEvaluators(new SilhouetteEvaluator(), new CalinskiHarabaszEvaluator())
    .WithScaler(new StandardScaler())
    .Run();

foreach (var r in experiment.Rankings)
    Console.WriteLine($"{r.AlgorithmName} K={r.ClusterCount} → {r.Scores["Silhouette"]:F3}");

var bestSilhouette = experiment.BestBy<SilhouetteEvaluator>();
var bestCH = experiment.BestBy<CalinskiHarabaszEvaluator>();
```

**Key points:**

* `TryClusterCounts(min, max)` auto-expands K for algorithms that accept it
* DBSCAN discovers K on its own — the range is **ignored**
* First evaluator added is the **primary** used for ranking
* `BestBy<T>()` retrieves the best result by a specific evaluator

---

### ➡️ Clustering Grid

Full hyperparameter grid search across multiple algorithms and scalers:

```csharp
var experiment = ClusteringExperiment
    .For(X)
    .WithGrid(new ClusteringGrid()
        .AddModel<KMeans>(g => g
            .Add("K", 2, 3, 4, 5, 6, 7, 8)
            .Add("InitMethod", KMeansInit.Random, KMeansInit.PlusPlus)
            .AddScaler<StandardScaler>(s => { }))
        .AddModel<DBSCAN>(g => g
            .Add("Epsilon", 0.3, 0.5, 0.8, 1.0, 1.5)
            .Add("MinPoints", 3, 5, 10)
            .AddScaler<MinMaxScaler>(s => { }))
        .AddModel<AgglomerativeClustering>(g => g
            .Add("K", 2, 3, 4, 5)
            .Add("Linkage", LinkageType.Ward, LinkageType.Complete)))
    .WithEvaluators(
        new SilhouetteEvaluator(),
        new CalinskiHarabaszEvaluator(),
        new DaviesBouldinEvaluator())
    .Run();

foreach (var r in experiment.Rankings.Take(10))
{
    Console.WriteLine(
        $"{r.AlgorithmName,-25} K={r.ClusterCount,2}  " +
        $"Silhouette={r.Scores["Silhouette"],7:F3}  " +
        $"CH={r.Scores["CalinskiHarabasz"],10:F1}  " +
        $"({r.Duration.TotalMilliseconds:F0} ms)");
}
```

---

### 🔵 Clustering Models

**KMeans**

Class: `KMeans`

Hyperparameters:

* `K`
* `MaxIterations`
* `Tolerance`
* `Seed`
* `InitMethod` (Random, PlusPlus)

Exposes after fit: `Centroids`, `Inertia`, `Iterations`

**DBSCAN**

Class: `DBSCAN`

Hyperparameters:

* `Epsilon`
* `MinPoints`

Discovers K automatically. Noise points labeled `-1`. Exposes `NoiseCount`.

**Agglomerative Clustering**

Class: `AgglomerativeClustering`

Hyperparameters:

* `K`
* `Linkage` (Single, Complete, Average, Ward)

Bottom-up hierarchical merging. Exposes `Dendrogram`.

---

### 📐 Clustering Evaluators

All evaluators implement `IClusteringEvaluator` where **higher score = better**.

| Evaluator | Class | Metric | Notes |
|-----------|-------|--------|-------|
| Silhouette | `SilhouetteEvaluator` | $s \in [-1, 1]$ | Higher = better |
| Inertia (Elbow) | `InertiaEvaluator` | $-W$ (negated) | Use `RawInertia()` for elbow curve |
| Davies-Bouldin | `DaviesBouldinEvaluator` | $-DB$ (negated) | Lower DB = better separation |
| Calinski-Harabasz | `CalinskiHarabaszEvaluator` | $CH$ | Higher = better, fast |

---

### 🎲 Monte Carlo Clustering (Uncertainty Estimation)

Use the `MonteCarloClustering` class to quantify **how stable** your clustering results are via bootstrap resampling. Two analysis modes are available:

#### Bootstrap — Consensus Matrix & Score Distribution

Runs the algorithm many times on bootstrap-resampled data. Produces a **consensus matrix**, per-point **stability scores**, and full **score distributions** with confidence intervals.

```csharp
var mc = new MonteCarloClustering { Iterations = 200, Seed = 42 };
var result = mc.RunBootstrap(
    data,
    new KMeans { K = 3 },
    new SilhouetteEvaluator(),
    new StandardScaler());   // optional

// Score uncertainty
var ci = result.ScoreConfidenceInterval();       // e.g. (0.68, 0.74)
double se = result.ScoreDistribution.StandardError;
var histogram = result.ScoreDistribution.Histogram(20);

// Consensus matrix (N × N) — fraction of times each pair co-clustered
Matrix consensus = result.ConsensusMatrix;

// Per-point stability [0, 1] — how consistently each point stays in its cluster
double[] stability = result.PointStability;
double[] convergence = result.ConvergenceCurve;  // running mean of score
```

#### Experiment — Optimal-K Distribution

Runs a full K-range experiment many times on bootstrap samples. Shows **how often each K value is selected as best**, revealing whether the optimal K is robust.

```csharp
var mc = new MonteCarloClustering { Iterations = 100, Seed = 42 };
var kResult = mc.RunExperiment(
    data,
    new KMeans(),
    new SilhouetteEvaluator(),
    minK: 2, maxK: 8);

// Which K values won across the 100 bootstrap runs?
foreach (var (k, count) in kResult.OptimalKDistribution.OrderByDescending(x => x.Value))
    Console.WriteLine($"K={k}: chosen {count}/100 times");

// Score distribution for the best K in each iteration
var ci = kResult.ScoreConfidenceInterval();
```

#### Fluent API Integration

Add Monte Carlo uncertainty with a single builder call:

```csharp
var result = ClusteringExperiment
    .For(data)
    .WithAlgorithm(new KMeans())
    .TryClusterCounts(2, 8)
    .WithEvaluator(new SilhouetteEvaluator())
    .WithScaler(new StandardScaler())
    .WithMonteCarloUncertainty(iterations: 200, seed: 42)
    .Run();

// Standard result
Console.WriteLine($"Best K = {result.BestClusterCount}");

// Monte Carlo result (populated automatically)
var mcResult = result.MonteCarloResult;
Console.WriteLine($"Score CI = {mcResult.ScoreConfidenceInterval()}");
Console.WriteLine($"K distribution: {string.Join(", ",
    mcResult.OptimalKDistribution.Select(kv => $"K={kv.Key}: {kv.Value}"))}");
```

**Key points:**

* Bootstrap uses **sampling with replacement** — each iteration sees ~63 % unique points
* Consensus matrix cell (i,j) = fraction of runs where points i and j co-clustered (normalized by co-occurrence)
* Point stability = average consensus with same-cluster neighbours; close to 1.0 = very stable
* `RunExperiment` requires a K-accepting algorithm (KMeans, AgglomerativeClustering) — uses reflection to set K
* All results include full `MonteCarloResult` from the statistics engine (Mean, StdDev, Percentile, Histogram, CI, StandardError)
* Reproducible when `Seed` is set

---

## 🔻 Dimensionality Reduction

CSharpNumerics includes **unsupervised dimensionality reduction** as an optional preprocessing step in both supervised and clustering pipelines. Reducers implement `IDimensionalityReducer` and slot into the pipeline between feature selection and scaling.

### Pipeline Order

| Pipeline | Order |
|----------|-------|
| **Supervised** | Selector → **Reducer** → Scaler → Model |
| **Clustering** | **Reducer** → Scaler → Model |

### ➡️ Interface

All reducers implement `IDimensionalityReducer`:

```csharp
public interface IDimensionalityReducer
{
    int NComponents { get; set; }
    Matrix FitTransform(Matrix X);
    Matrix Transform(Matrix X);
    IDimensionalityReducer Clone();
}
```

### 🔵 Algorithms

**Principal Component Analysis (PCA)**

Class: `PCA`

Projects data onto the top eigenvectors of the covariance matrix.
Uses power iteration with deflation for eigendecomposition.

Hyperparameters:

* `NComponents` — number of output dimensions
* `MaxIterations` — power iteration limit (default 1000)
* `Tolerance` — convergence threshold (default 1e-8)
* `Seed` — optional random seed

Exposes after fit: `Components`, `ExplainedVariance`, `ExplainedVarianceRatio`, `Mean`

### ➡️ Standalone Usage

```csharp
var pca = new PCA { NComponents = 2 };
Matrix reduced = pca.FitTransform(X);  // n × 2

// Inspect explained variance
for (int i = 0; i < pca.ExplainedVarianceRatio.Length; i++)
    Console.WriteLine($"PC{i + 1}: {pca.ExplainedVarianceRatio[i]:P1}");

// Transform new data
Matrix newReduced = pca.Transform(Xtest);
```

### ➡️ Clustering Pipeline Integration

```csharp
var experiment = ClusteringExperiment
    .For(X)
    .WithAlgorithm(new KMeans())
    .TryClusterCounts(2, 8)
    .WithEvaluator(new SilhouetteEvaluator())
    .WithReducer(new PCA { NComponents = 5 })
    .WithScaler(new StandardScaler())
    .Run();

Console.WriteLine(experiment.BestClusterCount);
```

### ➡️ Clustering Grid Integration

```csharp
var experiment = ClusteringExperiment
    .For(X)
    .WithGrid(new ClusteringGrid()
        .AddModel<KMeans>(g => g
            .Add("K", 2, 3, 4, 5)
            .AddReducer<PCA>(r => r.Add("NComponents", 2, 5, 10))
            .AddScaler<StandardScaler>(s => { })))
    .WithEvaluator(new SilhouetteEvaluator())
    .Run();
```

### ➡️ Supervised Pipeline Integration

```csharp
var result = SupervisedExperiment
    .For(X, y)
    .WithGrid(new PipelineGrid()
        .AddModel<KNearestNeighbors>(g => g
            .Add("K", 3, 5, 7)
            .AddReducer<PCA>(r => r.Add("NComponents", 2, 5))
            .AddScaler<StandardScaler>(s => { }))
        .AddModel<DecisionTree>(g => g
            .Add("MaxDepth", 3, 5, 10)))
    .WithCrossValidator(CrossValidatorConfig.KFold(folds: 5))
    .Run();
```

**Key points:**

* Reducers are **optional** — existing pipelines work unchanged
* PCA uses power iteration — no external dependencies
* `ExplainedVarianceRatio` shows how much variance each component captures
* Grid search over `NComponents` finds the optimal dimensionality automatically
* Works with both supervised and clustering pipelines
* Follows the same `FitTransform` / `Transform` / `Clone` pattern as scalers
* Implements `IHasHyperparameters` for grid search integration

