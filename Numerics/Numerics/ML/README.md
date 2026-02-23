
CSharpNumerics includes a **lightweight, fully numerical machine learning framework** designed for **research, experimentation, and educational use**.
The focus is on **transparency**, **mathematical clarity**, and **pipeline-based model evaluation** — not black-box automation.
## ➡️ Pipeline Grid
All models are implemented directly on top of the library’s **Matrix** and **Vector** primitives.

Models can be combined with:

* **Scalers** (e.g. `StandardScaler`)
* **Feature selectors** (e.g. `SelectKBest`)
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
| `StratifiedKFoldCrossValidator` | ❌             | ❌                  | Maintains class proportions; only for classification; useful for imbalanced datasets. |


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

