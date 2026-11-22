# Sannolikhetsklassificering & utvärdering — Roadmap

> **Status 2026-11-22:** alla faser klara på grenen `feat/ml-probabilistic-classification`.
> Flyttas till `docs/completed/` efter merge.

## Mål

Göra ML-delen användbar för klassificering där **sannolikheten** spelar roll, inte bara klassen, och där
datan har **grupper** som inte får läcka mellan tränings- och testdata. Drivande användningsfall:
vetting av exoplanetkandidater i ExoplanetHunter Transit Lab (ungefär 12 000 exempel från ~2 700
stjärnor, 15 features, en modell som ska ge "sannolikhet att signalen är en planet").

## Nulägesanalys

| Komponent | Status | Brist |
|---|---|---|
| `IClassificationModel` | ✓ | `Predict` ger bara klass. `Logistic`, `MLPClassifier` och sekvensklassificerarna räknar ut sannolikheter internt men kastar dem |
| `DecisionTree`, `RandomForest` | ✓ | Varje nod provar varje unikt värde som tröskel och kopierar matrisen per försök (O(n²·p) i roten). Ingen slumpad feature-delmängd, inga vikter, inga sannolikheter |
| Korsvalidering | ✓ KFold, Stratified, ShuffleSplit, MonteCarlo, LOO | Gruppering bara via LeaveOneOut (en fold per grupp) eller tidsgrupper. Ingen GroupKFold |
| Mått | Accuracy, R², förväxlingsmatris | Ingen ROC-AUC, PR-AUC, Brier, log loss eller kalibrering |
| Vikter | ✗ | Ingen modell tar provvikter; obalanserade klasser går inte att väga |

## Designval

- **`IProbabilisticClassifier : IClassificationModel`** med `Matrix PredictProbabilities(Matrix X)`
  (n × NumClasses). Ett nytt gränssnitt i stället för en ändring av `IClassificationModel`, så att
  befintliga och externa modeller fortsätter att kompilera.
- **`ISampleWeightedModel : IModel`** med `Fit(Matrix X, VectorN y, VectorN sampleWeights)`.
  `ClassWeights.Balanced(y)` ger vikter som väger upp små klasser.
- **Grupperad korsvalidering** får grupperna i konstruktorn (en grupp-id per rad), precis som
  pipelines, och implementerar befintliga `ICrossValidator`. Då fungerar de i `SupervisedExperiment`
  via `CrossValidatorConfig.GroupKFold(groups)`.
- **Mått i `CSharpNumerics.ML.Evaluation`** som statiska funktioner på `VectorN`, med valfria vikter.
- **Träd:** sorterad tröskelsökning över indexlistor (O(n log n) per feature och nod), vägd Gini,
  sannolikheter i löven, `MaxFeatures` och `MinSamplesLeaf`. Utan vikter och med standardvärden
  väljer trädet samma kandidattrösklar som tidigare.

## Phase 1 — Sannolikheter och vikter i grundmodellerna

- [x] `IProbabilisticClassifier`, `ISampleWeightedModel`, `ClassWeights.Balanced`
- [x] `DecisionTree`: indexbaserad, sorterad split-sökning, vägd Gini, lövsannolikheter,
      `MaxFeatures`, `MinSamplesLeaf`, `Seed`
- [x] `RandomForest`: vägd bootstrap utan matriskopior, parallell trädträning, medelvärde av
      lövsannolikheter, `MaxFeatures`, `MinSamplesLeaf`, `Seed`
- [x] `Logistic`: `PredictProbabilities`, vägd gradient
- [x] `MLPClassifier`: `PredictProbabilities`, vägd förlust
- [x] Sekvensklassificerarna (`CNN1D`, `LSTM`, `BiLSTM`, `TCN`): `PredictProbabilities`
- [x] `Pipeline.Fit(X, y, sampleWeights)` och `Pipeline.PredictProbabilities`
- [x] Tester

## Phase 2 — Mått

- [x] `ClassificationMetrics`: ROC-kurva och ROC-AUC (med lika poäng hanterade), PR-kurva och
      average precision, precision vid given recall, Brier-poäng, log loss, kalibreringskurva
- [x] Tester mot handräknade värden

## Phase 3 — Grupperad korsvalidering

- [x] `GroupKFoldCrossValidator`: grupper fördelas efter storlek på den minst fyllda folden
- [x] `StratifiedGroupKFoldCrossValidator`: grupper fördelas så att klassandelarna blir så lika
      som möjligt mellan folds
- [x] Valfria provvikter och poängsättning (`ClassificationScoring`: Accuracy, RocAuc,
      AveragePrecision, NegativeLogLoss, NegativeBrierScore)
- [x] `CrossValidationResult`: sannolikheter utanför folden (`PredictedProbabilities`) och
      radindex (`RowIndices`)
- [x] `CrossValidatorConfig.GroupKFold` / `StratifiedGroupKFold`
- [x] Tester: ingen grupp i både tränings- och testdata, balans, determinism

## Phase 4 — Dokumentation

- [x] ML-README: sannolikheter, vikter, grupperad korsvalidering, mått
