# ML-expansion — Gradient boosting, attention & ONNX

## Mål

Fylla de mest efterfrågade luckorna i ML-delen: **gradient boosting** (den dominerande algoritmfamiljen på tabulär data), **attention/Transformer-lager** som komplement till LSTM/TCN i sekvensmodulen, **ONNX-export/import** för produktionsinterop, samt **t-SNE/UMAP** för visualisering. Tillsammans med autodiff-arbetet (se AutodiffRoadMap) lyfter detta ML-delen från "bred" till "konkurrenskraftig".

---

## Nulägesanalys — Befintlig arkitektur

### Vad som finns idag

| Komponent | Status | Plats | Kommentar |
|-----------|--------|-------|-----------|
| `DecisionTree`, `RandomForest` | ✓ | `ML/Models/` | Trädinfrastruktur att återanvända för boosting |
| `IModel`/`IHasHyperparameters` | ✓ | `ML/Models/Interfaces/` | Grid-search-kompatibilitet |
| `SupervisedExperiment`, `Pipeline` | ✓ | `ML/Experiment/`, `ML/Pipeline.cs` | Nya modeller pluggar in direkt |
| Sekvenslager (LSTM, BiLSTM, Conv1D, TCN) | ✓ | `ML/Sequence/` | `ILayer` + `SequentialModel`-arkitektur |
| `ISequenceModel`, 3D-datakonvention | ✓ | `ML/Sequence/` | Etablerad i SequenceModelsRoadMap |
| `PCA` | ✓ | `ML/` | Enda dimensionsreduktionen |
| Gradient boosting | ✗ | — | Största luckan i ML-delen |
| Attention/Transformer | ✗ | — | Finns ej |
| ONNX-export/import | ✗ | — | Modeller är inlåsta i biblioteket |
| t-SNE/UMAP | ✗ | — | Finns ej |
| `NaiveBayes.NumClasses` | ⚠ | `ML/Models/` | Kastar `NotImplementedException` — städas |

---

## Del 1 — Gradient boosting

### Vad det gör

Bygger ett ensemble av grunda regressionsträd sekventiellt, där varje träd anpassas mot förlustens gradient (och Hessian) för nuvarande prediktion — XGBoost/LightGBM-familjen. Vinner fortfarande de flesta tabulära benchmarks och är den mest efterfrågade enskilda algoritmen.

### Vad som krävs att bygga

| Komponent | Beskrivning | Insats |
|-----------|-------------|--------|
| **`RegressionTreeBuilder`** | Träd som splittar på gradient/Hessian-statistik (inte Gini/varians). Kan dela kod med `DecisionTree` men behöver egen split-logik | **Medel** |
| **`GradientBoostingRegressor`** | Squared/absolute/Huber-loss, learning rate (shrinkage), subsampling, early stopping mot valideringsset | **Medel** |
| **`GradientBoostingClassifier`** | Log-loss (binär) + softmax (multiklass, ett träd per klass och runda) | **Medel** |
| **Histogram-baserade splits** | Binna features (256 bins) → O(bins) i stället för O(samples) per split — LightGBM-tricket som ger 10–100× | **Medel-hög** |
| **Regularisering** | `lambda` (L2 på lövvikter), `gamma` (min split gain), max depth/leaves — XGBoost-stil | **Låg** (del av split-formeln) |
| **Feature importance** | Gain-baserad + split-count | **Låg** |
| **`IModel`-integration** | `IHasHyperparameters` → grid-search över lr/depth/estimators fungerar direkt i `SupervisedExperiment` | **Låg** |

### Genomförbarhet

**HÖG.** Trädbyggnads-infrastruktur finns i `DecisionTree`/`RandomForest`. XGBoost-split-formeln (gain = gradientsummor²/Hessiansummor) är väldokumenterad. Börja exakt-greedy (enkel, korrekt), lägg histogram-optimering efteråt med exakt-varianten som facit.

---

## Del 2 — Attention & Transformer-lager

### Vad det gör

Scaled dot-product attention låter varje position i en sekvens vikta information från alla andra positioner:

$$\text{Attention}(Q,K,V) = \text{softmax}\!\left(\frac{QK^\top}{\sqrt{d_k}}\right)V$$

Multi-head attention + feedforward + residual/LayerNorm = Transformer-encoderblock. Komplement till LSTM/TCN för långa beroenden, och grunden för modern sekvensmodellering.

### Vad som krävs att bygga

| Komponent | Beskrivning | Beroenden | Insats |
|-----------|-------------|-----------|--------|
| **`SelfAttentionLayer`** | Q/K/V-projektioner, scaled dot-product, softmax över tidsaxeln. `ILayer`-implementering | `Matrix`-mult, softmax | **Medel** |
| **`MultiHeadAttentionLayer`** | h parallella huvuden, konkatenering + utprojektion | `SelfAttentionLayer` | **Medel** |
| **`LayerNorm`** | Normalisering per position (jfr befintlig `BatchNorm1D`) | — | **Låg** |
| **`PositionalEncoding`** | Sinus/cosinus-encoding adderad på input | — | **Trivial** |
| **`TransformerEncoderBlock`** | MHA → Add&Norm → FFN → Add&Norm, komposition av ovanstående | Allt ovan + `Residual` (finns) | **Låg-medel** |
| **`TransformerClassifier`/`TransformerRegressor`** | `ISequenceModel`-implementeringar i samma stil som `LSTMClassifier` | `SequentialModel` | **Låg** |
| **Backward genom attention** | Handskriven: softmax-Jacobian + matmul-kedjor. **Eller: vänta på autodiff (Phase 4 i AutodiffRoadMap) och slippa** | — | **Hög** (handskriven) / **Låg** (autodiff) |

### Sekvensering mot autodiff

Attention-backward för hand är den enskilt felbenägnaste biten. **Rekommendation:** implementera forward + handskriven backward för `LayerNorm`/`PositionalEncoding` (enkla), men lås attention-lagret till autodiff-träningsvägen om tidplanen tillåter — annars exakt-greedy handskriven backward med gradient-check-tester mot finita differenser.

---

## Del 3 — ONNX-export/import

### Varför

ONNX gör modeller tränade i CSharpNumerics körbara i ONNX Runtime (produktion, mobil, webb) och tvärtom — den största enskilda adoption-drivaren för ett .NET-ML-bibliotek. Utan interop är varje tränad modell inlåst.

### Vad som krävs att bygga

| Komponent | Beskrivning | Insats |
|-----------|-------------|--------|
| **Beslut: protobuf-beroende** | ONNX är protobuf-baserat. Alternativ A: `Google.Protobuf` + genererade ONNX-scheman i ett **separat paket** `CSharpNumerics.Onnx` (kärnbiblioteket förblir beroendefritt). Alternativ B: handskriven minimal protobuf-writer (inga beroenden, mer jobb). **Rekommendation: A, separat paket** | — |
| **Export: MLP/Sequential** | Dense/ReLU/Sigmoid/Tanh → Gemm/aktiverings-noder | **Medel** |
| **Export: Conv1D/pooling/LSTM** | Mappning till ONNX Conv/MaxPool/LSTM-operatorer (LSTM-vikternas layout är pillrig) | **Medel-hög** |
| **Export: träd-ensembler** | RandomForest/GradientBoosting → `ai.onnx.ml` TreeEnsemble-operatorer | **Medel** |
| **Export: skalers/pipeline** | MinMax/Standard → Scaler-noder så hela `Pipeline` exporteras, inte bara modellen | **Låg-medel** |
| **Import (begränsad)** | Läsa in MLP-vikter från ONNX för finjustering — full generell import är utom scope | **Medel** |
| **Verifiering** | Round-trip-tester: prediktioner i CSharpNumerics ≡ ONNX Runtime på samma input (tolerans 1e-6) | **Låg** (kräver ONNX Runtime som test-beroende) |

---

## Del 4 — t-SNE & UMAP

| Komponent | Beskrivning | Insats |
|-----------|-------------|--------|
| **`TSne`** | Barnes–Hut t-SNE (kvadtree-approximation) med perplexity-sökning | **Medel-hög** |
| **`Umap`** | k-NN-graf + fuzzy simplicial set + SGD-layout | **Hög** |
| **Gemensamt `IReducer`-interface** | Samma kontrakt som `PCA` → pluggar in i `Pipeline` | **Låg** |

Lägst prioritet av de fyra delarna — värdefullt men inte differentierande. Exakt (icke-Barnes–Hut) t-SNE är en acceptabel första version för dataset < 5 000 punkter.

---

## Implementationsplan — Faser

### Phase 1 — Städning & boosting-grund
- [ ] Fixa `NaiveBayes.NumClasses`-stubben
- [ ] `RegressionTreeBuilder` med gradient/Hessian-splits (exakt greedy)
- [ ] `GradientBoostingRegressor` (squared loss, shrinkage, subsampling, early stopping)
- [ ] Enhetstester + jämförelse mot RandomForest på syntetiska dataset

### Phase 2 — Boosting komplett
- [ ] `GradientBoostingClassifier` (log-loss binär → softmax multiklass)
- [ ] Regularisering (lambda, gamma) + feature importance
- [ ] Histogram-baserade splits, validerade mot exakt-greedy
- [ ] `SupervisedExperiment`-integration + grid-search-exempel

### Phase 3 — Transformer-byggstenar
- [ ] `LayerNorm` + `PositionalEncoding` (forward + backward + tester)
- [ ] `SelfAttentionLayer` forward + gradient-check-testad backward (eller autodiff-väg)
- [ ] `MultiHeadAttentionLayer` + `TransformerEncoderBlock`
- [ ] `TransformerClassifier`/`TransformerRegressor` som `ISequenceModel`
- [ ] Benchmark mot BiLSTM/TCN på befintliga sekvenstestfall

### Phase 4 — ONNX
- [ ] Skapa `CSharpNumerics.Onnx`-paketprojekt (protobuf-beroendet isolerat)
- [ ] Export: MLP + skalers/pipeline
- [ ] Round-trip-verifiering mot ONNX Runtime
- [ ] Export: träd-ensembler (`ai.onnx.ml`)
- [ ] Export: Conv1D/LSTM
- [ ] Begränsad import (MLP-vikter)

### Phase 5 — Dimensionsreduktion
- [ ] `IReducer`-interface, `PCA` implementerar det
- [ ] Exakt t-SNE → Barnes–Hut-optimering
- [ ] `Umap` (sträckmål)

---

## Sammanfattning

| Del | Genomförbarhet | Insats | Största risk |
|-----|---------------|--------|--------------|
| **Gradient boosting** | **Hög** | Medel | Histogram-optimeringens korrekthet — mitigeras av exakt-greedy som facit |
| **Transformer** | **Medel-hög** | Medel-hög | Attention-backward för hand — mitigeras av autodiff-spåret |
| **ONNX** | **Medel-hög** | Medel-hög | Operator-mappningsdetaljer (LSTM-layout); separat paket skyddar kärnan |
| **t-SNE/UMAP** | **Medel** | Medel-hög | Lägst prioritet — kan skjutas |

**Rekommendation:** Gradient boosting först — högst efterfrågan, lägst risk, återanvänder trädkoden. Transformer sekvenseras efter autodiff-Phase 4 om möjligt. ONNX kan löpa parallellt (oberoende av allt annat). t-SNE/UMAP tas om tid finns.
