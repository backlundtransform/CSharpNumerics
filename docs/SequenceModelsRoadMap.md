# Sequence Models (1D-CNN & Bi-LSTM) — Architecture Feasibility Study

## Mål

Utreda möjligheten att bygga ut `SupervisedExperiment`-pipelinen med **1D-CNN** och **Bi-LSTM**-modeller, med det långsiktiga målet att köra ML på **ljuskurvor för exoplanet-transitdata**.

---

## Nulägesanalys — Befintlig arkitektur

### Vad som finns idag

| Komponent | Status | Plats | Kommentar |
|-----------|--------|-------|-----------|
| `NeuralNetwork` | ✓ | `ML/NeuralNetwork/NeuralNetwork.cs` | Fullt feedforward-nät med backprop, mini-batch, early stopping |
| `Matrix` | ✓ | `Numerics/Objects/Matrix.cs` | Transpose, Slice, Multiply, + / - / * |
| `VectorN` | ✓ | `Numerics/Objects/VectorN.cs` | Dot, Hadamard (⊙), Outer, Norm, Slice |
| `IOptimizer` | ✓ | `Numerics/Optimization/` | Adam (med AdamW), GradientDescent (med Nesterov/Momentum) |
| `EarlyStopping` | ✓ | `Numerics/Optimization/Strategies/` | Patience-baserad tidig avslutning |
| `LearningRateSchedule` | ✓ | `Numerics/Optimization/Strategies/` | Decay-scheman |
| `ActivationType` | ✓ | `ML/Enums/` | ReLU, Sigmoid, Tanh, Linear |
| `IModel`/`IClassificationModel`/`IRegressionModel` | ✓ | `ML/Models/Interfaces/` | Fit(X,y) → Predict(X) → Clone() |
| `IHasHyperparameters` | ✓ | `ML/Models/Interfaces/` | Grid-search-kompatibilitet |
| `SupervisedExperiment` | ✓ | `ML/Experiment/` | Fluent API, grid-search, multipla cross-validators |
| `Pipeline` | ✓ | `ML/Pipeline.cs` | Scaler → Selector → Reducer → Model |
| `TimeSeries` | ✓ | `Statistics/Data/TimeSeries.cs` | DateTime[], double[][], ToMatrix() |
| `RollingCrossValidator` | ✓ | `ML/CrossValidators/` | Tidsserie-medveten korsvalidering |
| LSTM/GRU/RNN-lager | ✗ | — | Finns ej |
| Faltningslager (Conv) | ✗ | — | Finns ej |
| Layer-abstraktion | ✗ | — | NeuralNetwork är monolitisk; inga ILayer-interface |

### Nyckelidentifierade begränsningar

1. **Monolitisk NeuralNetwork** — Alla lager är implicita i `List<Matrix>` + `List<VectorN>`. Nätverket stödjer enbart fully-connected (Dense) lager. Det finns inget `ILayer`-interface för att plugga in nya lagertyper.

2. **Ingen sekvenshantering i Forward/Backward** — Forward tar en enda `VectorN` (inte en sekvens av tidssteg). Backprop saknar stöd för BPTT (Backpropagation Through Time).

3. **IModel-kontraktet är flat** — `Fit(Matrix X, VectorN y)` där varje rad i X är ett oberoende sample. Sekvensmodeller behöver 3D-indata: `(samples × timesteps × features)`.

---

## Alternativ 1 — 1D-CNN (Convolutional Neural Network)

### Vad det gör

En 1D-CNN applicerar lärda kärnor (kernels) som glider över en tidsserie (eller annan 1D-signal) och producerar ett featuremap per kanal. Typisk arkitektur:

```
Input [timesteps × features]
  ↓
Conv1D(filters=32, kernel_size=5, activation=ReLU)
  ↓
MaxPool1D(pool_size=2)
  ↓
Conv1D(filters=64, kernel_size=3, activation=ReLU)
  ↓
GlobalAveragePool / Flatten
  ↓
Dense → Output (classification/regression)
```

### Varför det är relevant för ljuskurvor

Exoplanet-transiter skapar **lokala dips** i ljuskurvan — korta, formberoende mönster inbäddade i en längre tidsserie. 1D-CNN:er excellerar på att detektera sådana **lokaliserade mönster** oavsett var i sekvensen de befinner sig (translationsinvarians).

### Vad som krävs att bygga

| Komponent | Beskrivning | Beroenden |
|-----------|-------------|-----------|
| **`ILayer` interface** | `VectorN[] Forward(VectorN[] input)`, `VectorN[] Backward(VectorN[] gradOutput)`, `(Matrix[], VectorN[]) GetGradients()`, `ApplyGradients(IOptimizer)` | Nytt abstrakt lager |
| **`Conv1DLayer`** | Per kanal: slid kernel-vektor över input, producera featuremap. Stöd: `filters`, `kernelSize`, `stride`, `padding` (same/valid). Vikter: `Matrix[filters × (kernelSize*inputChannels)]` + bias | `VectorN.Hadamard`, `VectorN.Dot`, `Matrix.Slice`, manuell faltningsloop |
| **`MaxPool1DLayer`** | Reducera sekvens via max-pooling. Spara argmax-index för backprop | Enkel indexering |
| **`GlobalAvgPool1DLayer`** | medelvärde per kanal → platt VectorN | `VectorN`-aggregering |
| **`FlattenLayer`** | Konkaternera multi-kanal → enda VectorN | — |
| **`DenseLayer`** | Wrapper runt befintlig NeuralNetwork för enstaka Dense-lager. Alternativt: ny klass med `Matrix * VectorN + bias` | Redan tillgängligt via `NeuralNetwork` |
| **`SequentialModel`** | Ny modellkompositionslogik: `Forward` propagerar genom `ILayer[]`, `Backward` i omvänd ordning, `ApplyGradients` itererar alla lager | Nytt — ersätter monolitisk `NeuralNetwork` för denna modelltyp |
| **`CNN1DClassifier` / `CNN1DRegressor`** | Implementerar `IClassificationModel` / `IRegressionModel` + `IHasHyperparameters`. Sätts ihop av `Conv1DLayer → Pool → Dense → Output`. Gör modellen kompatibel med `SupervisedExperiment` och grid-search | `SequentialModel`, alla lager ovan |

### Matematiska primitiver som redan finns

| Operation | Finns i | Användning |
|-----------|---------|------------|
| Elementvis multiplikation (⊙) | `VectorN.Hadamard()` | Gate-produkter, gradientskala i backprop |
| Matris × vektor | `Matrix * VectorN` | Dense-lager forward |
| Outer product | `VectorN.Outer()` | Gradientberäkning (dW = a ⊗ δ) |
| Slice | `Matrix.Slice()`, `VectorN.Slice()` | Fönsterextraktion för faltning |
| Sigmoid/Tanh/ReLU | `NeuralNetwork.Activate()` (privat) | Behöver göras publik eller dupliceras i `ILayer`-impl |

### Vad som saknas

| Operation | Beskrivning | Implementationsinsats |
|-----------|-------------|----------------------|
| **1D-faltning** | `∑ kernel[j] * input[i+j]` för varje stegposition | Enkel loop med `VectorN.Dot` per fönster — **låg insats** |
| **Padding** | Zero-padding av inputsekvens | `VectorN`-konkatenering med nollor — **låg insats** |
| **3D-datarepresentation** | `(samples × timesteps × features)` | `Matrix[]` (lista av matriser, en per sample) eller ny `Tensor3D`-klass — **medel insats** |
| **Aktivering som fristående klass** | Extern tillgång till ReLU, Sigmoid etc. med derivata | Refaktorera ut från `NeuralNetwork` — **låg insats** |

### Genomförbarhet

**HÖG.** 1D-faltning är i grunden en dot-product med sliding window — alla nödvändiga primitiver (dot, Hadamard, slice, matmul) finns redan. Huvudinsatsen är att införa ett `ILayer`-abstraktion och bygga den modulära lagerkompositionen (`SequentialModel`). Faltningsmatematiken i sig är trivial att implementera.

---

## Alternativ 2 — Bi-LSTM (Bidirectional Long Short-Term Memory)

### Vad det gör

En LSTM processar sekventiell data ett tidssteg i taget och upprätthåller ett cell-state $C_t$ och ett dolt tillstånd $h_t$ via fyra gates:

$$f_t = \sigma(W_f \cdot [h_{t-1}, x_t] + b_f) \quad \text{(forget gate)}$$
$$i_t = \sigma(W_i \cdot [h_{t-1}, x_t] + b_i) \quad \text{(input gate)}$$
$$\tilde{C}_t = \tanh(W_C \cdot [h_{t-1}, x_t] + b_C) \quad \text{(candidate)}$$
$$o_t = \sigma(W_o \cdot [h_{t-1}, x_t] + b_o) \quad \text{(output gate)}$$
$$C_t = f_t \odot C_{t-1} + i_t \odot \tilde{C}_t$$
$$h_t = o_t \odot \tanh(C_t)$$

En **Bi-LSTM** kör två LSTM:er parallellt — en framåt och en bakåt i tid — och konkatanerar deras output.

### Varför det är relevant för ljuskurvor

Transit-dips omges av kontext *på båda sidor* i tidsserien. En Bi-LSTM fångar temporala beroenden i båda tidsdirektionerna, vilket gör den kraftfull för att skilja äkta transiter från brusartefakter och variabla stjärnor.

### Vad som krävs att bygga

| Komponent | Beskrivning | Beroenden |
|-----------|-------------|-----------|
| **`ILayer` interface** | (Samma som för 1D-CNN — delas) | Nytt |
| **`LSTMCell`** | En enda LSTM-cell: Forward tar $(x_t, h_{t-1}, C_{t-1})$ → $(h_t, C_t)$. Backward beräknar gradienter via BPTT-steg. Vikter: $W_f, W_i, W_C, W_o$ (vardera `Matrix[(input+hidden) × hidden]`) + 4 bias-vektorer | `VectorN.Hadamard`, `Matrix * VectorN`, Sigmoid/Tanh-aktivering |
| **`LSTMLayer`** | Kör `LSTMCell` över hela sekvensen (T tidssteg). Cachar alla $(h_t, C_t)$ för BPTT. Truncated BPTT per konfiguration | `LSTMCell`, sekvensdatastruktur |
| **`BiLSTMLayer`** | Två `LSTMLayer` — en framåt, en bakåt. Konkatenerar `h_fwd[t] ∥ h_bwd[t]` per tidssteg → output-dimension = `2 × hidden` | `LSTMLayer` × 2 |
| **BPTT (Backpropagation Through Time)** | Specifik backprop-algoritm för recurrenta nät. Per tidssteg: propagera gradienterna bakåt genom gates med Hadamard-produkter. Hantera gradient clipping | `VectorN.Hadamard`, `Matrix.Transpose`, gradient-clipping-logik |
| **Gradient Clipping** | Klipp gradient-norm för att hantera exploding gradients (vanligt i RNN). `max_norm`-parameter | `VectorN.Norm()` — finns redan |
| **`BiLSTMClassifier` / `BiLSTMRegressor`** | Implementerar `IModel` + `IHasHyperparameters`. Komposition: BiLSTMLayer → Dense → Output. Kompatibel med `SupervisedExperiment` | `SequentialModel`, `BiLSTMLayer`, `DenseLayer` |

### Matematiska primitiver som redan finns

| Operation | Finns i | Användning i LSTM |
|-----------|---------|-------------------|
| Sigmoid $\sigma(x) = 1/(1+e^{-x})$ | `ActivationType.Sigmoid` | Forget, input, output gates |
| Tanh | `ActivationType.Tanh` | Kandidat-cell, output-gate |
| Hadamard (⊙) | `VectorN.Hadamard()` | $f_t \odot C_{t-1}$, $i_t \odot \tilde{C}_t$, $o_t \odot \tanh(C_t)$ |
| Matris × vektor | `Matrix * VectorN` | Gate-beräkningar $W \cdot [h, x]$ |
| Vektor-konkatenering | ✗ **finns ej** | $[h_{t-1}, x_t]$ behöver concat-operation |
| L2-norm | `VectorN.Norm()` | Gradient clipping |

### Vad som saknas

| Operation | Beskrivning | Implementationsinsats |
|-----------|-------------|----------------------|
| **VectorN.Concat()** | Konkatenera två vektorer $[h_{t-1} \| x_t]$ → ny VectorN | **Trivial** — ny metod, ~5 rader |
| **BPTT** | Backpropagation Through Time med gate-derivator | **Medel-hög insats** — korrekt gradient-flöde genom alla 4 gates × T tidssteg |
| **Gradient clipping** | `clip(‖g‖, max_norm)` → skala ned gradient om normen överskrider | **Låg insats** — `VectorN.Norm()` finns redan |
| **3D-datarepresentation** | `(samples × timesteps × features)` | Samma som för 1D-CNN |
| **Sekventiell batch-datahantering** | Packa/padda sekvenser till samma längd per batch | **Medel insats** |

### Genomförbarhet

**MEDEL-HÖG.** Alla gate-operationer (Hadamard, sigmoid, tanh, matmul) finns redan. Den största utmaningen är att korrekt implementera **BPTT** med alla fyra LSTM-gates — detta är matematiskt komplicerat men väl dokumenterat. Bi-LSTM adderar liten extra komplexitet utöver en vanlig LSTM. Gradient clipping behövs men är trivialt givet befintlig `Norm()`.

---

## Gemensam infrastruktur — Vad båda modellerna kräver

### 1. Layer-abstraktion (`ILayer` interface)

```csharp
namespace CSharpNumerics.ML.NeuralNetwork.Layers;

public interface ILayer
{
    // Forward pass: input activations → output activations
    // Formen beror på lagertypen (1D, sekventiell, etc.)
    VectorN[] Forward(VectorN[] input, bool training = true);
    
    // Backward pass: output-gradienter → input-gradienter
    VectorN[] Backward(VectorN[] gradOutput);
    
    // Tillämpa ackumulerade gradienter
    void ApplyGradients(IOptimizer weightOptimizer, IOptimizer biasOptimizer, int batchSize);
    
    // Antal lärda parametrar (för diagnostik/loggning)
    int ParameterCount { get; }
}
```

### 2. Sekventiell modellkomposition (`SequentialModel`)

En ny klass som ersätter behovet av monolitisk `NeuralNetwork` för komplexa arkitekturer:

```csharp
public class SequentialModel
{
    private readonly ILayer[] _layers;
    
    public VectorN Forward(VectorN[] input) { /* propagera framåt */ }
    public void Backward(VectorN[] lossGrad) { /* propagera bakåt */ }
    public void ApplyGradients(IOptimizer wOpt, IOptimizer bOpt, int batchSize) { /* alla lager */ }
}
```

### 3. 3D-datarepresentation

`IModel.Fit(Matrix X, VectorN y)` förväntar sig 2D-data. Sekvensmodeller behöver 3D:

**Alternativ A — Flatten-konvention:** Omforma `(samples × timesteps × features)` → `Matrix[samples × (timesteps*features)]` och låt modellen internt reshape:a. Fördel: inga ändringar i `IModel`-interface. Nackdel: implicit kontrakt.

**Alternativ B — Nytt `ISequenceModel`-interface:**
```csharp
public interface ISequenceModel : IModel
{
    int TimeSteps { get; set; }
    int Features { get; set; }
    // Fit/Predict ärver från IModel men tolkar Matrix-rader som flattened sekvenser
}
```
Fördel: explicit kontrakt, kompatibelt med `SupervisedExperiment` (som bara kräver `IModel`).

**Rekommendation:** Alternativ B — inför `ISequenceModel` som sub-interface. `SupervisedExperiment` behöver inga ändringar; alla befintliga cross-validators och grid-search-mekanismer fungerar oförändrat.

### 4. Aktiveringar som fristående funktioner

De nuvarande aktiveringarna (ReLU, Sigmoid, Tanh, Linear) och deras derivator är privata metoder i `NeuralNetwork`. Dessa behöver bli tillgängliga separat:

```csharp
namespace CSharpNumerics.ML.NeuralNetwork;

public static class Activations
{
    public static VectorN Apply(VectorN v, ActivationType type) { ... }
    public static VectorN Derivative(VectorN v, ActivationType type) { ... }
}
```

---

## Exoplanet-transit ljuskurvor — Domänkontext

### Dataformat

En typisk exoplanet-transitljuskurva:

```
tid (JD)    |  flux (normaliserat)  |  flux_err
2459000.00  |  1.0003               |  0.0002
2459000.01  |  0.9998               |  0.0003
...         |  (transit-dip ~0.01)  |  ...
```

- **Tidssteg:** 100–10 000 per observation (Kepler: ~70 000 per kvartal)
- **Features per tidssteg:** 1–3 (flux, flux_error, eventuellt centroid)
- **Labels:** Binär (transit / ingen transit) eller regression (transitdjup, period)

### Varför ML?

- **BLS (Box-fitting Least Squares)** är standardmetoden men missar grunda transiter och är känslig för stellar variability
- **1D-CNN** kan lära sig transitformens lokala signatur direkt från data
- **Bi-LSTM** kan modellera den temporala kontexten runt en transit
- **Kombination (CNN-LSTM)** — faltning för lokal featureextraktion → LSTM för sekventiell kontext — är state-of-the-art i litteraturen

### Referensarkitekturer (litteraturen)

1. **Shallue & Vanderburg (2018)** — AstroNet: 1D-CNN på foldade ljuskurvor (Kepler), >95% precision
2. **Ansdell et al. (2018)** — Random Forest + CNN ensemble
3. **Dattilo et al. (2019)** — AstroNet-K2: anpassad CNN för K2-data
4. **Osborn et al. (2020)** — CNN + Bi-GRU (GRU är en förenklad LSTM)

---

## Implementationsplan — Faser

### Phase 1 — Infrastruktur (gemensam)
- [x] Designa och implementera `ILayer` interface i `ML/NeuralNetwork/Layers/`
- [x] Implementera `DenseLayer` (wrapper runt befintlig Dense-logik)
- [x] Implementera fristående `Activations`-klass (refaktorera från NeuralNetwork)
- [x] Implementera `SequentialModel` med Forward/Backward/ApplyGradients
- [x] Implementera `VectorN.Concat()` i `Numerics/Objects/VectorN.cs`
- [x] Designa `ISequenceModel` sub-interface
- [x] Enhetstester för alla nya infrastrukturkomponenter

### Phase 2 — 1D-CNN
- [x] Implementera `Conv1DLayer` (forward + backward)
- [x] Implementera `MaxPool1DLayer` (forward + backward)
- [x] Implementera `GlobalAvgPool1DLayer`
- [x] Implementera `FlattenLayer`
- [x] Implementera `CNN1DClassifier` (IClassificationModel + IHasHyperparameters)
- [x] Implementera `CNN1DRegressor` (IRegressionModel + IHasHyperparameters)
- [x] Enhetstester: faltning, pooling, end-to-end-träning
- [x] Verifiera kompatibilitet med `SupervisedExperiment` grid-search

### Phase 3 — LSTM
- [x] Implementera `LSTMCell` (forward + backward med gate-derivator)
- [x] Implementera `LSTMLayer` (sekventiell forward + BPTT)
- [x] Implementera gradient clipping
- [x] Enhetstester: cell-level, gradient-flöde, enkel sekvensuppgift

### Phase 4 — Bi-LSTM
- [x] Implementera `BiLSTMLayer` (framåt + bakåt LSTM, konkatenering)
- [x] Implementera `BiLSTMClassifier` / `BiLSTMRegressor`
- [x] Enhetstester: bidirektionell korrekthet, end-to-end-träning
- [x] Verifiera kompatibilitet med `SupervisedExperiment`

### Phase 5 — Hybrid CNN-LSTM (valfri framtida fas)
- [ ] Implementera `CNN_BiLSTMClassifier` (Conv1D → BiLSTM → Dense)
- [ ] Benchmark mot ren CNN och ren Bi-LSTM

### Phase 6 — Integration & Dokumentation
- [x] Uppdatera `ML/README.md` med nya modelltyper
- [x] Integrera med befintlig `TimeSeries`-klass (utility: ljuskurva → flattened Matrix)
- [x] Exempelkod: exoplanet-transit klassificering

---

## Sammanfattning

| Modell | Genomförbarhet | Insats | Beroenden som saknas | Största risk |
|--------|---------------|--------|---------------------|--------------|
| **1D-CNN** | **Hög** | Medel | `ILayer`, `SequentialModel`, sliding-window conv | Liten — matematiken är trivial |
| **Bi-LSTM** | **Medel-Hög** | Medel-Hög | `ILayer`, `SequentialModel`, BPTT, `VectorN.Concat` | BPTT-korrekthet, vanishing/exploding gradients |
| **CNN-LSTM** | **Hög** (om ovan finns) | Låg (komposition) | Allt ovan | Arkitekturkoppling |

**Rekommendation:** Börja med **Phase 1 (infrastruktur)** som är gemensam för båda modellerna, sedan **Phase 2 (1D-CNN)** som är enklast att validera. LSTM (Phase 3–4) bygger vidare på samma infrastruktur och kan implementeras sekventiellt. Allt ryms inom befintlig `IModel`-arkitektur — inga breaking changes krävs.

Den befintliga `SupervisedExperiment`-pipelinen med grid-search och cross-validators kommer att fungera direkt, förutsatt att CNN/Bi-LSTM-modellerna implementerar `IModel` + `IHasHyperparameters`. Inget nytt experiment-ramverk behövs.
