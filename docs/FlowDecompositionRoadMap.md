# Flow Decomposition — CSharpNumerics byggstenar (för RAST.Decomp)

## Mål

Bygga ut `CSharpNumerics` med de **återanvändbara numeriska byggstenarna** (filter, state-estimation, TCN-lager, loss-funktioner, wavelets) som krävs för ett ML-drivet flödesdekompositionssystem — **utan externa Python-bibliotek**.

Själva **systemet `RAST.Decomp` byggs i RAST-repot**, inte här. RAST står för all orkestrering och konsumerar byggstenarna från denna roadmap. Systemet separerar totalflöde $Q_{total}(t)$ i delkomponenter:

$$Q_{total}(t) = Q_{WW}(t) + Q_{SRC}(t) + Q_{FRC}(t) + Q_{Spill}(t)$$

Arkitekturen är hybrid: deterministisk spillextrahering → statistisk WW-baseline → ML-separation av residual (FRC + SRC). Denna roadmap omfattar enbart byggstenarna i Phase 1–5; pipelinen som binder samman dem (Phase 6–7) hör till RAST.

> ### 🧭 Scope — CSharpNumerics vs RAST
>
> **All orkestrering sköts från RAST.** Denna roadmap omfattar enbart de *återanvändbara byggstenarna* (filter, state-estimation, TCN-lager, loss-funktioner, wavelets) som ska implementeras i `CSharpNumerics` som ett generellt numeriskt bibliotek.
>
> Den applikationsspecifika **dekompositionspipelinen** (allt under namespacet `CSharpNumerics.ML.FlowDecomposition` — `FlowDecomposer`, `WWBaselineExtractor`, `DecompositionTrainer`, m.fl.) byggs **INTE i denna kodbas**. Den hör hemma i RAST-repot och konsumerar CSharpNumerics-byggstenarna. Phase 6 och Phase 7 nedan är därför markerade **🚫 RAST** och behålls endast som referens.
>
> **Dokumentation:** varje fas som levererar byggstenar (Phase 1–5) ska även uppdatera modulens README (`Numerics`, `Statistics`, `ML`) — se `### Dokumentation`-listan per fas.

---

## Nulägesanalys — Vad som redan finns

| Komponent | Status | Plats | Relevans |
|-----------|--------|-------|----------|
| `Conv1DLayer` | ✅ | `ML/Sequence/Layers/Conv1DLayer.cs` | Bas för TCN — saknar causal padding & dilation |
| `LSTMLayer` / `BiLSTMLayer` | ✅ | `ML/Sequence/Layers/` | Alternativ sekvensmodell |
| `MaxPool1DLayer` / `GlobalAvgPool1DLayer` | ✅ | `ML/Sequence/Layers/` | Pooling-infrastruktur |
| `SequentialModel` / `ILayer` | ✅ | `ML/NeuralNetwork/` | Lagerstackning, forward/backward |
| `Adam` / `AdamW` / `GradientDescent` | ✅ | `Numerics/Optimization/` | Träningsoptimerare |
| `FourierSeries` | ✅ | `Numerics/SignalProcessing/` | WW diurnal-modellering |
| `TimeSeries` / `Series` | ✅ | `Statistics/Data/` | In/ut-datahantering |
| `TimeSeriesDetrending` | ✅ | `Statistics/TimeSeriesAnalysis/` | Stationaritetsförberedelse |
| `SlidingWindowStatistics` | ✅ | `Statistics/` | Rullande features |
| `PeakFitting` | ✅ | `Statistics/` | Händelsedetektering |
| `NonlinearLeastSquaresFitter` | ✅ | `Statistics/` | Recessionskurvpassning |
| `RobustFitter` (IRLS) | ✅ | `Statistics/` | Outlier-resistent WW-estimering |
| `RollingCrossValidator` | ✅ | `ML/CrossValidators/` | Tidsseriekorsvalidering |
| `MLPRegressor` | ✅ | `ML/` | Baseline-nätverksmodell |
| `RandomForest` | ✅ | `ML/` | Feature importance |
| `KMeans` / `DBSCAN` | ✅ | `ML/` | Händelsekluster |
| `PCA` | ✅ | `ML/` | Dimensionsreduktion |
| `LombScarglePeriodogram` | ✅ | `Statistics/` | Periodicitetsdetektion |
| `StandardScaler` / `MinMaxScaler` | ✅ | `ML/` | Feature-normalisering |
| `ODE Solvers` (RK4, Adaptive) | ✅ | `Numerics/` | Reservoirmodeller |
| `Integration` (Simpson/Trapezoidal) | ✅ | `Numerics/` | Volymberäkning |

### Identifierade luckor

| Saknas | Syfte i RAST.Decomp | Prioritet |
|--------|---------------------|-----------|
| Causal padding i Conv1D | TCN kräver kausal faltning | P1 |
| Dilated convolutions | TCN receptive field utan djupa nät | P1 |
| TCN-block (residual + dilation stack) | Primär FRC/SRC-separator | P1 |
| Dropout-lager | Regularisering under träning | P1 |
| Batch Normalization | Stabilisera djupa nätverksträning | P1 |
| Constrained Loss-funktioner | Icke-negativitet, bevarandelag | P1 |
| Kalman Filter | Spåra WW-baseline med brus | P1 |
| Savitzky-Golay Filter | SRC-utjämning med peak-bevarande | P1 |
| Butterworth IIR-filter | Låg/högpassfiltrering (SRC/FRC) | P1 |
| FIR-filter | Generalfiltrering | P2 |
| Holt-Winters (Triple Exp. Smoothing) | Lättvikts WW-modellering | P2 |
| Wavelet Transform (DWT) | Multi-resolution dekomposition | P2 |
| Attention Mechanism | Långdistansberoenden i flödeshändelser | P3 |
| GRU-cell | Lättviktsalternativ till LSTM | P3 |
| Constrained Optimization (QP) | Strikt $y_i \geq 0$, $\sum y_i = Q$ | P3 |

---

## Phase 1 — Signalbehandling & Filter

Grundläggande DSP-infrastruktur som krävs för pre-processing och feature engineering.

**Namespace:** `CSharpNumerics.Numerics.SignalProcessing`

- [x] `SavitzkyGolayFilter` — polynomutjämning med konfigurerbar ordning och fönsterbredd, bevarar peak-form
- [x] `ButterworthFilter` — IIR lågpass/högpass/bandpass design (ordning N, brytfrekvens), `Apply(double[] signal)` via biquad-kaskad
- [x] `FIRFilter` — finit impulsrespons-filter med godtyckliga koefficienter, `Apply(double[] signal)`
- [x] `FilterDesign` — statisk klass: `DesignLowpass(order, cutoff, sampleRate)`, `DesignHighpass(...)`, `DesignBandpass(...)`
- [x] `ZeroPhaseFiltFilt` — nollfasfiltrering (framåt + bakåt) för offline-analys

### Tests
- [x] Savitzky-Golay med ordning 0 motsvarar glidande medelvärde
- [x] Butterworth lågpass dämpar signal vid 2× brytfrekvens med > 12 dB/oktav
- [x] FiltFilt ger nollfasförskjutning (peak-position bevarad)
- [x] Butterworth högpass + lågpass på samma signal → rekonstruktion inom tolerans

### Dokumentation
Uppdatera **`Numerics/Numerics/Numerics/README.md`** (sektionen *Signal Processing*):
- [x] Lägg till sektion `Filtering` med `SavitzkyGolayFilter`, `ButterworthFilter`, `FIRFilter` (kodexempel: `Apply(double[])`)
- [x] Dokumentera `FilterDesign` (`DesignLowpass/Highpass/Bandpass`) och `ZeroPhaseFiltFilt` med ett before/after-exempel

---

## Phase 2 — Kalman Filter & Exponentiell Utjämning

State-space-modeller för adaptiv WW-baselinespårning.

**Namespace:** `CSharpNumerics.Statistics.StateEstimation`

- [x] `KalmanFilter` — linjärt state-space: Predict(F, Q) → Update(H, R, z), accessors för state `x̂` och kovarians `P`
- [x] `ExtendedKalmanFilter` — icke-linjär variant med Jacobian-funktioner `f(x)`, `h(x)`, `F(x)`, `H(x)`
- [x] `KalmanSmoother` — Rauch-Tung-Striebel bakåtpass för offline-utjämning
- [x] `HoltWintersSmoothing` — triple exponential smoothing: level + trend + säsongskomponent (additiv/multiplikativ), auto-initiering

### Tests
- [x] KalmanFilter på konstant signal med brus → konvergerar till sann nivå
- [x] KalmanFilter spårar linjär ramp med korrekt hastighet
- [x] ExtendedKalmanFilter estimerar sinusformad signal med känd modell
- [x] HoltWinters fångar 24h-säsong i syntetisk WW-serie, prediktionsfel < 5%
- [x] KalmanSmoother ger lägre varians än forward-only Kalman

### Dokumentation
Uppdatera **`Numerics/Numerics/Statistics/README.md`**:
- [x] Ny sektion `State Estimation` (`CSharpNumerics.Statistics.StateEstimation`) med `KalmanFilter`, `ExtendedKalmanFilter`, `KalmanSmoother` — Predict/Update-exempel
- [x] Lägg `HoltWintersSmoothing` under sektionen *Time Series Analysis* (additiv/multiplikativ, säsong)

---

## Phase 3 — TCN-infrastruktur (Causal Conv1D & Dilation)

Utöka befintlig `Conv1DLayer` med stöd för TCN-arkitekturen.

**Namespace:** `CSharpNumerics.ML.Sequence`

- [ ] Utöka `Conv1DLayer` med `PaddingMode.Causal` — left-padding = (kernelSize - 1) × dilation
- [ ] Utöka `Conv1DLayer` med `dilation`-parameter — dilated faltning i forward/backward
- [ ] `DropoutLayer` — stokastiskt nollställer element under träning (rate), passthrough vid inferens
- [ ] `BatchNorm1DLayer` — kanal-vis normalisering för sekvensdata, running mean/var vid inferens
- [ ] `ResidualBlock` — wrapper: Conv1D(causal, dilated) → BatchNorm → ReLU → Dropout → Conv1D → residual add (med 1×1 conv om kanaler ändras)
- [ ] `TCNBlock` — stack av `ResidualBlock` med exponentiellt ökande dilation [1, 2, 4, 8, 16, ...]
- [ ] `TCNRegressor` — TCN → Global pooling → Dense → output; implementerar `ISequenceModel`
- [ ] `TCNClassifier` — TCN → Global pooling → Dense → Softmax; implementerar `ISequenceModel`

### Tests
- [ ] Causal Conv1D: output beror inte på framtida tidssteg (strixt kausal)
- [ ] Dilated conv med dilation=4, kernel=3 har receptive field = 9
- [ ] TCNBlock med 8 lager (dilation 1→128) har receptive field ≥ 512 tidssteg
- [ ] Dropout sätter ~rate-andel av element till noll under träning
- [ ] BatchNorm normaliserar per kanal, running stats konvergerar
- [ ] TCNRegressor tränar och konvergerar på syntetisk sinusvåg
- [ ] ResidualBlock: gradient flödar genom skip-connection (ej vanishing)
- [ ] TCNRegressor: full forward+backward pass med dilation-stack utan krasch

### Dokumentation
Uppdatera **`Numerics/Numerics/ML/README.md`** (sektionerna *Sequence Models* / *Neural Network Building Blocks*):
- [ ] Dokumentera `Conv1DLayer`-utökningarna (`PaddingMode.Causal`, `dilation`)
- [ ] Lägg till `DropoutLayer`, `BatchNorm1DLayer`, `ResidualBlock`, `TCNBlock` i layer-listan
- [ ] Ny sektion `TCN Architecture` med `TCNRegressor` / `TCNClassifier` (hyperparametrar + `SupervisedExperiment`-exempel)

---

## Phase 4 — Constrained Training & Loss-funktioner

Fysik-informerad träning med bevarandelagar.

**Namespace:** `CSharpNumerics.ML.Training`

- [ ] `NonNegativityLoss` — $\lambda \sum \max(-y_i, 0)^2$, differentiabel penalty
- [ ] `ConservationLoss` — $\lambda \|Q_{total} - \sum y_i\|^2$, säkerställer dekomposition summerar
- [ ] `SmoothnessLoss` — $\lambda \|\nabla^2 y\|^2$, penaliserar icke-smooth SRC
- [ ] `CompositeLoss` — kombinerar multipla förlustfunktioner med vikter: `Add(ILoss, weight)`
- [ ] `ConstrainedTrainer` — träningsloop som applicerar `CompositeLoss` och stödjer curriculum learning (successivt ökande constraint-vikter)
- [ ] `SoftmaxConstraintHead` — output-lager som ger ratio ∈ [0,1] med summavillkor via softmax

### Tests
- [ ] NonNegativityLoss = 0 om alla predictions ≥ 0
- [ ] ConservationLoss = 0 om summan av komponenter == total
- [ ] CompositeLoss gradient = summa av viktade delgradienter
- [ ] ConstrainedTrainer med conservation → slutmodell ger summa inom 1% av total
- [ ] SoftmaxConstraintHead producerar output ∈ [0,1] och summa ≤ 1.0

### Dokumentation
Uppdatera **`Numerics/Numerics/ML/README.md`**:
- [ ] Ny sektion `Constrained Training` (`CSharpNumerics.ML.Training`) med `NonNegativityLoss`, `ConservationLoss`, `SmoothnessLoss`, `CompositeLoss`
- [ ] Dokumentera `ConstrainedTrainer` (curriculum learning) och `SoftmaxConstraintHead`

---

## Phase 5 — Wavelet Transform

Multi-resolution tidsfrekvens-dekomposition.

**Namespace:** `CSharpNumerics.Numerics.SignalProcessing.Wavelets`

- [ ] `DiscreteWaveletTransform` — DWT med konfigurerbara wavelet-familjer, N-nivå dekomposition → approximation + detail-koefficienter
- [ ] `WaveletFamily` — enum/klass: `Haar`, `Daubechies4`, `Daubechies8`, `Symlet4` med filterkoefficienter
- [ ] `InverseWaveletTransform` — rekonstruktion från koefficienter
- [ ] `WaveletDenoising` — soft/hard thresholding på detail-koefficienter (VisuShrink/BayesShrink)
- [ ] `MaximalOverlapDWT` (MODWT) — shift-invariant variant för tidsserieanalys

### Tests
- [ ] DWT → IDWT round-trip: rekonstruktionsfel < 1e-12
- [ ] Haar-wavelet på stegfunktion ger förväntat detail-mönster
- [ ] Daubechies-4 på sinusvåg: energi koncentrerad i rätt nivå
- [ ] WaveletDenoising minskar RMSE på brusig signal med > 50%
- [ ] MODWT-koefficienter är tidsinvariant (shiftad input → shiftad output)

### Dokumentation
Uppdatera **`Numerics/Numerics/Numerics/README.md`** (sektionen *Signal Processing*):
- [ ] Ny undersektion `Wavelets` (`CSharpNumerics.Numerics.SignalProcessing.Wavelets`) med `DiscreteWaveletTransform`, `WaveletFamily`, `InverseWaveletTransform`
- [ ] Dokumentera `WaveletDenoising` (soft/hard thresholding) och `MaximalOverlapDWT` (MODWT) — DWT→IDWT round-trip-exempel

---

## Phase 6 — Dekompositionspipeline & Modell &nbsp; 🚫 RAST

Sammanfoga alla delar till en komplett dekompositionsmodell.

> 🚫 **Byggs i RAST, inte i CSharpNumerics.** Hela detta namespace (`CSharpNumerics.ML.FlowDecomposition`) är applikationsspecifik orkestrering och implementeras i RAST-repot, som konsumerar byggstenarna från Phase 1–5. Listan nedan behålls endast som referens för RAST-arbetet. Inga README-tasks för CSharpNumerics.

**Namespace:** `RAST.Decomp` (tidigare planerat som `CSharpNumerics.ML.FlowDecomposition`)

- 🚫 `FlowDecomposer` — orkestreringsklass: konfigurerar pipeline (spill → WW → residual → ML → output)
- 🚫 `WWBaselineExtractor` — Fourier-harmonisk modell (N harmonics) + Kalman-spårning av drift
- 🚫 `ResidualComputer` — beräknar $Q_{res} = Q_{total} - Q_{spill} - Q_{WW}$
- 🚫 `TCNDecompositionModel` — tränad TCN som ger FRC/SRC-ratio givet residual + features
- 🚫 `DecompositionResult` — struct: WW, FRC, SRC, Spill som `double[]`, med kvalitetsmetrik
- 🚫 `DecompositionMetrics` — konserveringsfel, icke-negativitetscount, SRC-smoothness, recessionpassning
- 🚫 `RecessionAnalyzer` — identifiera och anpassa exponentiell recession $Q_0 e^{-t/\tau}$ på FRC/SRC
- 🚫 `FeatureEngineering` — extrahera temporal features: $dQ/dt$, rullande stats, hour_sin/cos, DOW-encoding

### Tests (i RAST)
- 🚫 FlowDecomposer på syntetisk data (känt WW + FRC + SRC) → rekonstruktion inom 2%
- 🚫 Conservation: $|Q_{total} - (WW + FRC + SRC + Spill)| < 0.01$ vid varje tidssteg
- 🚫 Alla komponenter icke-negativa
- 🚫 FRC-recession matchar exponentiell profil med $R^2 > 0.95$
- 🚫 SRC har bounded second derivative (smoothness-villkor)
- 🚫 DecompositionMetrics korrekt beräknade på kända testvektorer

---

## Phase 7 — Träning & Validering &nbsp; 🚫 RAST

Fullständig träningspipeline med DHI Sewdec som svag supervision.

> 🚫 **Byggs i RAST, inte i CSharpNumerics.** Träningsorkestreringen lever i RAST-repot. `HyperparameterSearch` återanvänder dock CSharpNumerics-primitiv (`PipelineGrid`, `RollingCrossValidator`) som redan finns. Listan nedan behålls endast som referens.

**Namespace:** `RAST.Decomp` (tidigare planerat som `CSharpNumerics.ML.FlowDecomposition`)

- 🚫 `DecompositionTrainer` — curriculum learning: Fas A (Sewdec-labels) → Fas B (+ fysik-constraints) → Fas C (self-supervised refinement)
- 🚫 `TrainingDataLoader` — läser parad data (totalflöde + Sewdec-output) från CSV, tidsjusterar, normaliserar
- 🚫 `CurriculumSchedule` — schema för constraint-vikter $\lambda_i(epoch)$: lineär/exponentiell upptrappning
- 🚫 `ModelSerializer` — spara/ladda TCN-vikter till binärfil (egna format, inga externa beroenden)
- 🚫 `HyperparameterSearch` — integrera med `PipelineGrid` + `RollingCrossValidator` för optimal TCN-arkitektur
- 🚫 `AblationRunner` — jämför TCN vs LSTM vs matematisk baseline vs Sewdec på samma data

### Tests (i RAST)
- 🚫 ModelSerializer save → load round-trip: identisk inferens
- 🚫 CurriculumSchedule producerar monotont ökande vikter
- 🚫 Tränad modell överpresterar naive baseline (all residual → FRC) på hållet-out-data
- 🚫 RollingCrossValidator undviker dataläckage (test alltid efter train i tid)

---

## Sammanfattning — Beroendegraf

```
╔════════════════════════════════╗
║  RAST-repot (orkestrering)      ║
║   Phase 7 (Träning)            ║
║       ↓                         ║
║   Phase 6 (Pipeline)           ║
╚════════════════════════════════╝
            ↑  konsumerar byggstenar
─────────────────────────────────────
  CSharpNumerics (byggstenar — denna roadmap)
    Phase 4        Phase 5
    (Constraints)  (Wavelets)
        ↓
    Phase 3 (TCN)
        ↓
    Phase 1 ← Phase 2
    (DSP)    (Kalman)
```

Phase 1 och 2 kan köras parallellt. Phase 3 kräver Phase 1 (filterinfrastruktur). Phase 4 och 5 kan köras parallellt efter Phase 3. **Phase 6 och 7 byggs i RAST** och kräver att Phase 1–5 är levererade som byggstenar i CSharpNumerics.

---

## Namespace-placering

| Komponent | Namespace | Motivering |
|-----------|-----------|-----------|
| SavitzkyGolayFilter, Butterworth, FIR | `Numerics.SignalProcessing` | Ren signalbehandling |
| FilterDesign, ZeroPhaseFiltFilt | `Numerics.SignalProcessing` | DSP-verktyg |
| KalmanFilter, EKF, Smoother | `Statistics.StateEstimation` | Stokastisk estimering |
| HoltWintersSmoothing | `Statistics.TimeSeriesAnalysis` | Tidsseriestatistik |
| DropoutLayer, BatchNorm1DLayer | `ML.Sequence.Layers` | Nätverkslager |
| ResidualBlock, TCNBlock | `ML.Sequence.Layers` | TCN-arkitekturblock |
| TCNRegressor, TCNClassifier | `ML.Sequence.Models` | Sekvensmodeller |
| NonNegativityLoss, CompositeLoss | `ML.Training` | Träningsinfrastruktur |
| ConstrainedTrainer | `ML.Training` | Träningsloop |
| DiscreteWaveletTransform, MODWT | `Numerics.SignalProcessing.Wavelets` | Wavelet-matematik |
| ~~FlowDecomposer, WWBaseline, etc.~~ | 🚫 `RAST.Decomp` | Applikationsspecifik pipeline — **byggs i RAST** |
| ~~DecompositionMetrics~~ | 🚫 `RAST.Decomp` | Kvalitetsmätning — **byggs i RAST** |
