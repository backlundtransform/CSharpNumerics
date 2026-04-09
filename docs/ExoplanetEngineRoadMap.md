# Exoplanet Engine Roadmap

## Mål

Bygga `CSharpNumerics.Engines.Exoplanet` — en analysmotor för exoplanet-transitdetektion som kombinerar befintliga moduler:

- **`CSharpNumerics.Physics.Optics`** — limbmörkning, transitgeometri
- **`CSharpNumerics.Statistics.Fitting`** — transitmodell-anpassning (Mandel-Agol, trapetsmodell)
- **`CSharpNumerics.Statistics.TimeSeriesAnalysis`** — BLS, Lomb-Scargle, avtrending, fasveckning
- **`CSharpNumerics.ML.Sequence`** — CNN1D/LSTM/BiLSTM-klassificering av ljuskurvefönster

**Slutmål:** Träna på märkt (labeled) transitdata → exportera tränad modell → prediktioner via webbtjänst.

API-integrationer mot NASA Exoplanet Archive, MAST/TESS, etc. sker **utanför** CSharpNumerics, men interna datamodeller ska matcha standardformaten (Kepler/TESS FITS-kolumner, KOI-tabeller).

---

## Nulägesanalys — Befintliga resurser

| Komponent | Status | Namespace / Klass |
|-----------|--------|-------------------|
| BLS-periodanalys | ✅ | `Statistics.TimeSeriesAnalysis.BoxFittingLeastSquares` |
| Lomb-Scargle-periodogram | ✅ | `Statistics.TimeSeriesAnalysis.LombScarglePeriodogram` |
| Fasveckning & binning | ✅ | `Statistics.TimeSeriesAnalysis.PhaseFolding` |
| Avtrending (polynom, SG, median) | ✅ | `Statistics.TimeSeriesAnalysis.TimeSeriesDetrending` |
| Toppdetektering | ✅ | `Statistics.TimeSeriesAnalysis.PeakFitting` |
| Icke-linjär minsta-kvadrat (LM) | ✅ | `Statistics.Fitting.NonlinearLeastSquaresFitter` |
| Robust anpassning (IRLS) | ✅ | `Statistics.Fitting.RobustFitter` |
| Bootstrap-konfidensintervall | ✅ | `Statistics.Fitting.ParameterEstimation` |
| Modellval (AIC/BIC) | ✅ | `Statistics.Fitting.GoodnessOfFit` |
| CNN1DClassifier | ✅ | `ML.Sequence.Models.Classification.CNN1DClassifier` |
| LSTMClassifier | ✅ | `ML.Sequence.Models.Classification.LSTMClassifier` |
| BiLSTMClassifier | ✅ | `ML.Sequence.Models.Classification.BiLSTMClassifier` |
| SupervisedExperiment (grid search, CV) | ✅ | `ML.Experiment.SupervisedExperiment` |
| SequenceDataHelper.CreateWindows | ✅ | `ML.Sequence.SequenceDataHelper` |
| TimeSeries-datastruktur | ✅ | `Statistics.Data.TimeSeries` |
| Ray-tracing, optiska ytor, sensor | ✅ | `Physics.Optics.*` |
| Astronomi (avstånd-konverteringar) | ⚠️ Minimal | `Physics.Astro.AstronomyExtensions` |
| Engine-infrastruktur (ISimulationEngine, EventBus, Clock) | ✅ | `Engines.Common.*` |
| Exoplanet-specifik transitfysik | ❌ | — |
| Ljuskurve-datamodell (TESS/Kepler-kompatibel) | ❌ | — |
| Feature-extraktion för transiter | ❌ | — |
| Inference-pipeline (träna → exportera → prediktera) | ❌ | — |

---

## Externt dataformat — Referens

Integrationer sker utanför CSharpNumerics, men interna typer ska mappas **1:1** mot standardformat. Nedan sammanfattas de viktigaste kolumnerna.

### Kepler/TESS ljuskurva (tidsserieformat)

| Kolumn | Typ | Enhet | Beskrivning |
|--------|-----|-------|-------------|
| `TIME` | `double` | BJD − 2457000 (TESS) / BJD − 2454833 (Kepler) | Barycentrisk Julian Date |
| `SAP_FLUX` | `double` | e⁻/s | Simple Aperture Photometry |
| `PDCSAP_FLUX` | `double` | e⁻/s | Pre-search Data Conditioned flux (systematikrensad) |
| `SAP_FLUX_ERR` | `double` | e⁻/s | Osäkerhet SAP |
| `PDCSAP_FLUX_ERR` | `double` | e⁻/s | Osäkerhet PDCSAP |
| `QUALITY` | `int` | flaggmask | Datakvalitetsflaggor (0 = ok) |
| `CADENCENO` | `int` | — | Sekvensnummer |

### KOI / TOI transit-parametrar

| Parameter | Kolumn (Kepler) | Kolumn (TESS TOI) | Enhet |
|-----------|-----------------|---------------------|-------|
| Orbital period | `koi_period` | `pl_orbper` | dagar |
| Transit epoch | `koi_time0bk` | `pl_tranmid` | BJD |
| Transit depth | `koi_depth` | `pl_trandep` | ppm |
| Transit duration | `koi_duration` | `pl_trandur` | timmar |
| Planet/star radius-ratio | `koi_ror` | `pl_ratror` | — |
| Impact parameter | `koi_impact` | `pl_imppar` | — |
| Disposition | `koi_disposition` | `tfopwg_disp` | CANDIDATE / CONFIRMED / FALSE POSITIVE |
| Disposition score | `koi_score` | — | 0–1 |
| Stellar Teff | `koi_steff` | `st_teff` | K |
| Stellar radius | `koi_srad` | `st_rad` | R☉ |
| Stellar logg | `koi_slogg` | `st_logg` | cgs |

---

## Phase 1 — Datamodell & ljuskurvetyper

Definiera de centrala datatyperna i `Engines/Exoplanet/` som matchar externt format. Inga externa beroenden — bara POCOs och konverteringslogik mot befintligt `TimeSeries`.

- [ ] Skapa `Engines/Exoplanet/` mappstruktur med `Data/`, `Pipeline/`, `Features/`, `Interfaces/`, `Enums/`
- [ ] `LightCurve` — klass med `double[] Time`, `double[] Flux`, `double[] FluxError`, `int[] QualityFlags`, `LightCurveMetadata Metadata`; factory-metod `FromTimeSeries(TimeSeries ts, string timeCol, string fluxCol, string errCol)`
- [ ] `LightCurveMetadata` — `string TargetId`, `string Mission` (Kepler/TESS/K2), `double TimeOffset` (BJD-offset), `CadenceType Cadence` (Short/Long)
- [ ] `TransitParameters` — record/klass: `double Period`, `double Epoch`, `double Depth`, `double Duration`, `double RadiusRatio`, `double ImpactParameter`, `double IngressDuration`
- [ ] `TransitCandidate` — `TransitParameters Parameters`, `double Score`, `TransitDisposition Disposition`, `LightCurve PhaseFoldedCurve`, `TransitFeatureSet Features`
- [ ] `StellarProperties` — `double EffectiveTemp`, `double Radius`, `double Mass`, `double SurfaceGravity`, `double Metallicity`, `SpectralType Type`
- [ ] `TransitDisposition` (enum) — `Unknown`, `Candidate`, `Confirmed`, `FalsePositive`, `AstrophysicalFalsePositive`
- [ ] `CadenceType` (enum) — `Short`, `Long`, `Fast` (TESS 20s)
- [ ] `LightCurveSanitizer` — statisk klass: `RemoveBadQuality(LightCurve, qualityMask)`, `RemoveOutliers(LightCurve, sigmaThreshold)`, `FillGaps(LightCurve, maxGapSize)`, `NormalizeFlux(LightCurve)` → returnerar ny `LightCurve`
- [ ] Enhetstester för samtliga datatyper och `LightCurveSanitizer`

---

## Phase 2 — Transitfysik (Physics.Astro)

Utöka `CSharpNumerics.Physics.Astro` med fysikaliska modeller för transit-geometri och limbmörkning. Dessa är domänoberoende fysik och tillhör Physics-sektionen.

- [ ] `TransitGeometry` — statisk klass i `Physics/Astro/`: `ImpactParameter(a, inclination, Rstar)`, `TransitProbability(a, Rstar, Rplanet)`, `TransitDuration(period, a, Rstar, Rplanet, inclination)`, `IngressDuration(...)`, `ContactTimes(...)` (T1–T4)
- [ ] `LimbDarkening` — statisk klass i `Physics/Astro/`: `Linear(mu, u1)`, `Quadratic(mu, u1, u2)`, `NonlinearFourParam(mu, c1..c4)`, `IntensityProfile(model, u_params, muArray)` — returnerar intensitetsprofil
- [ ] `TransitModel` — klass i `Physics/Astro/`: implementerar Mandel & Agol (2002) analytisk transitljuskurva. `Evaluate(double[] times, TransitParameters p, LimbDarkeningModel ld, StellarProperties star)` → `double[] modelFlux`. Stödjer cirkulär bana.
- [ ] `LimbDarkeningModel` (enum) i `Physics/Astro/` — `Uniform`, `Linear`, `Quadratic`, `NonlinearFourParam`
- [ ] `KeplerOrbit` — statisk klass i `Physics/Astro/`: `TrueAnomaly(meanAnomaly, eccentricity)` (Kepler-ekvationen, Newton-iteration), `SemiMajorAxis(period, stellarMass)` (Keplers tredje lag), `OrbitalVelocity(a, period)`
- [ ] Enhetstester: verifiera transitmodell mot kända Kepler-transitkurvor (analytiska testfall), limbmörkningsprofiler, Kepler-ekvationen

---

## Phase 3 — Transitdetektionspipeline

Kombinera befintliga `TimeSeriesAnalysis`-moduler till en sammanhållen detektionspipeline i `Engines/Exoplanet/Pipeline/`.

- [ ] `LightCurvePreprocessor` — orkestrerar: `Detrend(LightCurve, method)` → `RemoveOutliers()` → `Normalize()`. Delegerar till `TimeSeriesDetrending` och `LightCurveSanitizer`
- [ ] `PeriodSearcher` — wrapper kring BLS och Lomb-Scargle: `Search(LightCurve, minPeriod, maxPeriod, method)` → `PeriodSearchResult` (bästa period, FAP, spektrum). Stödjer `PeriodSearchMethod`-enum (`BLS`, `LombScargle`, `Both`)
- [ ] `TransitFitter` — anpassar `TransitModel` (Phase 2) till fasveckt ljuskurva via `NonlinearLeastSquaresFitter`. `Fit(LightCurve, trialPeriod, trialEpoch)` → `TransitFitResult` med `TransitParameters`, osäkerheter, χ²/BIC
- [ ] `TransitValidator` — utvärderar kandidater: kontrollerar SNR (transit-djup/brus), udda/jämn transit-djup-test, sekundär eklipsanalys, centrumfotoanalys. `Validate(TransitCandidate, LightCurve)` → `ValidationResult` med `IsValid`, `Warnings[]`, `Score`
- [ ] `TransitDetectionPipeline` — main orchestrator: `Detect(LightCurve, config)` → `TransitCandidate[]`. Steg: preprocessor → period search → phase fold → fit → validate → upprepa (multi-planet med iterativ subtraktion)
- [ ] `TransitDetectionConfig` — konfigurationsobjekt: `MinPeriodDays`, `MaxPeriodDays`, `MinTransitDepthPpm`, `SnrThreshold`, `MaxPlanets`, `DetrendingMethod`, `PeriodSearchMethod`
- [ ] Enhetstester: syntetisk ljuskurva med 1–3 injicerade transiter → pipeline ska återfinna samtliga med korrekt period (±1%)

---

## Phase 4 — Feature-extraktion

Extrahera transitspecifika features för ML-klassificering. Placeras i `Engines/Exoplanet/Features/`.

- [ ] `TransitFeatureExtractor` — extraherar features från en `TransitCandidate` + rå `LightCurve`: transit-djup, duration, ingress/egress-ratio, periodik-styrka (BLS SNR), udda/jämn djupskillnad, sekundäreklips-djup, ljuskurve-scatter, fasvecknings-χ², limbmörkningskoefficienter (anpassade), transitform-metric (V-shape vs U-shape)
- [ ] `TransitFeatureSet` — resultatobjekt: `Dictionary<string, double> Features`, med standardnamn: `Depth`, `Duration`, `Period`, `SnrBls`, `OddEvenRatio`, `VShapeMetric`, `SecondaryDepth`, `ScatterInTransit`, `ScatterOutTransit`, `LimbDarkeningU1`, `LimbDarkeningU2`, `IngressEgressRatio`
- [ ] `WindowedFeatureExtractor` — skapar `(Matrix X, VectorN y)` för ML-träning: kombinerar raw fasveckt ljuskurva (tidsfönster) med extraherade features som extra kolumner. Kompatibel med `SequenceDataHelper.CreateWindows`
- [ ] Enhetstester: verifiera extraherade features mot kända transiter

---

## Phase 5 — ML-träning & inference

Bygg tränings- och inference-workflow i `Engines/Exoplanet/Pipeline/`. Målet: träna en modell på labeled data → serialisera → ladda för prediktion.

- [ ] `TransitClassifierTrainer` — façade som sammanfogar feature-extraktion + `SupervisedExperiment`. `Train(LightCurve[] curves, TransitDisposition[] labels, TrainerConfig config)` → `TrainedTransitModel`. Grid-search över CNN1D, LSTM, BiLSTM med konfigurerbar hyperparameterrymd
- [ ] `TrainerConfig` — konfiguration: `int WindowSize`, `int Stride`, `ModelType[] CandidateModels`, `CrossValidatorConfig Cv`, `int Epochs`, `double ValidationSplit`
- [ ] `TrainedTransitModel` — wrapper: innehåller den tränade `IClassificationModel`, `TransitDetectionConfig`, feature-scaler-parametrar, träningsmetrik (accuracy, precision, recall, F1, confusion matrix). Stödjer `Predict(LightCurve)` → `TransitPrediction[]`
- [ ] `TransitPrediction` — `TransitCandidate Candidate`, `double Probability`, `TransitDisposition PredictedDisposition`
- [ ] `ModelSerializer` — serialisering/deserialisering av `TrainedTransitModel` till/från byte-array eller fil. Använder `FieldSerializer` från `Engines.Common` som bas. Format: JSON-metadata + binär viktdata
- [ ] `TransitInferencePipeline` — stateless inference: `LoadModel(byte[] serialized)` → `Predict(LightCurve)` → `TransitPrediction[]`. Designad för att anropas från webbtjänst. Trådsäker (inga mutable state)
- [ ] Enhetstester: full roundtrip — träna på syntetisk data → serialisera → deserialisera → prediktera → verifiera accuracy, precision, recall

---

## Phase 6 — Statistiska tillägg (Statistics)

Nya statistiska metoder som behövs för transitanalys men är generellt användbara. Placeras i respektive Statistics-namespace.

- [ ] `SigmaClipping` i `Statistics/Robust/` (`CSharpNumerics.Statistics.Robust`) — iterativ sigma-clipping: `Clip(double[] data, double sigmaLow, double sigmaHigh, int maxIter)` → `ClipResult` (mask, mean, std). Används av `LightCurveSanitizer.RemoveOutliers`
- [ ] `FalseAlarmProbability` i `Statistics/TimeSeriesAnalysis/` — FAP-beräkning för periodogram-toppar: analytisk (Baluev 2008) och bootstrap-baserad. `AnalyticalFAP(peakPower, numFrequencies, numDataPoints)`, `BootstrapFAP(times, values, peakPower, nBootstrap)`
- [ ] `SlidingWindowStatistics` i `Statistics/Robust/` (`CSharpNumerics.Statistics.Robust`) — löpande medelvärde, standardavvikelse, MAD (median absolute deviation) med konfigurbar fönsterstorlek. Används för scatter-beräkning i transit-features
- [ ] Enhetstester för samtliga statistiska tillägg

---

## Phase 7 — Engine-integration & demo

Integrera allt till `ExoplanetEngine : ISimulationEngine` och bygg en end-to-end demo.

- [ ] `ExoplanetEngine` — implementerar `ISimulationEngine`. `Init()` laddar konfiguration och eventuellt tränad modell. `Step(dt)` processar nästa ljuskurva-segment. Publicerar `TransitDetectedEvent` via `EventBus` vid fynd
- [ ] `TransitDetectedEvent` — eventtyp: `TransitCandidate Candidate`, `double Timestamp`
- [ ] `ExoplanetEngineConfig` — konfiguration: `TransitDetectionConfig Detection`, `TrainerConfig Training`, `string ModelPath` (för förtränad modell)
- [ ] End-to-end integrationstest: syntetisera TESS-liknande ljuskurvor (realistisk kadens, brus, systematics) → kör full pipeline (avtrend → BLS → fit → ML-klassificering) → verifiera att korrekt antal transiterande planeter detekteras
- [ ] Demo-test som visar hela workflow: data-inläsning → träning → inference → resultatutskrift. Dokumenterar det format som en webbtjänst behöver skicka/ta emot
- [ ] Uppdatera `Engines/README.md` med Exoplanet Engine-dokumentation: arkitekturdiagram, API-referens, kodexempel, dataformatbeskrivning

---

## Sammanfattning — Berörda namespaces

| Namespace | Typ av ändring |
|-----------|---------------|
| `Engines.Exoplanet` | **Ny** — motor, pipeline, features, datamodell |
| `Engines.Exoplanet.Data` | **Ny** — `LightCurve`, `TransitCandidate`, `StellarProperties`, etc. |
| `Engines.Exoplanet.Pipeline` | **Ny** — detektion, träning, inference |
| `Engines.Exoplanet.Features` | **Ny** — feature-extraktion |
| `Physics.Astro` | **Utökad** — `TransitGeometry`, `LimbDarkening`, `TransitModel`, `KeplerOrbit`, `LimbDarkeningModel` (enum) |
| `Statistics.Robust` | **Ny** — `SigmaClipping` |
| `Statistics` | **Ny** — `SlidingWindowStatistics` |
| `Statistics.TimeSeriesAnalysis` | **Utökad** — `FalseAlarmProbability` |
| `ML.Sequence` | Befintlig — används som är |
| `Statistics.Fitting` | Befintlig — används som är |
| `Engines.Common` | Befintlig — `ISimulationEngine`, `EventBus`, `FieldSerializer` |

---

## Dataflöde — Arkitekturskiss

```
               ┌─────────────────────────────────────────────┐
               │           EXTERN WEBBTJÄNST                 │
               │  (hämtar data från NASA/MAST/TESS API)      │
               └──────────────────┬──────────────────────────┘
                                  │ LightCurve (JSON/CSV)
                                  ▼
               ┌─────────────────────────────────────────────┐
               │         Engines.Exoplanet                   │
               │                                             │
               │  ┌─────────────┐   ┌────────────────────┐  │
               │  │ LightCurve  │──▶│ LightCurvePreproc  │  │
               │  │ Sanitizer   │   │ (Detrend/Normalize)│  │
               │  └─────────────┘   └────────┬───────────┘  │
               │                             │               │
               │                ┌────────────▼────────────┐  │
               │                │    PeriodSearcher       │  │
               │                │  (BLS / Lomb-Scargle)   │  │
               │                └────────────┬────────────┘  │
               │                             │               │
               │                ┌────────────▼────────────┐  │
               │                │    TransitFitter        │  │
               │                │  (Mandel-Agol via LM)   │  │
               │                └────────────┬────────────┘  │
               │                             │               │
               │       ┌────────────────┐    │               │
               │       │FeatureExtractor│◀───┘               │
               │       └───────┬────────┘                    │
               │               │                             │
               │  ┌────────────▼──────────────────────┐      │
               │  │  ML Inference (CNN1D/LSTM/BiLSTM) │      │
               │  │  via TrainedTransitModel           │      │
               │  └────────────┬──────────────────────┘      │
               │               │                             │
               │  ┌────────────▼────────────┐                │
               │  │  TransitValidator       │                │
               │  │  (SNR, odd/even, etc.)  │                │
               │  └────────────┬────────────┘                │
               │               │                             │
               │               ▼                             │
               │       TransitPrediction[]                   │
               │       (→ JSON till webbtjänst)              │
               └─────────────────────────────────────────────┘

Beroenden:
  Engines.Exoplanet → Physics.Astro (TransitModel, LimbDarkening, KeplerOrbit)
  Engines.Exoplanet → Statistics.Fitting (NonlinearLeastSquaresFitter)
  Engines.Exoplanet → Statistics.TimeSeriesAnalysis (BLS, LombScargle, Detrending)
  Engines.Exoplanet → ML.Sequence (CNN1DClassifier, LSTMClassifier, BiLSTMClassifier)
  Engines.Exoplanet → Engines.Common (ISimulationEngine, EventBus, FieldSerializer)
```
