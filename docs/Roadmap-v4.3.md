# v4.3 — Fundamentet betalar av sig

> **Status (2026-09-29):** Planerad, ej påbörjad. Utgår från `master` efter `v4.2.0`.

## Mål

v4.1–4.2 byggde **dekompositionsfundamentet** (LU, Cholesky, QR, egendekomposition) och levererade en
stor fysikutbyggnad. Fundamentet står — men **vinsten togs aldrig ut**: `Matrix.Inverse` går via LU, medan
resten av kodbasen fortfarande löser ekvationssystem på egen hand i åtta filer. Löftet i
[LinearAlgebraRoadMap](LinearAlgebraRoadMap.md) — "nästan allt annat blir bättre av det" — är alltså
ännu inte infriat.

v4.3 infriar det. Releasen har tre spår: **migrera anropsställena** till dekompositionerna (bättre
numerik, mindre kod), **bygga rotfinnarmodulen** (ny publik yta som ersätter en svag befintlig
implementation), och **lägga mätgrunden** inför prestandaarbetet i v4.4.

Inga breaking changes. Allt sker bakom befintliga API:er.

---

## Nulägesanalys — Befintlig arkitektur

### Vad som finns idag

| Komponent | Status | Plats | Kommentar |
|-----------|--------|-------|-----------|
| `LuDecomposition`, `CholeskyDecomposition` | ✓ | `Numerics/LinearAlgebra/Decompositions/` | Levererade i v4.1 |
| `QrDecomposition`, `EigenDecomposition` | ✓ | `Numerics/LinearAlgebra/Decompositions/` | Levererade i v4.1 |
| `Matrix.Inverse`, `LinearSystemSolver` | ✓ | `Numerics/Objects/`, `DifferentialEquationExtensions.cs` | Migrerade till LU |
| Kvantmodulens egenlösare | ✓ | `Physics/Quantum/SchrodingerExtensions.cs` | Migrerad — mönstret att följa |
| `SparseMatrix.SolvePCG` | ✓ | `Numerics/Objects/SparseMatrix.cs` | Preconditionerad CG; driver `Assembler2D` |
| Handskriven gausselimination | ⚠ | 8 filer — se tabell nedan | Duplicerar `LuDecomposition` |
| `CoupledOscillators.JacobiEigen` | ⚠ | `Physics/Mechanics/Oscillations/` | Duplicerar `EigenDecomposition` |
| `PCA` | ⚠ | `ML/DimensionalityReduction/Algorithms/` | Potensmetod med deflation |
| `NewtonRaphson` | ⚠ | `Numerics/NumericExtensions.cs` | Enda rotfinnaren — se begränsning 2 |
| `NaiveBayes.NumClasses` | ⚠ | `ML/Models/Classification/` | Kastar `NotImplementedException` |
| Bisection / Secant / Brent | ✗ | — | Finns ej |
| Benchmark-projekt | ✗ | — | Finns ej — inga mätningar alls |
| SVD | ✗ | — | Finns ej (planerad v4.4) |

### Nyckelidentifierade begränsningar

**1. Åtta filer löser täta system med egen kod.** Dekompositionerna finns, men anropsställena migrerades
aldrig. Utöver duplicerad kod innebär det att varje lösning räknas om från grunden — ingen faktorisering
återanvänds mellan högerled.

| Fil | Vad den gör idag | Åtgärd |
|-----|------------------|--------|
| `Statistics/Fitting/FittingSolver.cs` | Normalekvationer (AᵀA) + Gauss-Jordan-invers | → QR |
| `Statistics/InferentialStatisticsExtensions.cs` | Augmenterad matris, egen pivotering | → LU |
| `Numerics/Interpolation/MultivariateInterpolation.cs` | RBF-system Φw = f | → LU |
| `Physics/FluidDynamics/Aerodynamics/PanelMethod.cs` | Tät influensmatris | → LU (faktorisera en gång) |
| `Numerics/FiniteElement/Assembler1D.cs` | Egen gausselimination | → LU |
| `Numerics/Interpolation/CubicSplineInterpolation.cs` | Dense-fallback för n ≥ 4 | → LU (Thomas-grenen **behålls**) |
| `Numerics/DifferentialEquationExtensions.cs` | 3×3-gausselimination vid sidan av LU-fasaden | → LU |
| `Physics/Mechanics/Oscillations/CoupledOscillators.cs` | Privat Jacobi-egenlösare | → `EigenDecomposition` |

**2. Bibliotekets enda rotfinnare är opålitlig.** `NewtonRaphson` (`NumericExtensions.cs:156`) kör exakt
100 iterationer varje gång — ingen konvergenskontroll, inget tidigt avbrott, inget skydd mot f′(x) ≈ 0
(division ger tyst `Inf`/`NaN`), ingen tolerans- eller maxiterationsparameter, och ingen signal om att
den *inte* konvergerade. Derivatan hämtas via `DerivativeExtensions.Derivate`, som bygger en ny
Pascal-matris per anrop — alltså 100 matrisallokeringar per rotsökning oavsett om problemet löstes på
iteration 4. Dessutom finns egna Newton-loopar i `KeplerOrbit` och `Gravitation/LagrangePoints`.

**3. Normalekvationer i fitting-vägen.** `FittingSolver` bygger AᵀA och inverterar. Konditionstalet för
AᵀA är kvadraten på A:s — för Vandermonde-designmatriser (polynomanpassning) är det en reell
noggrannhetsförlust. QR löser samma problem utan att kvadrera konditionen. Påverkar `LeastSquaresFitter`,
`WeightedLeastSquaresFitter`, `NonlinearLeastSquaresFitter`, `RobustFitter` och `ParameterEstimation`.

**4. Kalman-vinsten beräknas via explicit invers.** `KalmanFilter`, `ExtendedKalmanFilter` och
`KalmanSmoother` gör alla `S.Inverse()` på innovationskovariansen. S är symmetrisk positivt definit —
Cholesky-lösning är både snabbare och stabilare, och bevarar symmetrin som explicit invertering kan bryta.

**5. Ingenting är mätt.** Det finns inget benchmark-projekt. Spår 1 ändrar heta kodvägar och v4.4 planerar
SIMD — utan baseline går ingen förbättring att bevisa. Bibliotekets egen princip är *mät före optimering*.

---

## Del 1 — Ta ut vinsten från dekompositionerna

### Vad som krävs att bygga

| Komponent | Beskrivning | Insats |
|-----------|-------------|--------|
| **Migrering av 8 anropsställen** | Ersätt handskriven eliminering med `LuDecomposition`/`QrDecomposition` enligt tabellen ovan | **Låg-medel** |
| **`FittingSolver` → QR** | Ersätt normalekvationer + Gauss-Jordan med QR-lösning. Behåll `(XᵀX)⁻¹`-diagonalen för standardfel (via R) | **Medel** |
| **Kalman → Cholesky** | `K = P·Hᵀ·S⁻¹` löses som `S·Kᵀ = H·Pᵀ` i stället för explicit invers | **Låg** |
| **`PCA` → `EigenDecomposition`** | Ersätt potensmetod + deflation. Behåll dual-PCA-grenen för n < d | **Låg-medel** |
| **`CoupledOscillators` → `EigenDecomposition`** | Ta bort privat `JacobiEigen`, samma mönster som kvantmodulen | **Låg** |
| **Avveckla dubbla egenvärdes-API:er** | `EigenValues`/`DominantEigenVector`/`EigenVector` i `DifferentialEquationExtensions` delegerar till `EigenDecomposition` | **Låg** |

### Designprinciper

- **Inga signaturändringar.** Varje migrering sker inuti befintliga metoder; publikt API är oförändrat.
- **Faktorisera en gång.** Där flera högerled löses mot samma matris (PanelMethod, FEM, Kalman-loopar)
  ska faktoriseringen cachas — det är hela poängen med dekompositionsklasserna.
- **Numeriskt befogade undantag behålls.** `CubicSpline`s Thomas-algoritm för tridiagonala system är
  O(n) mot LU:s O(n³) och ska inte migreras. Bara dense-fallbacken byts.
- **Regressionstester först.** Varje migrering föregås av ett test som låser nuvarande resultat, så att
  bytet bevisligen inte ändrar utfallet (utöver förväntad noggrannhetsförbättring).

### Direkta vinster i befintlig kod

| Befintlig funktion | Förbättring |
|--------------------|-------------|
| Polynomanpassning (`LeastSquaresFitter`) | QR undviker kvadrerat konditionstal på Vandermonde-matriser |
| `Ridge`, `Linear`, `ElasticNet` | Stabilare lösning på illa skalad data |
| `KalmanFilter` m.fl. | Garanterad symmetri i kovariansuppdateringen, färre flops |
| `PanelMethod` | Faktorisering återanvänds över flera anfallsvinklar |
| `PCA` | Deterministisk och exakt i stället för iterativ med toleranströskel |
| Hela kodbasen | ~8 kopior av samma algoritm försvinner |

---

## Del 2 — Rotfinnarmodulen

### Vad som krävs att bygga

| Komponent | Beskrivning | Insats |
|-----------|-------------|--------|
| **`Bisection`** | Garanterad konvergens givet teckenväxling på `[a, b]` | **Trivial** |
| **`Secant`** | Derivatafri Newton-variant | **Trivial** |
| **`Brent`** | Industristandard: bisektion + sekant + invers kvadratisk interpolation | **Låg-medel** |
| **`Newton`** | Riktig implementation: tolerans, maxiterationer, skydd mot f′ ≈ 0, valfri analytisk derivata | **Låg** |
| **`RootResult`** | Returtyp med `Value`, `Converged`, `Iterations`, `Residual` — samma mönster som `Optimization/Strategies/` | **Låg** |
| **`FindRoot`-fasad** | `Func<double,double>.FindRoot(a, b)` i samma stil som befintliga extensions | **Låg** |

Placeras i `Numerics/RootFinding/`.

### Designprinciper

- **`NewtonRaphson` blir kvar som wrapper** över nya `Newton` med samma signatur och defaultbeteende —
  ingen breaking change, men den slutar allokera 100 Pascal-matriser per anrop.
- **Analytisk derivata som overload.** Nuvarande beteende (finita differenser) behålls som default,
  men `Newton(f, df, x0)` låter anroparen slippa derivataapproximationen helt.
- **Konvergens rapporteras, tystnas inte.** `RootResult.Converged` gör det möjligt för anropare att
  upptäcka misslyckanden — idag är det omöjligt.

### Direkta vinster i befintlig kod

| Befintlig funktion | Förbättring |
|--------------------|-------------|
| `Physics/Astro/KeplerOrbit` | Egen Newton-loop för Keplers ekvation → `Newton`/`Brent` |
| `Physics/Gravitation/LagrangePoints` | Egen Newton-metod för kollineära punkter → samma fasad |
| `NewtonRaphson`-anropare | Tidigt avbrott vid konvergens i stället för 100 fasta iterationer |

---

## Del 3 — Mätgrunden

### Vad som krävs att bygga

| Komponent | Beskrivning | Insats |
|-----------|-------------|--------|
| **`Numerics.Benchmarks`-projekt** | BenchmarkDotNet, eget projekt i lösningen, ingår ej i NuGet-paketet | **Låg** |
| **Baseline: linjär algebra** | Matris×matris, matris×vektor, LU/QR/Cholesky-faktorisering, `Solve` | **Låg** |
| **Baseline: före/efter Del 1** | Samma mätning körd före och efter migreringen — belägger att bytet inte kostar prestanda | **Låg** |
| **Baseline: ML-träningsloop** | En epok MLP-träning, som referenspunkt inför v4.4 | **Låg-medel** |
| **Math.NET-jämförelse** | Isolerad i benchmark-projektet — kärnbiblioteket förblir beroendefritt | **Låg** |

### Designprinciper

- **Beroenden isoleras.** BenchmarkDotNet och Math.NET hamnar enbart i benchmark-projektet
  (princip 3 i höstroadmapen).
- **Benchmarks körs inte i CI per commit.** De är för långsamma; körs manuellt inför release och
  resultaten checkas in som markdown.

---

## Del 4 — Kvarvarande städning

Punkter som stod i v4.1-scopet och ännu inte är gjorda:

| Punkt | Plats |
|-------|-------|
| `NaiveBayes.NumClasses` kastar `NotImplementedException` | `ML/Models/Classification/NaiveBayes.cs:18` — enda i hela kodbasen; sätts i `Fit` som alla andra klassificerare |
| Tempfilsreferens i csproj | `<None Remove="Numerics\NumericExtensions.cs~RF43b9e4.TMP" />` |
| Dubblerade roadmaps | `AdvancedGameEngineRoadMap`, `ExoplanetEngineRoadMap`, `Multiphysics-Roadmap`, `TerrainSpreadRoadMap` finns i både `docs/` och `docs/completed/`. **Obs:** Multiphysics-kopian i `completed/` är den *äldre* (markerar FEM som "deferred") — `docs/`-versionen är den korrekta |
| Oanvänt fält | `SARSA._pendingNextAction` (CS0169-varning) |

---

## Implementationsplan — Faser

### Phase 1 — Rotfinnare ✔ klar
- [x] Skapa `Numerics/RootFinding/`-struktur + `RootResult`
- [x] Implementera `Bisection`, `Secant`
- [x] Implementera `Brent`
- [x] Implementera `Newton` med tolerans, maxiter, f′-skydd och analytisk-derivata-overload
- [x] `NewtonRaphson` blir wrapper över `Newton` (signatur oförändrad)
- [x] Migrera `KeplerOrbit` och `LagrangePoints` till fasaden
- [x] Enhetstester: patologiska funktioner, platta derivator, ingen teckenväxling, konvergensrapportering

> **Noterat under Phase 1:** `TimeserieValidationTests` har tre fel som *inte* rör rotfinnarna.
> `TestData/CsvTestDataGenerator` skriver decimaltal med aktuell kultur (`$"{v:F2}"`), så på en svensk
> maskin blir `6.44` till `6,44` och kolliderar med CSV-avgränsaren. Testerna passerar med
> `DOTNET_SYSTEM_GLOBALIZATION_INVARIANT=1` och på CI (Linux). Fixen är `CultureInfo.InvariantCulture`
> i generatorn — ligger utanför v4.3-scopet men bör tas någon gång.

### Phase 2 — Benchmark-baseline
- [x] Skapa `Numerics.Benchmarks`-projekt (BenchmarkDotNet), lägg till i `.sln`, exkludera från paketering
- [x] Benchmarks för matmul, SpMV, LU/QR/Cholesky, `Solve`
- [x] Benchmark för en MLP-träningsepok
- [x] Kör och checka in baseline **före** Del 1-migreringen → [baseline-2026-10-01](benchmarks/baseline-2026-10-01.md)

> **Noterat under Phase 2:** Windows app control-policyn (0x800711C7) blockerar den DLL-kopia som
> BenchmarkDotNet lägger i sin genererade per-benchmark-mapp, så standardtoolchainen (en process per
> benchmark) ger `No Workload Results` på den här maskinen. Körningarna görs därför med `--inProcess`.
> Mätningarna är giltiga men något mindre isolerade än med processeparation — jämför alltid baslinje
> mot omkörning gjord på samma sätt.

### Phase 3 — Migrering till dekompositionerna
- [x] Regressionstester för anropsställena — nya sviter för `PanelMethod` och `PCA`, som saknade
      täckning helt; övriga sites täcks av befintliga tester
- [x] Migrera `MultivariateInterpolation`, `Assembler1D`, `CubicSpline`-fallback,
      `InferentialStatisticsExtensions`, `DifferentialEquationExtensions` till LU
- [x] Migrera `FittingSolver` till QR + verifiera standardfelen mot nuvarande värden — gjort för
      `LeastSquaresFitter`, `WeightedLeastSquaresFitter`, `RobustFitter`, `ParameterEstimation` och
      `NonlinearLeastSquaresFitter`s kovarians. LM-iterationen lämnas på dämpade normalekvationer
      (dämpningen läggs på `JtJ`s diagonal; QR kräver den utvidgade formuleringen [J; √(λD)], en
      algoritmändring)
- [x] Migrera `KalmanFilter`/`ExtendedKalmanFilter`/`KalmanSmoother` till Cholesky-lösning
- [x] Migrera `CoupledOscillators` till `EigenDecomposition`
- [x] Migrera `PCA` till `EigenDecomposition`
- [x] Cacha faktoriseringar där flera högerled löses mot samma matris — gäller Kalman-vinsten,
      där `Cholesky.Solve(Matrix)` nu löser alla kolumner mot en faktorisering

**Två punkter ströks efter att koden lästs:**

- **`PanelMethod` migrerades inte.** Den har varken tester eller anropare. Karaktäriseringstesterna
  avslöjade att lösaren ger Cp ≈ 1 på varje panel för en cylinder vid noll anfallsvinkel, där
  potentialflöde ger ett intervall över [−3, 1]. Om det är en bugg eller om en cylinder ligger
  utanför dess domän (den påtvingar ett Kutta-villkor vid en bakkant cylindern inte har) är
  ouppklarat — testet ligger `[Ignore]`-markerat med frågan nedskriven. Utan tester och utan
  anropare är dess duplicering den minst skadliga, och att utreda lösaren är separat arbete.
- **`EigenValues`/`DominantEigenVector`/`EigenVector` delegerar inte.** Planpunkten antog att de var
  duplicerade generella egenlösare. De är i stället lågnoggranna hjälpfunktioner vars publika
  kontrakt är avrundade heltalskvoter — `Math.Abs(Math.Round(c / min))` — asserterat av befintliga
  tester (`result[0] == 2`) och konsumerat av `OdeSolver`. Att delegera dem byter ut kontraktet mot
  normaliserade egenvektorer, vilket är en API-ändring med omvalidering av `OdeSolver`, inte en
  refaktorering bakom befintligt API. Kräver ett eget beslut.

### Phase 4 — Städning
- [ ] `NaiveBayes.NumClasses` sätts i `Fit`
- [ ] Ta bort tempfilsreferensen i csproj
- [ ] Ta bort dubblerade roadmaps ur `docs/` respektive `docs/completed/`
- [ ] Ta bort `SARSA._pendingNextAction`

### Phase 5 — Verifiering & release
- [ ] Kör om benchmarks och jämför mot baseline från Phase 2
- [ ] Uppdatera README med rotfinnar-exempel
- [ ] Uppdatera `LinearAlgebraRoadMap` Phase 5 → klar
- [ ] Versionsbump till 4.3.0 + tagg

---

## Sammanfattning

| Del | Genomförbarhet | Insats | Största risk |
|-----|---------------|--------|--------------|
| **Rotfinnare** | **Hög** | Låg | Ingen — fristående ny kod |
| **Benchmarks** | **Hög** | Låg | Ingen — påverkar inte biblioteket |
| **Migrering till dekompositioner** | **Hög** | Medel | Beteendeförändring i befintliga resultat — mitigeras av regressionstester först |
| **Städning** | **Hög** | Trivial | Ingen |

**Rekommendation:** Kör faserna i ordning. Rotfinnarna först — de är fristående, ger snabb synlig
leverans och rör ingen befintlig kod. Benchmarks näst, så att baseline finns *innan* migreringen ändrar
heta kodvägar. Migreringen sist och filvis, med regressionstester före varje byte: det är den enda delen
som kan ändra numeriska resultat, och den ska kunna backas per fil om något ser fel ut.

**Vad som medvetet skjuts till v4.4:** SVD (med pseudoinvers och konditionstal), SIMD-kärnor,
dual numbers. SVD hör tematiskt hemma här men är det enskilt största arbetet och tjänar på att
QR-maskineriet är inkört. Gradient boosting är oberoende av allt annat och kan flyttas in om v4.3
behöver en mer användarsynlig funktion — se [MLExpansionRoadMap](MLExpansionRoadMap.md).
