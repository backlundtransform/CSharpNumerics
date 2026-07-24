# Linjär algebra — Dekompositioner, glesa lösare & rotfinnare

## Mål

Bygga ut linjär algebra-fundamentet med **matrisdekompositioner** (LU, QR, Cholesky, SVD, egendekomposition), **iterativa glesa lösare** (CG, BiCGSTAB, GMRES) samt ett komplett **rotfinnar-modul**. Detta är den enskilt viktigaste investeringen för att göra CSharpNumerics till ett ledande ramverk — nästan alla andra delar (PCA, Ridge, minsta kvadrat, FEM, Kalman, kvantmodulen) blir bättre av den.

---

## Nulägesanalys — Befintlig arkitektur

### Vad som finns idag

| Komponent | Status | Plats | Kommentar |
|-----------|--------|-------|-----------|
| `Matrix` | ✓ | `Numerics/Objects/Matrix.cs` | Inverse, Determinant, Transpose, Slice, operatorer |
| `SparseMatrix` | ✓ | `Numerics/Objects/SparseMatrix.cs` | Lagring finns, men inga iterativa lösare |
| `GaussElimination` / `LinearSystemSolver` | ✓ | `Numerics/DifferentialEquationExtensions.cs` | Direkt lösning av täta system |
| Egenvärden (potensmetod) | ✓ | `Numerics/DifferentialEquationExtensions.cs` | Endast dominant + iterativ full |
| Symmetrisk egenlösare | ✓ | `Physics/Quantum/` | Inlåst i kvantmodulen — bör lyftas ut |
| `NewtonRaphson` | ✓ | `Numerics/NumericExtensions.cs` | Enda rotfinnaren |
| LU-dekomposition | ✗ | — | Finns ej |
| QR-dekomposition | ✗ | — | Finns ej |
| Cholesky | ✗ | — | Finns ej |
| SVD | ✗ | — | Finns ej |
| Explicit egendekomposition | ✗ | — | Finns ej som fristående API |
| CG/BiCGSTAB/GMRES | ✗ | — | Finns ej |
| Brent/bisektion/sekant | ✗ | — | Finns ej |

### Nyckelidentifierade begränsningar

1. **`Matrix.Inverse` via kofaktorer/eliminering utan pivotering-API** — Ingen möjlighet att återanvända faktorisering för flera högerled; `Solve(A, b)` beräknas från scratch varje gång.
2. **PCA använder inte SVD** — PCA via kovariansmatris + potensmetod är numeriskt sämre än SVD-baserad PCA.
3. **`SparseMatrix` är en ren datastruktur** — FEM-assemblern bygger glesa system men löser dem tätt, vilket sätter ett hårt tak på problemstorlek.
4. **Konditionstal, rank, pseudoinvers saknas** — kräver SVD.

---

## Del 1 — Dekompositioner (täta matriser)

### Vad som krävs att bygga

| Komponent | Beskrivning | Algoritm | Insats |
|-----------|-------------|----------|--------|
| **`LuDecomposition`** | `A = P·L·U` med partiell pivotering. `Solve(b)`, `Determinant`, `Inverse` | Doolittle med radpivotering | **Låg-medel** |
| **`CholeskyDecomposition`** | `A = L·Lᵀ` för symmetriska positivt definita. `Solve(b)`, `IsPositiveDefinite` | Cholesky–Banachiewicz | **Låg** |
| **`QrDecomposition`** | `A = Q·R`. `Solve(b)` (minsta kvadrat för överbestämda system) | Householder-reflektioner | **Medel** |
| **`EigenDecomposition`** | Egenvärden + egenvektorer. Symmetrisk: Jacobi eller QR med Wilkinson-shift. Osymmetrisk: Hessenberg + QR-iteration | QR-algoritm | **Medel-hög** |
| **`SvdDecomposition`** | `A = U·Σ·Vᵀ`. Ger `Rank`, `ConditionNumber`, `PseudoInverse` | Golub–Kahan bidiagonalisering + QR | **Hög** |
| **`Matrix.Decompose()`-fasad** | Extension-metoder: `matrix.Lu()`, `matrix.Qr()`, `matrix.Cholesky()`, `matrix.Svd()`, `matrix.Eigen()` | — | **Låg** |

### Designprinciper

- Placeras i `Numerics/LinearAlgebra/Decompositions/`.
- Varje dekomposition är en **klass som cachar faktoriseringen** — `var lu = A.Lu(); lu.Solve(b1); lu.Solve(b2);` återanvänder O(n³)-arbetet.
- `Matrix.Inverse` och `LinearSystemSolver` refaktoreras internt till att använda LU (ingen breaking change i publikt API).
- Symmetriska egenlösaren i `Physics/Quantum/` ersätts med den nya `EigenDecomposition` (kvantmodulen blir konsument, inte ägare).

### Direkta vinster i befintlig kod

| Befintlig funktion | Förbättring |
|--------------------|-------------|
| `PCA` (`ML/`) | SVD-baserad → numeriskt stabil, hanterar fler features än samples |
| `Ridge`/`LinearRegression` | Cholesky/QR-lösning i stället för normalekvationer + invers |
| `KalmanFilter` (`Statistics/`) | Cholesky för kovariansuppdatering → garanterad symmetri |
| `NonlinearFitting` | QR för Jacobian-system i Levenberg–Marquardt-stil |
| `Schrodinger`-lösaren | Generell symmetrisk egenlösare |

---

## Del 2 — Iterativa glesa lösare

### Vad som krävs att bygga

| Komponent | Beskrivning | Insats |
|-----------|-------------|--------|
| **`ConjugateGradient`** | För symmetriska positivt definita glesa system (FEM-styvhetsmatriser!) | **Låg-medel** |
| **`BiCgStab`** | För osymmetriska system | **Medel** |
| **`Gmres`** | Restarted GMRES(m) för generella system | **Medel-hög** |
| **Preconditioners** | Jacobi (diagonal), ILU(0), IC(0) | **Medel** |
| **`SparseMatrix` CSR-format** | Nuvarande lagring kompletteras med Compressed Sparse Row för snabb SpMV (matris×vektor är den heta loopen i alla iterativa lösare) | **Medel** |
| **`IIterativeSolverOptions`** | Tolerans, max-iterationer, konvergenshistorik — samma mönster som `Optimization/Strategies/` | **Låg** |

### Varför det är relevant

`FiniteElement/`-assemblern (1D/2D) bygger glesa styvhetsmatriser som idag löses med tät Gauss-eliminering — O(n³) och O(n²) minne. Med CG + ILU-preconditioner skalar samma problem till 100 000+ frihetsgrader. Samma sak gäller `FiniteDifference/`-diskretiseringarna (Poisson, värmeledning, Navier–Stokes-trycksteg).

---

## Del 3 — Rotfinnare

### Vad som krävs att bygga

| Komponent | Beskrivning | Insats |
|-----------|-------------|--------|
| **`Bisection`** | Garanterad konvergens givet teckenväxling | **Trivial** |
| **`Secant`** | Derivatafri Newton-variant | **Trivial** |
| **`Brent`** | Industristandard: kombinerar bisektion, sekant, invers kvadratisk interpolation | **Låg-medel** |
| **`PolynomialRoots`** | Companion-matris + egendekomposition (återanvänder Del 1!) | **Låg** (givet `EigenDecomposition`) |
| **`RootFinder`-fasad** | `Func<double,double>.FindRoot(a, b)` extension i samma stil som `NewtonRaphson` | **Låg** |

Placeras i `Numerics/RootFinding/`. Direkt användbart i: `KeplerOrbit` (Keplers ekvation löses idag med handrullad iteration), `TransitGeometry`, optikmodulens stråle–yta-skärningar, `Bernoulli`-lösningar i fluiddynamiken.

---

## Implementationsplan — Faser

### Phase 1 — Dekompositionsgrund
- [ ] Skapa `Numerics/LinearAlgebra/Decompositions/`-struktur
- [ ] Implementera `LuDecomposition` med partiell pivotering + `Solve`/`Determinant`/`Inverse`
- [ ] Implementera `CholeskyDecomposition` + `IsPositiveDefinite`
- [ ] Refaktorera `Matrix.Inverse` och `LinearSystemSolver` till LU internt (inga API-ändringar)
- [ ] Enhetstester: kända faktoriseringar, singulära matriser, round-trip `A ≈ P·L·U`

### Phase 2 — QR & egendekomposition
- [ ] Implementera `QrDecomposition` (Householder) + minsta kvadrat-`Solve`
- [ ] Implementera symmetrisk `EigenDecomposition` (QR med shift)
- [ ] Implementera osymmetrisk egenlösare (Hessenberg + QR-iteration)
- [ ] Migrera kvantmodulens symmetriska egenlösare till den nya
- [ ] Enhetstester: ortogonalitet `QᵀQ = I`, egenpar-residualer `‖Av − λv‖`

### Phase 3 — SVD
- [ ] Implementera `SvdDecomposition` (Golub–Kahan)
- [ ] `Rank`, `ConditionNumber`, `PseudoInverse`, `Nullspace`
- [ ] Refaktorera `PCA` till SVD-baserad
- [ ] Enhetstester: rekonstruktion `A ≈ UΣVᵀ`, jämförelse mot kända referensvärden

### Phase 4 — Glesa lösare
- [ ] CSR-lagring i `SparseMatrix` + snabb SpMV
- [ ] Implementera `ConjugateGradient` + Jacobi-preconditioner
- [ ] Implementera `BiCgStab` och `Gmres(m)`
- [ ] ILU(0)/IC(0)-preconditioners
- [ ] Koppla in i `FiniteElement/`-lösningsvägen (opt-in via options)
- [ ] Enhetstester + konvergenstester på FEM-genererade system

### Phase 5 — Rotfinnare
- [ ] Implementera `Bisection`, `Secant`, `Brent`
- [ ] Implementera `PolynomialRoots` via companion-matris
- [ ] Migrera `KeplerOrbit`s ekvationslösning till `Brent`/`NewtonRaphson`-fasaden
- [ ] Enhetstester: Wilkinson-polynom, patologiska funktioner

### Phase 6 — Integration & dokumentation
- [ ] Uppdatera README med dekompositions-exempel
- [ ] Benchmarks mot Math.NET Numerics (se PerformanceRoadMap)
- [ ] Exempelkod: minsta kvadrat via QR, PCA via SVD, FEM med CG

---

## Sammanfattning

| Del | Genomförbarhet | Insats | Största risk |
|-----|---------------|--------|--------------|
| **LU/Cholesky/QR** | **Hög** | Medel | Liten — väldokumenterade algoritmer |
| **Eigen/SVD** | **Medel-hög** | Hög | Konvergens- och stabilitetsdetaljer i QR-iterationen |
| **Glesa lösare** | **Hög** | Medel | Preconditioner-kvalitet avgör praktisk nytta |
| **Rotfinnare** | **Hög** | Låg | Ingen |

**Rekommendation:** Kör faserna i ordning — LU/Cholesky ger snabba vinster och bygger vana vid mönstret; SVD är svårast och tjänar på att QR-maskineriet redan finns. Del 2 och Del 3 är oberoende av varandra och kan paralleliseras eller flyttas om.
