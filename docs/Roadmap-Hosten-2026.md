# CSharpNumerics — Övergripande roadmap hösten 2026

> **Status (2026-08-29):** Ur v4.1-scopet är LinearAlgebra Phase 1–2 klart (LU/Cholesky/QR/egendekomposition,
> refaktorerad `Matrix.Inverse`/`LinearSystemSolver`, kvantmodulen migrerad) — se statusnoten i
> [LinearAlgebraRoadMap](LinearAlgebraRoadMap.md). Övriga v4.1-punkter (benchmarks, SIMD, rotfinnare,
> dual numbers, städning) är ej påbörjade. Dokumentet beskriver den ursprungliga höstplanen; vid nya mål
> för biblioteket, uppdatera eller ersätt planen härifrån.

## Vision

CSharpNumerics ska bli ett av de främsta ramverken för **maskininlärning och numerisk analys i C#**. Bredden finns redan (numerik, statistik, ML inkl. RL och sekvensmodeller, samt en unik fysikmodul). Höstens arbete fokuserar på det som skiljer ett bra bibliotek från ett ledande ramverk: **fundament, prestanda och differentiering**.

Den strategiska tesen: inget annat .NET-bibliotek har både fysiken och ML:en i samma paket. Autodiff + fysikmodulen → **Physics-Informed Neural Networks** som flaggskepps-feature gör CSharpNumerics till det självklara valet för vetenskaplig ML i C#.

---

## Delroadmaps

| Roadmap | Innehåll | Roll |
|---------|----------|------|
| [LinearAlgebraRoadMap](LinearAlgebraRoadMap.md) | LU/QR/Cholesky/SVD/eigen, glesa iterativa lösare (CG/BiCGSTAB/GMRES), rotfinnare (Brent m.fl.) | **Fundamentet** — nästan allt annat blir bättre av det |
| [PerformanceRoadMap](PerformanceRoadMap.md) | BenchmarkDotNet-projekt, SIMD/TensorPrimitives, allokeringsfria kärnor, generic math (`INumber<T>`) | **Trovärdigheten** — ledande ramverk visar siffror |
| [AutodiffRoadMap](AutodiffRoadMap.md) | Dual numbers (forward mode), tape (reverse mode), NN på autodiff, PINN-showcase | **Differentieringen** — nischen ingen annan i .NET fyller |
| [MLExpansionRoadMap](MLExpansionRoadMap.md) | Gradient boosting, Transformer/attention, ONNX-export, t-SNE/UMAP | **Efterfrågan** — de mest saknade algoritmerna |

---

## Releaseplan

### v4.1 — "Fundamentet" (september–oktober)

| Innehåll | Källa |
|----------|-------|
| LU, Cholesky, QR + refaktorerad `Matrix.Inverse`/`LinearSystemSolver` | LinearAlgebra Phase 1–2 |
| Egendekomposition (symmetrisk + osymmetrisk) | LinearAlgebra Phase 2 |
| `Numerics.Benchmarks`-projekt + baseline + Math.NET-jämförelser | Performance Phase 1 |
| SIMD-vektoriserade kärnor (`VectorKernels`, `MatMulKernel`) | Performance Phase 2 |
| Rotfinnare: Bisection, Secant, Brent | LinearAlgebra Phase 5 |
| Dual numbers + `AutoDerivative`-extensions | Autodiff Phase 1 |
| Städning: `.sln`-referenser till borttagna Engines-projekt, `NaiveBayes`-stub, temp-filsreferens i csproj | — |

### v4.2 — "Motorn" (november–december)

| Innehåll | Källa |
|----------|-------|
| SVD + SVD-baserad PCA, pseudoinvers, konditionstal | LinearAlgebra Phase 3 |
| Glesa lösare (CG + preconditioners) kopplade till FEM | LinearAlgebra Phase 4 |
| Reverse mode-autodiff (tape) + `GradientTape`-API | Autodiff Phase 3 |
| Allokeringsfria träningsloopar | Performance Phase 3 |
| Gradient boosting (regressor + classifier) | MLExpansion Phase 1–2 |

### v5.0 — "Flaggskeppet" (runt årsskiftet)

| Innehåll | Källa |
|----------|-------|
| Neuronnät på autodiff (nya lager kräver bara forward-kod) | Autodiff Phase 4 |
| **PINN-showcase**: värmeledning validerad mot bibliotekets egen FiniteDifference-lösare | Autodiff Phase 5 |
| Transformer-lager (attention, LayerNorm, encoderblock) | MLExpansion Phase 3 |
| ONNX-export (`CSharpNumerics.Onnx`-paket) | MLExpansion Phase 4 |
| Generic math / `float`-träning (om spiken faller väl ut) | Performance Phase 4 |

---

## Vägledande principer

1. **Mät före optimering** — benchmarks skrivs innan SIMD-arbetet börjar, annars kan förbättringen inte bevisas.
2. **Inga breaking changes utan major-bump** — dekompositioner och kärnor införs bakom befintliga API:er; nya förmågor är opt-in.
3. **Kärnbiblioteket förblir beroendefritt** — externa beroenden (protobuf för ONNX, ONNX Runtime för tester, Math.NET för benchmarks) isoleras i separata paket/projekt.
4. **Befintlig kod som facit** — autodiff verifieras mot `DerivativeExtensions`, PINN mot `FiniteDifference`, histogram-boosting mot exakt-greedy. Biblioteket testar sig självt.
5. **Sekvensering framför parallellitet** — LinearAlgebra före SVD-PCA, autodiff före Transformer-backward, benchmarks före allt prestandaarbete.

---

## Beroendekarta

```
Benchmarks (P1) ──────────────► SIMD-kärnor (P2) ──► Allokeringsfria loopar (P3)
                                      │
LU/Cholesky/QR (LA1–2) ──► SVD (LA3) ─┼──► SVD-PCA
        │                             │
        ├──► Rotfinnare/polynomrötter │
        └──► Glesa lösare (LA4) ──► FEM i skala
                                      │
Dual numbers (AD1) ──► Tape (AD3) ────┴──► NN på autodiff (AD4) ──► PINN (AD5)
                                                    │
Gradient boosting (ML1–2)  [oberoende]              └──► Transformer (ML3)
ONNX (ML4)                 [oberoende]
```
