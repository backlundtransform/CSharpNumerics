# Prestanda — SIMD, benchmarks & generic math

## Mål

Göra CSharpNumerics **mätbart snabbt** och kunna bevisa det: ett BenchmarkDotNet-projekt med publicerade jämförelser mot Math.NET Numerics, SIMD-vektorisering av kärnoperationerna, och `INumber<T>`-baserade generiska typer så att ML-delen kan träna i `float`. Prestanda är trovärdighetsfrågan — ledande ramverk visar siffror.

---

## Nulägesanalys — Befintlig arkitektur

### Vad som finns idag

| Komponent | Status | Kommentar |
|-----------|--------|-----------|
| Multi-targeting `net10.0; net8.0; netstandard2.1` | ✓ | `CSharpNumerics.csproj` — bra grund för villkorlig SIMD |
| `Matrix`/`VectorN`/`Tensor`-operationer | ✓ | Skalära loopar, `double[,]`/`double[]` |
| Benchmark-projekt | ✗ | Finns ej — inga mätningar alls |
| SIMD (`Vector<T>`, `TensorPrimitives`) | ✗ | Används ej |
| `Span<T>`-API:er | ✗ | Alla operationer allokerar nya objekt |
| `float`-stöd i ML | ✗ | Allt är `double` — dubbelt minne, halv throughput för NN-träning |
| Parallellisering | ✗ | Ingen `Parallel.For` i matmul/faltning/träningsloopar |

### Nyckelidentifierade begränsningar

1. **Matrismultiplikation är naiv trippel-loop** — ingen cache-blocking, ingen SIMD, ingen parallellism. Detta är den heta loopen i NN-träning, FEM och dekompositioner.
2. **`VectorN`-operationer allokerar** — varje `+`, `Hadamard`, `Dot` skapar nya arrayer; i träningsloopar med miljontals iterationer dominerar GC-trycket.
3. **`netstandard2.1`-target saknar `TensorPrimitives` och `INumber<T>`** — kräver `#if`-strategi eller separat kodväg.
4. **Inga regressionstester för prestanda** — en långsam ändring märks aldrig.

---

## Del 1 — BenchmarkDotNet-projekt

### Vad som krävs att bygga

| Komponent | Beskrivning | Insats |
|-----------|-------------|--------|
| **`Numerics.Benchmarks`-projekt** | Nytt konsol-projekt i solutionen, `net10.0`, BenchmarkDotNet | **Låg** |
| **Kärnbenchmarks** | MatMul (64/256/1024), dekompositioner, FFT, Dot/Hadamard, ODE-steg | **Låg** |
| **ML-benchmarks** | MLP-epoch, Conv1D forward/backward, KMeans-iteration, RandomForest-fit | **Låg-medel** |
| **Konkurrent-jämförelser** | Samma operationer i Math.NET Numerics (referens i benchmark-projektet, ej i biblioteket) | **Låg** |
| **CI-integration** | GitHub Action som kör benchmarks på PR-label och kommenterar resultat | **Medel** |
| **README-sektion** | Publicerade tabeller + grafer, uppdateras per release | **Låg** |

Benchmarks skrivs **före** optimeringarna i Del 2 — annars går det inte att visa förbättringen.

---

## Del 2 — SIMD & allokeringsfria kärnor

### Strategi

Inför ett internt lager `Numerics/LinearAlgebra/Kernels/` med statiska lågnivåoperationer på `Span<double>`/`ReadOnlySpan<double>`. Publika typer (`Matrix`, `VectorN`, `Tensor`) behåller sina API:er men delegerar till kärnorna — inga breaking changes.

| Komponent | Beskrivning | Insats |
|-----------|-------------|--------|
| **`VectorKernels`** | Add, Subtract, Multiply, Dot, Hadamard, Norm, Scale på `Span<double>` via `TensorPrimitives` (net8+) med skalär fallback (netstandard2.1) | **Låg-medel** |
| **`MatMulKernel`** | Cache-blockad matmul + `Vector<double>`-inre loop + `Parallel.For` över radblock för stora matriser | **Medel-hög** |
| **`#if NET8_0_OR_GREATER`-strategi** | SIMD-vägen för moderna targets, skalär för netstandard2.1 — samma tester körs mot båda | **Låg** |
| **In-place-överlagringar** | `VectorN.AddInPlace`, `Matrix.MultiplyInto(result)` för heta träningsloopar | **Medel** |
| **Träningsloop-genomgång** | Ersätt allokerande operatorer med in-place-kärnor i `NeuralNetwork`, `SequentialModel`, LSTM/Conv-lager | **Medel** |

### Förväntad effekt

| Operation | Förväntad speedup | Mekanism |
|-----------|-------------------|----------|
| Dot/Hadamard/Add | 4–8× | AVX2/AVX-512 via TensorPrimitives |
| MatMul (1024×1024) | 10–50× | Blocking + SIMD + parallellism |
| NN-träning (epoch) | 3–10× | Ovanstående + eliminerad GC-press |

---

## Del 3 — Generic math (`INumber<T>`)

### Vad som krävs att bygga

| Komponent | Beskrivning | Insats |
|-----------|-------------|--------|
| **`Matrix<T>` / `VectorN<T>` / `Tensor<T>`** | Generiska varianter med `where T : INumber<T>` (net8+ only). Befintliga `Matrix` etc. blir alias/wrapper för `T = double` | **Hög** |
| **`float`-träningsväg i ML** | `SequentialModel` och lager parametriserade på `T` → halverat minne, dubblad SIMD-bredd | **Medel-hög** |
| **Precision-tester** | Samma testsvit körs för `double` och `float` med anpassade toleranser | **Medel** |

### Avgränsning

Generiska typer görs **endast för net8.0/net10.0** — netstandard2.1 behåller `double`-typerna. Det är ett medvetet val: konsumenter på moderna targets får full funktionalitet, äldre targets fryses funktionsmässigt.

**Alternativ att utreda först:** i stället för fullt generiska publika typer kan enbart ML-delens interna beräkningsväg parametriseras (`float`-vikter i lager). Mindre yta, snabbare vinst — rekommenderad startpunkt.

---

## Implementationsplan — Faser

### Phase 1 — Mät först
- [ ] Skapa `Numerics.Benchmarks`-projekt med BenchmarkDotNet
- [ ] Kärnbenchmarks: MatMul, Dot, Hadamard, FFT, Simpson, RK4
- [ ] ML-benchmarks: MLP-epoch, Conv1D, KMeans
- [ ] Jämförelsebenchmarks mot Math.NET Numerics
- [ ] Baseline-resultat sparas i `docs/benchmarks/baseline.md`

### Phase 2 — Vektoriserade kärnor
- [ ] `Numerics/LinearAlgebra/Kernels/` med `VectorKernels` (TensorPrimitives + fallback)
- [ ] Delegera `VectorN`-operationer till kärnorna
- [ ] Cache-blockad + parallell `MatMulKernel`, delegera `Matrix`-multiplikation
- [ ] Verifiera identiska resultat mot skalära vägen (bitvis-tolerans-tester)
- [ ] Kör om benchmarks — dokumentera speedup

### Phase 3 — Allokeringsfria träningsloopar
- [ ] In-place-API:er (`AddInPlace`, `MultiplyInto`)
- [ ] Genomgång av `NeuralNetwork`/`SequentialModel`/sekvenslager: eliminera allokeringar i inner loops
- [ ] `MemoryDiagnoser`-benchmarks före/efter

### Phase 4 — Generic math (utredning → implementation)
- [ ] Spike: `float`-parametriserad `DenseLayer` — mät faktisk vinst
- [ ] Beslut: full `Matrix<T>` eller enbart ML-intern `T`-väg
- [ ] Implementera vald väg för net8+
- [ ] Precision-testsvit för `float`

### Phase 5 — CI & publicering
- [ ] Benchmark-workflow i GitHub Actions (label-triggad)
- [ ] README-sektion med resultattabeller
- [ ] Regressionströskel: PR varnas om kärnbenchmark försämras >10 %

---

## Sammanfattning

| Del | Genomförbarhet | Insats | Största risk |
|-----|---------------|--------|--------------|
| **Benchmarks** | **Hög** | Låg | Ingen — ren infrastruktur |
| **SIMD-kärnor** | **Hög** | Medel | netstandard2.1-fallback dubblar testytan |
| **Generic math** | **Medel** | Hög | API-yta exploderar om allt generifieras — börja smalt |

**Rekommendation:** Phase 1 först, alltid. Sedan Phase 2 som ger störst vinst per timme. Generic math (Phase 4) startas som en avgränsad spike innan beslut om full utrullning — det är lätt att gräva ner sig där.
