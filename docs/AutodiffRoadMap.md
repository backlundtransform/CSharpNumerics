# Automatisk differentiering — Dual numbers, reverse mode & PINN

## Mål

Införa **automatisk differentiering (autodiff)** som förstklassig modul: forward mode via dual numbers och reverse mode via tape. Långsiktigt blir reverse mode motorn bakom neuronnätens backprop (ersätter handskriven lager-backprop), och kombinationen autodiff + fysikmodulen möjliggör **Physics-Informed Neural Networks (PINNs)** — en nisch där CSharpNumerics kan bli bäst i .NET-världen, eftersom inget annat bibliotek har både fysiken och ML:en i samma paket.

---

## Nulägesanalys — Befintlig arkitektur

### Vad som finns idag

| Komponent | Status | Plats | Kommentar |
|-----------|--------|-------|-----------|
| Finita differenser | ✓ | `Numerics/DerivativeExtensions.cs` | Derivator av godtycklig ordning, multivariat, komplex |
| Handskriven backprop | ✓ | `ML/NeuralNetwork/`, `ML/Sequence/` | Varje lager implementerar sin egen `Backward` |
| Gradientbaserade optimerare | ✓ | `Numerics/Optimization/` | Adam, GradientDescent — tar färdiga gradienter |
| Physics-informed losses | ✓ | `ML/Losses/` | Conservation, Smoothness, NonNegativity, ConstrainedTrainer |
| `ComplexNumber`, `Quaternion` | ✓ | `Numerics/Objects/` | Mönster för nya taltyper finns etablerat |
| Dual numbers | ✗ | — | Finns ej |
| Tape/beräkningsgraf | ✗ | — | Finns ej |
| Autodiff-gradienter till optimerare | ✗ | — | Optimerarna matas med finita differenser eller handskrivna gradienter |

### Nyckelidentifierade begränsningar

1. **Finita differenser är O(n) funktionsevalueringar per gradient** och lider av trunkerings-/avrundningsfel — dåligt för optimering i hög dimension.
2. **Handskriven backprop skalar inte** — varje nytt lager (attention, nya aktiveringar, custom losses) kräver manuellt härledda och implementerade gradienter; en felräkning ger tyst fel modell.
3. **Physics-informed losses begränsas** — `ConservationLoss` m.fl. behöver derivator av nätverkets output m.a.p. *input* (inte vikter), vilket dagens backprop inte exponerar. Det är exakt vad PINNs kräver.

---

## Del 1 — Forward mode: Dual numbers

### Vad det gör

Ett dualt tal $a + b\varepsilon$ (där $\varepsilon^2 = 0$) propagerar värde och derivata samtidigt genom varje operation:

$$f(a + b\varepsilon) = f(a) + f'(a)\,b\,\varepsilon$$

Kör funktionen med `Dual`-argument → få exakt derivata (maskinprecision), ingen steglängd, ingen trunkering.

### Vad som krävs att bygga

| Komponent | Beskrivning | Insats |
|-----------|-------------|--------|
| **`Dual` (struct)** | `(double Value, double Derivative)`, operatorer `+ − × ÷`, jämförelser | **Låg** — samma mönster som `Quaternion` |
| **`DualMath`** | Sin, Cos, Exp, Log, Pow, Sqrt, Tanh, … med kedjeregel | **Låg** |
| **`HyperDual`** | Andraderivator exakt (för Hessianer/krökning) | **Låg-medel** |
| **`DualVector`** | Multivariat forward mode: gradient av `Func<DualVector, Dual>` i n pass | **Låg** |
| **Extension-fasad** | `Func<double,double>.AutoDerivative(x)` som spegel av befintliga `Derivative`-extensions | **Låg** |

### Direkta vinster i befintlig kod

| Befintlig funktion | Förbättring |
|--------------------|-------------|
| `NewtonRaphson` | Exakt `f'` i stället för finita differenser → färre iterationer, robustare |
| `Optimization/Minimizer` | Exakta gradienter för `Func`-baserade mål |
| `NonlinearFitting` | Exakta Jacobianer i kurvanpassning |
| `KeplerOrbit`, GR-tester | Exakta derivator i fysikberäkningar |

---

## Del 2 — Reverse mode: Tape

### Vad det gör

Reverse mode bygger en beräkningsgraf (tape) under forward-passet och propagerar sedan adjoints bakåt — **hela gradienten till kostnaden av ~2 funktionsevalueringar**, oavsett antal parametrar. Det är detta som är "backpropagation", generaliserad till godtyckliga program.

### Vad som krävs att bygga

| Komponent | Beskrivning | Insats |
|-----------|-------------|--------|
| **`AdVariable`** | Nod med värde + adjoint + referens till tape-position | **Låg** |
| **`Tape`** | Append-only lista av operationer `(op, parents, lokala derivator)`. `Backward()` itererar baklänges | **Medel** |
| **Operator-overloads + `AdMath`** | Samma yta som `DualMath`, men registrerar på tape | **Medel** |
| **Vektor-/matrisoperationer på tape** | MatMul, Dot, Hadamard, Sum, aktiveringar som *enskilda tape-noder* (inte elementvis — annars exploderar tapen) | **Medel-hög** |
| **`GradientTape`-API** | `using var tape = new GradientTape(); … var grads = tape.Gradient(loss, parameters);` — bekant mönster från TF/PyTorch | **Låg** |
| **Checkpointing/tape-reset** | Minneshantering för långa träningsloopar | **Medel** |

### Designprinciper

- Placeras i `Numerics/AutoDiff/` — det är en numerik-primitiv, inte en ML-feature.
- Tape-operationerna delegerar till samma `Kernels` som SIMD-arbetet (se PerformanceRoadMap) — autodiff får vektorisering gratis.
- **Verifiering:** varje tape-operation testas mot finita differenser (`DerivativeExtensions` blir facit-generator — befintlig kod som testinfrastruktur).

---

## Del 3 — Neuronnät på autodiff + PINN

### Steg 1: NN-migrering (opt-in)

`SequentialModel` får en autodiff-driven träningsväg vid sidan av den handskrivna: lagren definierar bara `Forward` i tape-operationer, gradienterna kommer från `tape.Gradient()`. Handskrivna `Backward` behålls tills paritet är bevisad (samma gradienter, jämförbar hastighet), sedan fasas de ut.

Vinst: nya lager (attention! — se MLExpansionRoadMap) kräver bara forward-kod. Custom losses blir triviala.

### Steg 2: PINN-showcase

En PINN tränar ett nät $u_\theta(x, t)$ att uppfylla en PDE genom att lägga PDE-residualen i lossen:

$$\mathcal{L} = \underbrace{\|u_\theta - u_{data}\|^2}_{\text{data}} + \lambda \underbrace{\|\partial_t u_\theta + \mathcal{N}[u_\theta]\|^2}_{\text{fysik (kräver } \partial u/\partial x \text{ via autodiff)}}$$

| Komponent | Beskrivning | Beroenden |
|-----------|-------------|-----------|
| **Input-gradienter** | `tape.Gradient(output, inputs)` — derivator m.a.p. *input*, inte bara vikter | Del 2 |
| **`PinnTrainer`** | Kollokationspunkter + data-loss + PDE-residual-loss, bygger på `ConstrainedTrainer`-mönstret | `ML/Losses/` |
| **Showcase 1: värmeledning 1D** | Jämför PINN mot befintlig `FiniteDifference`-lösning — facit finns redan i biblioteket! | `FiniteDifference/` |
| **Showcase 2: harmonisk oscillator** | ODE-PINN mot `Mechanics/Oscillators` | `Physics/` |
| **Showcase 3 (sträckmål): Burgers ekvation** | Klassiskt PINN-benchmark från litteraturen (Raissi et al. 2019) | Fluid-modulen som referens |

Att facit-lösarna redan finns i samma bibliotek är unikt — valideringen blir en testfil, inte ett forskningsprojekt.

---

## Implementationsplan — Faser

### Phase 1 — Dual numbers
- [ ] Implementera `Dual` struct + `DualMath` i `Numerics/AutoDiff/`
- [ ] `DualVector` för multivariata gradienter
- [ ] `AutoDerivative`-extensions som spegel av `DerivativeExtensions`
- [ ] Koppla in i `NewtonRaphson` och `Minimizer` (opt-in-överlagringar)
- [ ] Enhetstester: alla `DualMath`-funktioner mot analytiska derivator

### Phase 2 — HyperDual & Jacobianer
- [ ] `HyperDual` för exakta andraderivator
- [ ] Jacobian-/Hessian-hjälpare för `Func<VectorN, VectorN>`
- [ ] Koppla in i `NonlinearFitting`

### Phase 3 — Reverse mode-kärna
- [ ] `Tape`, `AdVariable`, skalära operationer + `AdMath`
- [ ] `GradientTape`-API
- [ ] Testsvit: varje operation verifieras mot finita differenser
- [ ] Vektor-/matrisoperationer som tape-noder (MatMul, Dot, Hadamard, Sum, aktiveringar)

### Phase 4 — NN på autodiff
- [ ] Autodiff-träningsväg i `SequentialModel` (opt-in)
- [ ] Paritetstester: autodiff-gradienter ≡ handskrivna för Dense, Conv1D, LSTM
- [ ] Benchmark: autodiff vs handskriven backprop (målsättning: inom 1.5×)
- [ ] Beslut om utfasning av handskrivna `Backward`

### Phase 5 — PINN
- [ ] Input-gradienter genom tapen
- [ ] `PinnTrainer` + kollokationspunkt-sampling
- [ ] Showcase: värmeledning 1D validerad mot `FiniteDifference`
- [ ] Showcase: harmonisk oscillator
- [ ] Dokumentation + exempel i README (flaggskepps-feature!)

---

## Sammanfattning

| Del | Genomförbarhet | Insats | Största risk |
|-----|---------------|--------|--------------|
| **Dual numbers** | **Hög** | Låg | Ingen — perfekt höststart |
| **Reverse mode** | **Medel-hög** | Hög | Tape-prestanda; elementvisa noder får inte explodera |
| **NN-migrering** | **Medel** | Medel | Gradientparitet med befintliga lager måste bevisas |
| **PINN** | **Medel-hög** (givet Del 2) | Medel | Träningsstabilitet — välkänt känsligt, men facit-lösare finns |

**Rekommendation:** Phase 1 är liten, elegant och ger omedelbar nytta — börja där. Reverse mode är höstens tyngsta men mest strategiska bygge; designa tape-API:t tidigt och verifiera mot finita differenser hela vägen. PINN-showcasen är det som syns utåt — den motiverar allt annat.
