# Project Guidelines

## Architecture

CSharpNumerics is a pure C# scientific computing library organised into **five sections**. Each section maps to a top-level namespace under `CSharpNumerics`.

| Section | Namespace | Responsibility |
|---------|-----------|----------------|
| **Numerical Analysis** | `CSharpNumerics.Numerics` | Pure mathematical abstractions — linear algebra, transforms, ODEs, vector fields, numerical integration, optimization, complex numbers, interpolation, finite differences. If it is a mathematical object or algorithm with no domain-specific semantics, it belongs here. |
| **Statistics** | `CSharpNumerics.Statistics` | Stochastic analysis, probability theory, statistical methods, distributions, hypothesis testing, Monte Carlo simulation, and data processing such as time series. |
| **Machine Learning** | `CSharpNumerics.ML` | Supervised and unsupervised learning, neural networks, reinforcement learning, cross-validation, pipelines, scalers, feature selection, and dimensionality reduction. |
| **Physics** | `CSharpNumerics.Physics` | Physics models and domain-specific equations — heat transfer, electromagnetics, environmental modelling, orbital mechanics, oscillations, radiation, astronomy. A vector field is a mathematical object (→ Numerics); the heat equation is physics (→ Physics). |
| **Simulation Engines** | `CSharpNumerics.Engines` | Simulator-level logic built on top of the other sections. Sub-engines: `Audio`, `Game`, `GIS`. Shared engine infrastructure (event bus, clocks, adapters) lives in `Engines.Common`. |

### Cross-section dependencies

Sections may depend on each other. Typical flow:

```
Engines  →  ML / Physics / Statistics / Numerics
ML       →  Numerics / Statistics
Physics  →  Numerics
Statistics → Numerics
Numerics →  (no internal deps)
```

### Boundary rules

- **Numerics** owns all math primitives: `Matrix`, `VectorN`, `Vector`, `ComplexNumber`, `ScalarField`, `VectorField`, transforms, ODE solvers, optimization (`Optimization/`).
- **Statistics** owns stochastic methods, distributions, Monte Carlo, time series, hypothesis tests.
- **ML** owns learning algorithms and pipelines. It may consume `Numerics.Optimization` but must not reimplement optimizers.
- **Physics** owns physical models and constants. It uses `Numerics` math objects but adds physical semantics.
- **Engines** orchestrate simulations. `Engines.Common` holds shared infrastructure (interfaces, event bus, clock) that may be used by any sub-engine. Engine-specific logic stays in its own sub-folder (`Audio/`, `Game/`, `GIS/`).

## Build and Test

```powershell
# Build all targets (net10.0, net8.0, netstandard2.1)
dotnet build Numerics/Numerics.sln

# Run all tests
dotnet test Numerics/Numerics.sln

# Run tests for a specific section
dotnet test --filter "FullyQualifiedName~MachineLearningTests"
dotnet test --filter "FullyQualifiedName~OptimizationTests"
```

## Documentation

- The library has **five sections**: Numerics, Statistics, ML, Physics, Engines. Each section has its own `README.md` inside its top-level folder.
- **Sub-sections belong to their parent section's README.** For example, `Optimization/` is a sub-section of Numerics — its documentation goes in `Numerics/README.md`, not in a separate `Optimization/README.md`. Do **not** create per-sub-section READMEs.
- **Update the section README after completing work** that adds or changes public API in that section.
- The **root `README.md`** is marketing-facing and should only be updated by a developer — never auto-update it.
- Section READMEs follow an example-driven style with code samples, tables, and architecture diagrams. See `Engines/GIS/README.md` or `ML/README.md` for reference.

### Section README locations

| Section | README path |
|---------|-------------|
| Numerics | `Numerics/Numerics/Numerics/README.md` |
| Statistics | `Numerics/Numerics/Statistics/README.md` |
| ML | `Numerics/Numerics/ML/README.md` |
| Physics | `Numerics/Numerics/Physics/README.md` |
| Engines | `Numerics/Numerics/Engines/README.md` (sub-engines may have their own, e.g. `Engines/GIS/README.md`) |

## Feature Roadmaps

Feature work follows a **research → plan → implement → complete** lifecycle tracked via roadmap files in `docs/`.

### Workflow

1. **Research phase** — investigate the feature, gather context, and produce a `docs/<FeatureName>RoadMap.md`.
2. **Planning phase** — divide the roadmap into numbered **phases** with checkable tasks (`- [ ]`).
3. **Implementation phases** — work through each phase, checking off tasks (`- [x]`) as they are completed.
4. **Completion** — when every task is checked, move the file to `docs/completed/`.

### Rules

- File naming: `docs/<FeatureName>RoadMap.md` (PascalCase feature name).
- Multiple roadmaps may be active in `docs/` at the same time.
- Each roadmap must be divided into **phases** (e.g., Phase 1, Phase 2, …). Phases group related tasks and are implemented in order.
- Use Markdown checkboxes (`- [ ]` / `- [x]`) inside each phase to track individual tasks.
- A roadmap is considered **complete** only when all checkboxes across all phases are checked.
- Completed roadmaps are moved to `docs/completed/` — do not delete them.

### Example structure

```
docs/
  QuaternionAlgebraRoadMap.md   ← active
  AudioDSPRoadMap.md            ← active
  completed/
    MonteCarloRoadMap.md        ← finished
```

### Example roadmap template

```markdown
# <Feature Name> Roadmap

## Phase 1 — Research & Design
- [ ] Investigate existing API surface
- [ ] Define public types and interfaces
- [ ] Draft architecture

## Phase 2 — Core Implementation
- [ ] Implement primary classes
- [ ] Add unit tests

## Phase 3 — Integration & Docs
- [ ] Wire into consuming sections
- [ ] Update section README
```

## Code Style

- File-scoped namespaces (`namespace X;`).
- Target: `net10.0;net8.0;netstandard2.1` — code must compile on all three.
- No external dependencies in the core library (`CSharpNumerics.csproj`). All algorithms are implemented from scratch.
- Interfaces go in an `Interfaces/` sub-folder within their section.
- Enums go in an `Enums/` sub-folder within their section.
- Tests use MSTest (`[TestClass]`, `[TestMethod]`) in the `NumericTest` project.

## Conventions

- **Namespace matches folder path**: `CSharpNumerics.Numerics.Optimization.SingleObjective` → `Numerics/Optimization/SingleObjective/`.
- **One public class per file**, file name matches class name.
- Prefer `VectorN` and `Matrix` from `CSharpNumerics.Numerics.Objects` for all linear algebra.
- When adding a new capability that could be shared (e.g., an optimizer or convergence criterion), place it in the most general section that owns the abstraction, not in the consuming section.
