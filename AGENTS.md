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

- Each section has its own `README.md` inside its folder. **Update the section README after completing work** that adds or changes public API in that section.
- The **root `README.md`** is marketing-facing and should only be updated by a developer — never auto-update it.
- Section READMEs follow an example-driven style with code samples, tables, and architecture diagrams. See `Engines/GIS/README.md` or `ML/README.md` for reference.

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
