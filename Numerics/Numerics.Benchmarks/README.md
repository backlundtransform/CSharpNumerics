# Numerics.Benchmarks

Measurement project for CSharpNumerics, built on [BenchmarkDotNet](https://benchmarkdotnet.org/).

It exists so that performance claims can be checked rather than asserted. The guiding
principle in the roadmap is *measure before optimising*: the SIMD and allocation-free work
planned for v4.4 is judged against the baseline recorded here, and the v4.3 migration of
call sites onto the shared decompositions is expected to leave these numbers no worse.

## Why it is a separate project

The core library has no external dependencies, and it stays that way. BenchmarkDotNet lives
here and nowhere else. The project is `IsPackable=false`, so it never reaches the NuGet
package, and CI runs tests against `NumericTest` directly, so benchmarks never run there —
they are far too slow for per-commit execution.

## Running

Release configuration is required; BenchmarkDotNet refuses to measure a Debug build.

```powershell
# Pick from a menu
dotnet run -c Release --project Numerics/Numerics.Benchmarks

# Everything (the full baseline — takes a while)
dotnet run -c Release --project Numerics/Numerics.Benchmarks -- --filter *

# One group
dotnet run -c Release --project Numerics/Numerics.Benchmarks -- --filter *Dense*
dotnet run -c Release --project Numerics/Numerics.Benchmarks -- --filter *Sparse*
dotnet run -c Release --project Numerics/Numerics.Benchmarks -- --filter *MachineLearning*

# List without running
dotnet run -c Release --project Numerics/Numerics.Benchmarks -- --list flat

# Smoke test that every benchmark executes (one cold iteration each, timings meaningless)
dotnet run -c Release --project Numerics/Numerics.Benchmarks -- --job Dry --filter *
```

Reports land in `docs/benchmarks/artifacts/` (git-ignored). To record a baseline, copy the
generated markdown report into `docs/benchmarks/` under a name that says what it measured and
at which commit.

## What is measured

| Group | Benchmarks | Why |
|-------|-----------|-----|
| `DenseLinearAlgebraBenchmarks` | Matrix×Matrix, Matrix×Vector, LU / QR / Cholesky factorization, LU factor-and-solve, LU solve reusing a factorization, `Matrix.Inverse` | The triple-loop multiply is the hot path in neural network training, FEM and every decomposition. The two LU solve variants show what reusing a factorization is worth. |
| `SparseLinearAlgebraBenchmarks` | SpMV, PCG solve | The 5-point Laplacian on an interior grid — the shape a 2-D FEM or finite difference problem assembles. SpMV is the inner loop of every iterative solver. |
| `MachineLearningBenchmarks` | One MLP training epoch | Reference point for the allocation-free training loops planned for v4.4. |

Every benchmark carries `[MemoryDiagnoser]`. For this library the **Allocated** column is often
more informative than the timing: the current kernels allocate a fresh `Matrix` or `VectorN`
for every intermediate result, and that is precisely what later work aims to remove.

## Writing a benchmark

- Seed any randomness with a fixed constant in `[GlobalSetup]` so runs are comparable across
  commits.
- Return the result from the benchmark method. A method that returns `void` and computes into
  a local can be optimised away entirely.
- Keep setup out of the measured method unless the setup *is* what you are measuring — the
  two LU benchmarks are a deliberate pair that separates factoring from solving.
- Build the problem correctly. A conjugate gradient benchmark on a matrix that is accidentally
  non-symmetric measures a solver that has no reason to converge, not the solver's real cost.
