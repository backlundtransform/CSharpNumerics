using System.Reflection;
using BenchmarkDotNet.Running;

// Entry point for the benchmark suite.
//
//   dotnet run -c Release --project Numerics/Numerics.Benchmarks            (menu)
//   dotnet run -c Release --project Numerics/Numerics.Benchmarks -- --filter *Dense*
//   dotnet run -c Release --project Numerics/Numerics.Benchmarks -- --list flat
//
// Results are written to BenchmarkDotNet.Artifacts/; copy the markdown report into
// docs/benchmarks/ when recording a baseline.
BenchmarkSwitcher.FromAssembly(Assembly.GetExecutingAssembly()).Run(args);
