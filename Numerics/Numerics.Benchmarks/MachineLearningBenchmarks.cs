using BenchmarkDotNet.Attributes;
using CSharpNumerics.ML.Models.Regression;
using CSharpNumerics.Numerics.Objects;

namespace Numerics.Benchmarks;

/// <summary>
/// Baseline for the neural network training loop — the reference point for the
/// allocation-free training work planned for v4.4.
/// </summary>
/// <remarks>
/// One epoch over a fixed synthetic regression set. The allocation column matters more
/// here than the timing: every forward and backward pass currently allocates new
/// <see cref="Matrix"/> and <see cref="VectorN"/> instances.
/// </remarks>
[MemoryDiagnoser]
public class MachineLearningBenchmarks
{
    [Params(256, 2048)]
    public int Samples;

    private const int Features = 16;

    private Matrix inputs = default!;
    private VectorN targets = default!;

    [GlobalSetup]
    public void Setup()
    {
        var random = new Random(20260930);

        var x = new double[Samples, Features];
        var y = new double[Samples];

        for (var i = 0; i < Samples; i++)
        {
            var sum = 0.0;
            for (var j = 0; j < Features; j++)
            {
                var value = random.NextDouble() * 2.0 - 1.0;
                x[i, j] = value;
                sum += value * (j + 1);
            }

            y[i] = sum / Features;
        }

        inputs = new Matrix(x);
        targets = new VectorN(y);
    }

    /// <summary>
    /// A single training epoch of a 16 → 32 → 16 → 1 network. Early stopping is held off
    /// so the measurement is one full pass rather than a variable number of them.
    /// </summary>
    [Benchmark(Description = "MLP regressor, one training epoch")]
    public MLPRegressor TrainOneEpoch()
    {
        var model = new MLPRegressor
        {
            HiddenLayers = new[] { 32, 16 },
            Epochs = 1,
            BatchSize = 32,
            LearningRate = 0.01,
            Patience = int.MaxValue
        };

        model.Fit(inputs, targets);

        return model;
    }
}
