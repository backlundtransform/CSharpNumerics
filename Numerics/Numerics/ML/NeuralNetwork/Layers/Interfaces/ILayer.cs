using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;

namespace CSharpNumerics.ML.NeuralNetwork.Layers.Interfaces;

public interface ILayer
{
    VectorN[] Forward(VectorN[] input, bool training = true);

    VectorN[] Backward(VectorN[] gradOutput);

    void ApplyGradients(IOptimizer weightOptimizer, IOptimizer biasOptimizer, int batchSize);

    int ParameterCount { get; }
}