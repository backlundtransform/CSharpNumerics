using CSharpNumerics.ML.Enums;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Linq;

namespace CSharpNumerics.ML.NeuralNetwork;

public static class Activations
{
    public static VectorN Apply(VectorN vector, ActivationType activation)
    {
        return activation switch
        {
            ActivationType.ReLU => ApplyElementwise(vector, x => Math.Max(0.0, x)),
            ActivationType.Sigmoid => ApplyElementwise(vector, x => 1.0 / (1.0 + Math.Exp(-x))),
            ActivationType.Tanh => ApplyElementwise(vector, Math.Tanh),
            ActivationType.Linear => new VectorN(vector.Values),
            _ => throw new ArgumentOutOfRangeException(nameof(activation))
        };
    }

    public static VectorN Derivative(VectorN activatedVector, ActivationType activation)
    {
        return activation switch
        {
            ActivationType.ReLU => ApplyElementwise(activatedVector, x => x > 0.0 ? 1.0 : 0.0),
            ActivationType.Sigmoid => ApplyElementwise(activatedVector, x => x * (1.0 - x)),
            ActivationType.Tanh => ApplyElementwise(activatedVector, x => 1.0 - (x * x)),
            ActivationType.Linear => ApplyElementwise(activatedVector, _ => 1.0),
            _ => throw new ArgumentOutOfRangeException(nameof(activation))
        };
    }

    public static VectorN Softmax(VectorN vector)
    {
        double max = vector.Values.Max();
        var exp = new double[vector.Length];
        double sum = 0.0;

        for (int i = 0; i < vector.Length; i++)
        {
            exp[i] = Math.Exp(vector[i] - max);
            sum += exp[i];
        }

        for (int i = 0; i < exp.Length; i++)
        {
            exp[i] /= sum;
        }

        return new VectorN(exp);
    }

    private static VectorN ApplyElementwise(VectorN vector, Func<double, double> transform)
    {
        var result = new double[vector.Length];
        for (int i = 0; i < vector.Length; i++)
        {
            result[i] = transform(vector[i]);
        }

        return new VectorN(result);
    }
}