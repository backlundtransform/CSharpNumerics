using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Objects;
using Numerics;
using Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Text;

namespace CSharpNumerics.ML.Models.Regression;

public class MLPRegressor : IModel,IHasHyperparameters, IRegressionModel
{
    public int[] HiddenLayers { get; set; } = new[] { 32, 16 };
    public double LearningRate { get; set; } = 0.01;
    public int Epochs { get; set; } = 500;
    public int BatchSize { get; set; } = 16;
    public double L2 { get; set; } = 0.0;

    public ActivationType Activation { get; set; } = ActivationType.ReLU;


    private List<Matrix> _weights;
    private List<VectorN> _biases;

    private void Initialize(int inputSize)
    {
        _weights = new();
        _biases = new();

        var layers = new List<int> { inputSize };
        layers.AddRange(HiddenLayers);
        layers.Add(1); // output

        var rnd = new Random(123);

        for (int i = 0; i < layers.Count - 1; i++)
        {
            _weights.Add(RandomMatrix(layers[i], layers[i + 1], rnd));
            _biases.Add(new VectorN(layers[i + 1]));
        }
    }

    private VectorN Forward(VectorN x, out List<VectorN> activations)
    {
        activations = new() { x };

        for (int i = 0; i < _weights.Count; i++)
        {
            var z = (_weights[i].Transpose() * activations[^1]) + _biases[i];
            var a = (i == _weights.Count - 1)
                ? z
                : Activate(z);
            activations.Add(a);
        }

        return activations[^1];
    }

    private void Backward(
    VectorN y,
    List<VectorN> activations)
    {
        var deltas = new List<VectorN>();

       
        var error = activations[^1] - y;
        deltas.Add(error);

    
        for (int i = _weights.Count - 2; i >= 0; i--)
        {
            var w = _weights[i + 1];
            var delta = (w * deltas[^1]).Hadamard(ActivationDerivative(activations[i + 1]));
            deltas.Add(delta);
        }

        deltas.Reverse();

        for (int i = 0; i < _weights.Count; i++)
        {
            _weights[i] -= LearningRate *
                (activations[i].Outer(deltas[i]) + L2 * _weights[i]);

            _biases[i] -= LearningRate * deltas[i];
        }
    }

    public void Fit(Matrix X, VectorN y)
    {
        Initialize(X.columnLength);

        for (int epoch = 0; epoch < Epochs; epoch++)
        {
            for (int i = 0; i < X.rowLength; i++)
            {
                var xi = X.RowSlice(i);
                var yi = new VectorN(new[] { y[i] });

                var yhat = Forward(xi, out var acts);
                Backward(yi, acts);
            }
        }
    }

    public VectorN Predict(Matrix X)
    {
        var result = new double[X.rowLength];

        for (int i = 0; i < X.rowLength; i++)
        {
            var yhat = Forward(X.RowSlice(i), out _);
            result[i] = yhat[0];
        }

        return new VectorN(result);
    }


    public void SetHyperParameters(Dictionary<string, object> p)
    {
        if (p.TryGetValue("HiddenLayers", out var h)) HiddenLayers = (int[])h;
        if (p.TryGetValue("LearningRate", out var lr)) LearningRate = Convert.ToDouble(lr);
        if (p.TryGetValue("Epochs", out var e)) Epochs = Convert.ToInt32(e);
        if (p.TryGetValue("BatchSize", out var b)) BatchSize = Convert.ToInt32(b);
        if (p.TryGetValue("L2", out var l2)) L2 = Convert.ToDouble(l2);
        if (p.TryGetValue("Activation", out var a)) Activation = (ActivationType)a;

    }

    private static Matrix RandomMatrix(int rows, int cols, Random rnd)
    {
        var values = new double[rows, cols];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                values[i, j] = rnd.NextDouble() * 2 - 1; // Uniform [-1, 1]
        return new Matrix { values = values, rowLength = rows, columnLength = cols };
    }


    private VectorN ActivationDerivative(VectorN v)
    {
        var result = new double[v.Length];
        for (int i = 0; i < v.Length; i++)
            result[i] = Derivative(v.Values[i]);
        return new VectorN(result);
    }

    public double Derivative(double x)
    {
        return Activation switch
        {
            ActivationType.ReLU => x > 0 ? 1.0 : 0.0,

           
            ActivationType.Sigmoid => Math.Exp(-x) / Math.Pow(1.0 + Math.Exp(-x), 2),

          
            ActivationType.Tanh => 1.0 - Math.Pow(Math.Abs(Math.Tanh(x)), 2),

            ActivationType.Linear => 1.0,
            _ => 1.0
        };
    }
    private VectorN Activate(VectorN v)
    {
        switch (Activation)
        {
            case ActivationType.ReLU:
                return ReLU(v);
            case ActivationType.Sigmoid:
                return Sigmoid(v);
            case ActivationType.Tanh:
                return Tanh(v);
            case ActivationType.Linear:
                return Linear(v);
            default:
                throw new ArgumentOutOfRangeException();
        }
    }

    private static VectorN ReLU(VectorN v)
    {
        var result = new double[v.Length];
        for (int i = 0; i < v.Length; i++)
            result[i] = Math.Max(0, v.Values[i]);
        return new VectorN(result);
    }
    private static VectorN Sigmoid(VectorN v)
    {
        var result = new double[v.Length];
        for (int i = 0; i < v.Length; i++)
            result[i] = 1.0 / (1.0 + Math.Exp(-v.Values[i]));
        return new VectorN(result);
    }
    private static VectorN Tanh(VectorN v)
    {
        var result = new double[v.Length];
        for (int i = 0; i < v.Length; i++)
            result[i] = Math.Tanh(v.Values[i]);
        return new VectorN(result);
    }
    private static VectorN Linear(VectorN v)
    {
        return v;
    }
}
