using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;
using NN = CSharpNumerics.ML.NeuralNetwork.NeuralNetwork;

namespace CSharpNumerics.ML.Models.Regression;

public class MLPRegressor : IHasHyperparameters, IRegressionModel
{
    public int[] HiddenLayers { get; set; } = new[] { 32, 16 };
    public double LearningRate { get; set; } = 0.01;
    public int Epochs { get; set; } = 500;
    public int BatchSize { get; set; } = 16;
    public double L2 { get; set; } = 0.0;

    public double ValidationSplit { get; set; } = 0.2;
    public int Patience { get; set; } = 10;
    public double MinDelta { get; set; } = 1e-4;

    public ActivationType Activation { get; set; } = ActivationType.ReLU;

    private NN _network;

    public void Fit(Matrix X, VectorN y)
    {
        _network = new NN(X.columnLength, HiddenLayers, 1,
            Activation, NN.OutputMode.Linear);

        int[] indices = Enumerable.Range(0, X.rowLength).ToArray();
        var rnd = new Random(42);

        int trainSize = (ValidationSplit > 0 && X.rowLength > 5)
            ? (int)(X.rowLength * (1 - ValidationSplit))
            : X.rowLength;

        int valSize = X.rowLength - trainSize;

        double bestValLoss = double.MaxValue;
        int patienceCounter = 0;

        var (bestWeights, bestBiases) = _network.SnapshotWeights();

        for (int epoch = 0; epoch < Epochs; epoch++)
        {
            Shuffle(indices, rnd);

            for (int i = 0; i < trainSize; i += BatchSize)
            {
                int currentBatchSize = Math.Min(BatchSize, trainSize - i);
                var weightGrads = _network.InitWeightGrads();
                var biasGrads = _network.InitBiasGrads();

                for (int b = 0; b < currentBatchSize; b++)
                {
                    int idx = indices[i + b];
                    var xi = X.RowSlice(idx);
                    var yi = new VectorN(new[] { y[idx] });

                    var yhat = _network.Forward(xi, out var acts);
                    _network.ComputeGradients(yi, acts, out var dw, out var db);
                    _network.AccumulateGradients(weightGrads, biasGrads, dw, db);
                }

                _network.ApplyGradients(weightGrads, biasGrads, LearningRate, currentBatchSize, L2);
            }

            if (valSize > 0)
            {
                double currentValLoss = CalculateValidationLoss(X, y, indices, trainSize);

                if (currentValLoss < bestValLoss - MinDelta)
                {
                    bestValLoss = currentValLoss;
                    patienceCounter = 0;
                    (bestWeights, bestBiases) = _network.SnapshotWeights();
                }
                else
                {
                    patienceCounter++;
                }

                if (patienceCounter >= Patience)
                {
                    _network.RestoreWeights(bestWeights, bestBiases);
                    return;
                }
            }
        }

        if (valSize > 0)
        {
            _network.RestoreWeights(bestWeights, bestBiases);
        }
    }

    private double CalculateValidationLoss(Matrix X, VectorN y, int[] indices, int trainSize)
    {
        double totalLoss = 0;
        int valCount = indices.Length - trainSize;
        if (valCount <= 0) return 0;

        for (int i = trainSize; i < indices.Length; i++)
        {
            int idx = indices[i];
            var pred = _network.Forward(X.RowSlice(idx));
            totalLoss += Math.Pow(pred[0] - y[idx], 2);
        }
        return totalLoss / valCount;
    }

    public VectorN Predict(Matrix X)
    {
        var result = new double[X.rowLength];

        for (int i = 0; i < X.rowLength; i++)
        {
            var yhat = _network.Forward(X.RowSlice(i));
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

    public IModel Clone()
    {
        return new MLPRegressor
        {
            HiddenLayers = HiddenLayers == null ? null : (int[])HiddenLayers.Clone(),
            LearningRate = LearningRate,
            Epochs = Epochs,
            BatchSize = BatchSize,
            L2 = L2,
            ValidationSplit = ValidationSplit,
            Patience = Patience,
            MinDelta = MinDelta,
            Activation = Activation
        };
    }

    private void Shuffle(int[] array, Random rnd)
    {
        for (int i = array.Length - 1; i > 0; i--)
        {
            int j = rnd.Next(i + 1);
            (array[i], array[j]) = (array[j], array[i]);
        }
    }
}
