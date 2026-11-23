using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Strategies;
using System;
using System.Collections.Generic;
using System.Linq;
using NN = CSharpNumerics.ML.NeuralNetwork.NeuralNetwork;

namespace CSharpNumerics.ML.Models.Classification;

    public class MLPClassifier : IProbabilisticClassifier, ISampleWeightedModel, IHasHyperparameters
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

        public int NumClasses { get; private set; }

        private NN _network;

        public void Fit(Matrix X, VectorN y) => FitWeighted(X, y, null);

        /// <summary>
        /// Trains with a weight per row: each row's loss (and gradient) is scaled by its weight.
        /// The weights are rescaled to average 1, so the learning rate keeps its meaning.
        /// </summary>
        public void Fit(Matrix X, VectorN y, VectorN sampleWeights) => FitWeighted(X, y, sampleWeights);

        private void FitWeighted(Matrix X, VectorN y, VectorN? sampleWeights)
        {
            double[] weights = ClassWeights.Validate(sampleWeights, y.Length);
            double mean = weights.Average();
            for (int i = 0; i < weights.Length; i++) weights[i] /= mean;

            var numClasses = (int)y.Values.Max() + 1;
            NumClasses = numClasses;
            _network = new NN(X.columnLength, HiddenLayers, numClasses,
                Activation, NN.OutputMode.Softmax);

            int[] indices = Enumerable.Range(0, X.rowLength).ToArray();
            var rnd = new Random(42);

            int trainSize = (ValidationSplit > 0 && X.rowLength > 5)
                ? (int)(X.rowLength * (1 - ValidationSplit))
                : X.rowLength;
            int valSize = X.rowLength - trainSize;

            var (bestWeights, bestBiases) = _network.SnapshotWeights();
            var earlyStopping = new EarlyStopping(Patience, MinDelta);

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

                        var yi_oh = new double[numClasses];
                        yi_oh[(int)y[idx]] = 1.0;
                        var target = new VectorN(yi_oh);

                        var yhat = _network.Forward(xi, out var acts);
                        _network.ComputeGradients(target, acts, out var dw, out var db);
                        if (weights[idx] != 1.0)
                        {
                            for (int l = 0; l < dw.Count; l++)
                            {
                                dw[l] = weights[idx] * dw[l];
                                db[l] = weights[idx] * db[l];
                            }
                        }
                        _network.AccumulateGradients(weightGrads, biasGrads, dw, db);
                    }

                    _network.ApplyGradients(weightGrads, biasGrads, LearningRate, currentBatchSize, L2);
                }

                if (valSize > 0)
                {
                    double currentValLoss = CalculateValidationLoss(X, y, weights, indices, trainSize);

                    if (!earlyStopping.ShouldStop && currentValLoss < earlyStopping.BestValue - MinDelta)
                        (bestWeights, bestBiases) = _network.SnapshotWeights();

                    if (earlyStopping.Check(currentValLoss))
                    {
                        _network.RestoreWeights(bestWeights, bestBiases);
                        return;
                    }
                }
            }
            if (valSize > 0) { _network.RestoreWeights(bestWeights, bestBiases); }
        }

        /// <summary>Softmax output: one probability per class for each row.</summary>
        public Matrix PredictProbabilities(Matrix X)
        {
            if (_network == null)
                throw new InvalidOperationException("Model has not been fitted.");

            var result = new Matrix(X.rowLength, NumClasses);
            for (int i = 0; i < X.rowLength; i++)
            {
                var probs = _network.Forward(X.RowSlice(i));
                for (int c = 0; c < NumClasses; c++)
                    result.values[i, c] = probs[c];
            }
            return result;
        }

        public VectorN Predict(Matrix X)
        {
            double[] result = new double[X.rowLength];
            for (int i = 0; i < X.rowLength; i++)
            {
                var probs = _network.Forward(X.RowSlice(i));
                result[i] = Array.IndexOf(probs.Values, probs.Values.Max());
            }
            return new VectorN(result);
        }

        private double CalculateValidationLoss(Matrix X, VectorN y, double[] weights, int[] indices, int trainSize)
        {
            double totalLoss = 0, totalWeight = 0;
            for (int i = trainSize; i < indices.Length; i++)
            {
                int idx = indices[i];
                var probs = _network.Forward(X.RowSlice(idx));
                totalLoss -= weights[idx] * Math.Log(probs.Values[(int)y[idx]] + 1e-15);
                totalWeight += weights[idx];
            }
            return totalWeight > 0 ? totalLoss / totalWeight : 0.0;
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
            return new MLPClassifier
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
    
