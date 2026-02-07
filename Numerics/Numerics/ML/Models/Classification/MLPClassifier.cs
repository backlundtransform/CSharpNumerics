using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Objects;
using Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.Models.Classification;

    public class MLPClassifier : IClassificationModel, IHasHyperparameters
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
  

    private List<Matrix> _weights;
    private List<VectorN> _biases;
        

        private void Initialize(int inputSize, int outputSize)
        {
            _weights = new();
            _biases = new();
       

            var layers = new List<int> { inputSize };
            layers.AddRange(HiddenLayers);
            layers.Add(outputSize); // K klasser

            var rnd = new Random(123);
            for (int i = 0; i < layers.Count - 1; i++)
            {
                _weights.Add(RandomMatrix(layers[i], layers[i + 1], rnd));
                _biases.Add(new VectorN(layers[i + 1]));
            }
        }

        private VectorN Softmax(VectorN z)
        {
            var values = z.Values;
            double max = values.Max();
            double sum = 0;
            double[] expValues = new double[values.Length];
            for (int i = 0; i < values.Length; i++)
            {
                expValues[i] = Math.Exp(values[i] - max);
                sum += expValues[i];
            }
            for (int i = 0; i < values.Length; i++) expValues[i] /= sum;
            return new VectorN(expValues);
        }

        private VectorN Forward(VectorN x, out List<VectorN> activations)
        {
            activations = new() { x };
            for (int i = 0; i < _weights.Count; i++)
            {
                var z = (_weights[i].Transpose() * activations[^1]) + _biases[i];
               
                var a = (i == _weights.Count - 1) ? Softmax(z) : Activate(z);
                activations.Add(a);
            }
            return activations[^1];
        }

        public void Fit(Matrix X, VectorN y)
        {
            var numClasses = (int)y.Values.Max() + 1;
            NumClasses = numClasses;
            Initialize(X.columnLength, numClasses);

            int[] indices = Enumerable.Range(0, X.rowLength).ToArray();
            var rnd = new Random(42);

            int trainSize = (ValidationSplit > 0 && X.rowLength > 5)
                ? (int)(X.rowLength * (1 - ValidationSplit))
                : X.rowLength;
            int valSize = X.rowLength - trainSize;

            List<Matrix> bestWeights = _weights.Select(m => new Matrix(m.values)).ToList();
            List<VectorN> bestBiases = _biases.Select(v => new VectorN(v.Values)).ToList();
            double bestValLoss = double.MaxValue;
            int patienceCounter = 0;

            for (int epoch = 0; epoch < Epochs; epoch++)
            {
                Shuffle(indices, rnd);

                for (int i = 0; i < trainSize; i += BatchSize)
                {
                    int currentBatchSize = Math.Min(BatchSize, trainSize - i);
                    var weightGrads = InitGradsLike(_weights);
                    var biasGrads = InitGradsLike(_biases);

                    for (int b = 0; b < currentBatchSize; b++)
                    {
                        int idx = indices[i + b];
                        var xi = X.RowSlice(idx);

                       
                        var yi_oh = new double[numClasses];
                        yi_oh[(int)y[idx]] = 1.0;
                        var target = new VectorN(yi_oh);

                        var yhat = Forward(xi, out var acts);
                        ComputeGradients(target, acts, out var dw, out var db);

                        for (int l = 0; l < _weights.Count; l++)
                        {
                            weightGrads[l] += dw[l];
                            biasGrads[l] += db[l];
                        }
                    }

                    for (int l = 0; l < _weights.Count; l++)
                    {
                        _weights[l] -= (LearningRate / currentBatchSize) * (weightGrads[l] + L2 * _weights[l]);
                        _biases[l] -= (LearningRate / currentBatchSize) * biasGrads[l];
                    }
                }

                if (valSize > 0)
                {
                    double currentValLoss = CalculateValidationLoss(X, y, indices, trainSize, numClasses);
                    if (currentValLoss < bestValLoss - MinDelta)
                    {
                        bestValLoss = currentValLoss;
                        patienceCounter = 0;
                        bestWeights = _weights.Select(m => new Matrix(m.values)).ToList();
                        bestBiases = _biases.Select(v => new VectorN(v.Values)).ToList();
                    }
                    else { patienceCounter++; }

                    if (patienceCounter >= Patience)
                    {
                        _weights = bestWeights; _biases = bestBiases;
                        return;
                    }
                }
            }
            if (valSize > 0) { _weights = bestWeights; _biases = bestBiases; }
        }

        public VectorN Predict(Matrix X)
        {
            double[] result = new double[X.rowLength];
            for (int i = 0; i < X.rowLength; i++)
            {
                var probs = Forward(X.RowSlice(i), out _);
                result[i] = Array.IndexOf(probs.Values, probs.Values.Max()); 
            }
            return new VectorN(result);
        }

        private double CalculateValidationLoss(Matrix X, VectorN y, int[] indices, int trainSize, int numClasses)
        {
            double totalLoss = 0;
            int valCount = indices.Length - trainSize;
            for (int i = trainSize; i < indices.Length; i++)
            {
                int idx = indices[i];
                var probs = Forward(X.RowSlice(idx), out _);
              
                totalLoss -= Math.Log(probs.Values[(int)y[idx]] + 1e-15);
            }
            return totalLoss / valCount;
        }

        private void ComputeGradients(VectorN target, List<VectorN> activations, out List<Matrix> dW, out List<VectorN> dB)
        {
            var deltas = new List<VectorN>();
            var error = activations[^1] - target; 
            deltas.Add(error);

            for (int i = _weights.Count - 2; i >= 0; i--)
            {
                var w = _weights[i + 1];
                var delta = (w * deltas[^1]).Hadamard(ActivationDerivative(activations[i + 1]));
                deltas.Add(delta);
            }
            deltas.Reverse();

            dW = new List<Matrix>();
            dB = new List<VectorN>();
            for (int i = 0; i < _weights.Count; i++)
            {
                dW.Add(activations[i].Outer(deltas[i]));
                dB.Add(deltas[i]);
            }
        }

        private List<Matrix> InitGradsLike(List<Matrix> template) =>
        template.Select(m => new Matrix { values = new double[m.rowLength, m.columnLength], rowLength = m.rowLength, columnLength = m.columnLength }).ToList();

        private List<VectorN> InitGradsLike(List<VectorN> template) =>
            template.Select(v => new VectorN(v.Length)).ToList();
     

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
       
        double limit = Math.Sqrt(6.0 / (rows + cols));
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                values[i, j] = (rnd.NextDouble() * 2 - 1) * limit;

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

                ActivationType.Sigmoid => x * (1.0 - x),

                ActivationType.Tanh => 1.0 - x * x,

                ActivationType.Linear => 1.0,
                _ => 1.0
            };
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

        private void Shuffle(int[] array, Random rnd)
        {
            for (int i = array.Length - 1; i > 0; i--)
            {
                int j = rnd.Next(i + 1);
                (array[i], array[j]) = (array[j], array[i]);
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
    
