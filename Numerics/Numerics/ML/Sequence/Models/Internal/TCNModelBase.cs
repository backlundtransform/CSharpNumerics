using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.ML.NeuralNetwork;
using CSharpNumerics.ML.NeuralNetwork.Layers;
using CSharpNumerics.ML.NeuralNetwork.Layers.Interfaces;
using CSharpNumerics.ML.Sequence.Interfaces;
using CSharpNumerics.ML.Sequence.Layers;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.SingleObjective;
using CSharpNumerics.Numerics.Optimization.Strategies;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.Sequence.Models.Internal;

/// <summary>
/// Shared training/inference scaffolding for Temporal Convolutional Network models.
/// Builds a network of <c>TCNBlock → GlobalAvgPool → Dense(hidden) → Dense(output)</c> and
/// trains it with mini-batch SGD and early stopping, mirroring <see cref="CNN1DModelBase"/>.
/// </summary>
public abstract class TCNModelBase : ISequenceModel, IHasHyperparameters
{
    private SequentialModel _network;

    public int TimeSteps { get; set; }
    public int Features { get; set; } = 1;
    public int Channels { get; set; } = 16;
    public int KernelSize { get; set; } = 3;
    public int Levels { get; set; } = 4;
    public double DropoutRate { get; set; }
    public int HiddenUnits { get; set; } = 16;
    public ActivationType Activation { get; set; } = ActivationType.ReLU;
    public double LearningRate { get; set; } = 0.01;
    public int Epochs { get; set; } = 200;
    public int BatchSize { get; set; } = 16;
    public double L2 { get; set; }
    public double ValidationSplit { get; set; } = 0.2;
    public int Patience { get; set; } = 15;
    public double MinDelta { get; set; } = 1e-4;
    public int Seed { get; set; } = 123;

    public void Fit(Matrix X, VectorN y)
    {
        if (X.rowLength != y.Length)
            throw new ArgumentException("X.Rows must match y.Length.");

        ResolveSequenceDimensions(X.columnLength);
        _network = BuildNetwork(ResolveOutputSize(y));

        int[] indices = Enumerable.Range(0, X.rowLength).ToArray();
        var random = new Random(42);

        int trainSize = (ValidationSplit > 0 && X.rowLength > 5)
            ? (int)(X.rowLength * (1 - ValidationSplit))
            : X.rowLength;

        int valSize = X.rowLength - trainSize;
        var earlyStopping = new EarlyStopping(Patience, MinDelta);

        for (int epoch = 0; epoch < Epochs; epoch++)
        {
            Shuffle(indices, random);

            for (int batchStart = 0; batchStart < trainSize; batchStart += BatchSize)
            {
                int currentBatchSize = Math.Min(BatchSize, trainSize - batchStart);

                for (int offset = 0; offset < currentBatchSize; offset++)
                {
                    int sampleIndex = indices[batchStart + offset];
                    var output = _network.ForwardSingle(ToSequence(X.RowSlice(sampleIndex)));
                    var lossGradient = ComputeLossGradient(output, y[sampleIndex]);
                    _network.Backward(new[] { lossGradient });
                }

                var weightOptimizer = new GradientDescent(LearningRate, l2: L2);
                var biasOptimizer = new GradientDescent(LearningRate);
                _network.ApplyGradients(weightOptimizer, biasOptimizer, currentBatchSize);
            }

            if (valSize > 0)
            {
                double validationLoss = CalculateValidationLoss(X, y, indices, trainSize);
                if (earlyStopping.Check(validationLoss))
                    return;
            }
        }
    }

    public VectorN Predict(Matrix X)
    {
        if (_network == null)
            throw new InvalidOperationException("Model has not been fitted.");

        ResolveSequenceDimensions(X.columnLength);
        var predictions = new double[X.rowLength];

        for (int row = 0; row < X.rowLength; row++)
        {
            var output = _network.ForwardSingle(ToSequence(X.RowSlice(row)), training: false);
            predictions[row] = MapPrediction(output);
        }

        return new VectorN(predictions);
    }

    public virtual void SetHyperParameters(Dictionary<string, object> parameters)
    {
        if (parameters == null)
            return;

        if (parameters.TryGetValue("TimeSteps", out var timeSteps)) TimeSteps = Convert.ToInt32(timeSteps);
        if (parameters.TryGetValue("Features", out var features)) Features = Convert.ToInt32(features);
        if (parameters.TryGetValue("Channels", out var channels)) Channels = Convert.ToInt32(channels);
        if (parameters.TryGetValue("KernelSize", out var kernelSize)) KernelSize = Convert.ToInt32(kernelSize);
        if (parameters.TryGetValue("Levels", out var levels)) Levels = Convert.ToInt32(levels);
        if (parameters.TryGetValue("DropoutRate", out var dropoutRate)) DropoutRate = Convert.ToDouble(dropoutRate);
        if (parameters.TryGetValue("HiddenUnits", out var hiddenUnits)) HiddenUnits = Convert.ToInt32(hiddenUnits);
        if (parameters.TryGetValue("Activation", out var activation)) Activation = (ActivationType)activation;
        if (parameters.TryGetValue("LearningRate", out var learningRate)) LearningRate = Convert.ToDouble(learningRate);
        if (parameters.TryGetValue("Epochs", out var epochs)) Epochs = Convert.ToInt32(epochs);
        if (parameters.TryGetValue("BatchSize", out var batchSize)) BatchSize = Convert.ToInt32(batchSize);
        if (parameters.TryGetValue("L2", out var l2)) L2 = Convert.ToDouble(l2);
        if (parameters.TryGetValue("ValidationSplit", out var validationSplit)) ValidationSplit = Convert.ToDouble(validationSplit);
        if (parameters.TryGetValue("Patience", out var patience)) Patience = Convert.ToInt32(patience);
        if (parameters.TryGetValue("MinDelta", out var minDelta)) MinDelta = Convert.ToDouble(minDelta);
        if (parameters.TryGetValue("Seed", out var seed)) Seed = Convert.ToInt32(seed);
    }

    protected abstract int ResolveOutputSize(VectorN y);

    protected abstract VectorN ComputeLossGradient(VectorN output, double target);

    protected abstract double ComputeValidationLoss(VectorN output, double target);

    protected abstract double MapPrediction(VectorN output);

    public abstract IModel Clone();

    protected void CopySharedParametersTo(TCNModelBase clone)
    {
        clone.TimeSteps = TimeSteps;
        clone.Features = Features;
        clone.Channels = Channels;
        clone.KernelSize = KernelSize;
        clone.Levels = Levels;
        clone.DropoutRate = DropoutRate;
        clone.HiddenUnits = HiddenUnits;
        clone.Activation = Activation;
        clone.LearningRate = LearningRate;
        clone.Epochs = Epochs;
        clone.BatchSize = BatchSize;
        clone.L2 = L2;
        clone.ValidationSplit = ValidationSplit;
        clone.Patience = Patience;
        clone.MinDelta = MinDelta;
        clone.Seed = Seed;
    }

    private SequentialModel BuildNetwork(int outputSize)
    {
        var layers = new List<ILayer>
        {
            new TCNBlock(Features, Channels, KernelSize, Levels, DropoutRate, Activation, Seed),
            new GlobalAvgPool1DLayer()
        };

        int seed = Seed + 1000;
        int currentDenseInput = Channels;   // global average pool emits one vector of Channels

        if (HiddenUnits > 0)
        {
            layers.Add(new DenseLayer(currentDenseInput, HiddenUnits, Activation, seed++));
            currentDenseInput = HiddenUnits;
        }

        layers.Add(new DenseLayer(currentDenseInput, outputSize, ActivationType.Linear, seed));
        return new SequentialModel(layers.ToArray());
    }

    private VectorN[] ToSequence(VectorN flattened)
    {
        var sequence = new VectorN[TimeSteps];
        int index = 0;

        for (int timestep = 0; timestep < TimeSteps; timestep++)
        {
            var values = new double[Features];
            for (int feature = 0; feature < Features; feature++)
                values[feature] = flattened[index++];
            sequence[timestep] = new VectorN(values);
        }

        return sequence;
    }

    private void ResolveSequenceDimensions(int flattenedWidth)
    {
        if (TimeSteps <= 0 && Features <= 0)
        {
            Features = 1;
            TimeSteps = flattenedWidth;
        }
        else if (TimeSteps <= 0)
        {
            if (flattenedWidth % Features != 0)
                throw new ArgumentException("Flattened feature width must be divisible by the configured feature count.");
            TimeSteps = flattenedWidth / Features;
        }
        else if (Features <= 0)
        {
            if (flattenedWidth % TimeSteps != 0)
                throw new ArgumentException("Flattened feature width must be divisible by the configured timestep count.");
            Features = flattenedWidth / TimeSteps;
        }

        if ((TimeSteps * Features) != flattenedWidth)
            throw new ArgumentException("Flattened sequence width must equal TimeSteps * Features.");
    }

    private double CalculateValidationLoss(Matrix X, VectorN y, int[] indices, int trainSize)
    {
        int validationCount = indices.Length - trainSize;
        if (validationCount <= 0)
            return 0.0;

        double totalLoss = 0.0;
        for (int i = trainSize; i < indices.Length; i++)
        {
            int sampleIndex = indices[i];
            var output = _network.ForwardSingle(ToSequence(X.RowSlice(sampleIndex)), training: false);
            totalLoss += ComputeValidationLoss(output, y[sampleIndex]);
        }

        return totalLoss / validationCount;
    }

    private static void Shuffle(int[] values, Random random)
    {
        for (int i = values.Length - 1; i > 0; i--)
        {
            int j = random.Next(i + 1);
            (values[i], values[j]) = (values[j], values[i]);
        }
    }
}
