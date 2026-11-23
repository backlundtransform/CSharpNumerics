using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.ML.NeuralNetwork;
using CSharpNumerics.ML.Sequence.Models.Internal;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.Sequence.Models.Classification;

public class BiLSTMClassifier : BiLSTMModelBase, IProbabilisticClassifier
{
    public int NumClasses { get; private set; } = 2;

    public override void SetHyperParameters(Dictionary<string, object> parameters)
    {
        base.SetHyperParameters(parameters);
        if (parameters != null && parameters.TryGetValue("NumClasses", out var numClasses))
            NumClasses = Convert.ToInt32(numClasses);
    }

    protected override int ResolveOutputSize(VectorN y)
    {
        NumClasses = Math.Max(NumClasses, (int)y.Values.Max() + 1);
        return NumClasses;
    }

    protected override VectorN ComputeLossGradient(VectorN output, double target)
    {
        var probabilities = Activations.Softmax(output);
        var targetVector = new double[NumClasses];
        targetVector[(int)target] = 1.0;
        return probabilities - new VectorN(targetVector);
    }

    protected override double ComputeValidationLoss(VectorN output, double target)
    {
        var probabilities = Activations.Softmax(output);
        return -Math.Log(probabilities[(int)target] + 1e-15);
    }

    /// <summary>Softmax of the network output: one probability per class for each row.</summary>
    public Matrix PredictProbabilities(Matrix X)
    {
        var outputs = ForwardRows(X);
        var result = new Matrix(outputs.Length, NumClasses);
        for (int row = 0; row < outputs.Length; row++)
        {
            var probabilities = Activations.Softmax(outputs[row]);
            for (int c = 0; c < NumClasses; c++)
                result.values[row, c] = probabilities[c];
        }
        return result;
    }

    protected override double MapPrediction(VectorN output)
    {
        var probabilities = Activations.Softmax(output);
        return Array.IndexOf(probabilities.Values, probabilities.Values.Max());
    }

    public override IModel Clone()
    {
        var clone = new BiLSTMClassifier { NumClasses = NumClasses };
        CopySharedParametersTo(clone);
        return clone;
    }
}
