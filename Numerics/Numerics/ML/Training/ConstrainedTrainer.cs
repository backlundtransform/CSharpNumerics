using CSharpNumerics.ML.NeuralNetwork;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.SingleObjective;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.Training;

/// <summary>
/// Mini-batch SGD trainer that fits a model under physics constraints. The total objective is a
/// data-fidelity term (squared error by default) plus a constraint loss (typically a
/// <see cref="CompositeLoss"/> of non-negativity / conservation / smoothness). It supports
/// <em>curriculum learning</em>: the constraint contribution is ramped from 0 up to its full
/// weight over <see cref="CurriculumWarmupEpochs"/>, so the model first learns to fit the data
/// and only gradually has the constraints tightened.
/// <para>
/// The model maps each input row to a component vector: a row is presented as a single-timestep
/// sequence and the network must emit a single output vector.
/// </para>
/// </summary>
public class ConstrainedTrainer
{
    private readonly ILoss _constraintLoss;
    private readonly ILoss _dataLoss;

    /// <param name="constraintLoss">The constraint term (e.g. a <see cref="CompositeLoss"/>).</param>
    /// <param name="dataLoss">Optional data-fidelity loss; squared error is used when null.</param>
    public ConstrainedTrainer(ILoss constraintLoss, ILoss dataLoss = null)
    {
        _constraintLoss = constraintLoss ?? throw new ArgumentNullException(nameof(constraintLoss));
        _dataLoss = dataLoss;   // null → built-in squared error
    }

    public int Epochs { get; set; } = 200;
    public double LearningRate { get; set; } = 0.01;
    public int BatchSize { get; set; } = 16;
    public double L2 { get; set; }

    /// <summary>
    /// Epochs over which the constraint weight ramps linearly from 0 to 1. Zero disables the
    /// curriculum, applying full constraints from the first epoch.
    /// </summary>
    public int CurriculumWarmupEpochs { get; set; }

    public int Seed { get; set; } = 42;

    /// <summary>
    /// Trains <paramref name="model"/> to map rows of <paramref name="X"/> to the corresponding
    /// rows of <paramref name="targetComponents"/>. Returns the average total loss per epoch.
    /// </summary>
    public IReadOnlyList<double> Train(SequentialModel model, Matrix X, Matrix targetComponents)
    {
        if (model == null) throw new ArgumentNullException(nameof(model));
        if (X.rowLength != targetComponents.rowLength)
            throw new ArgumentException("X and targetComponents must have the same number of rows.");

        int n = X.rowLength;
        int[] indices = Enumerable.Range(0, n).ToArray();
        var random = new Random(Seed);
        var history = new List<double>(Epochs);

        for (int epoch = 0; epoch < Epochs; epoch++)
        {
            double curriculum = CurriculumWarmupEpochs <= 0
                ? 1.0
                : Math.Min(1.0, (epoch + 1.0) / CurriculumWarmupEpochs);

            Shuffle(indices, random);
            double epochLoss = 0.0;

            for (int batchStart = 0; batchStart < n; batchStart += BatchSize)
            {
                int currentBatchSize = Math.Min(BatchSize, n - batchStart);

                for (int offset = 0; offset < currentBatchSize; offset++)
                {
                    int i = indices[batchStart + offset];
                    var input = new[] { new VectorN(X.RowSlice(i).Values) };
                    var target = new VectorN(targetComponents.RowSlice(i).Values);

                    VectorN prediction = model.ForwardSingle(input);

                    VectorN dataGrad = DataGradient(prediction, target);
                    VectorN constraintGrad = _constraintLoss.Gradient(prediction, target);

                    var grad = new double[prediction.Length];
                    for (int k = 0; k < grad.Length; k++)
                        grad[k] = dataGrad[k] + (curriculum * constraintGrad[k]);

                    model.Backward(new[] { new VectorN(grad) });

                    epochLoss += DataValue(prediction, target)
                        + (curriculum * _constraintLoss.Compute(prediction, target));
                }

                var weightOptimizer = new GradientDescent(LearningRate, l2: L2);
                var biasOptimizer = new GradientDescent(LearningRate);
                model.ApplyGradients(weightOptimizer, biasOptimizer, currentBatchSize);
            }

            history.Add(epochLoss / n);
        }

        return history;
    }

    private VectorN DataGradient(VectorN prediction, VectorN target)
    {
        if (_dataLoss != null)
            return _dataLoss.Gradient(prediction, target);

        // Squared error: ∂/∂pred ½‖pred − target‖² = pred − target.
        var grad = new double[prediction.Length];
        for (int i = 0; i < grad.Length; i++)
            grad[i] = prediction[i] - target[i];
        return new VectorN(grad);
    }

    private double DataValue(VectorN prediction, VectorN target)
    {
        if (_dataLoss != null)
            return _dataLoss.Compute(prediction, target);

        double sum = 0.0;
        for (int i = 0; i < prediction.Length; i++)
        {
            double d = prediction[i] - target[i];
            sum += d * d;
        }
        return 0.5 * sum;
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
