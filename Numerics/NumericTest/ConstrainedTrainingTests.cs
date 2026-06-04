using System;
using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.NeuralNetwork;
using CSharpNumerics.ML.NeuralNetwork.Layers;
using CSharpNumerics.ML.Training;
using CSharpNumerics.Numerics.Objects;

namespace NumericsTests;

[TestClass]
public class ConstrainedTrainingTests
{
    #region NonNegativityLoss

    [TestMethod]
    public void NonNegativityLoss_IsZero_WhenAllPredictionsNonNegative()
    {
        var loss = new NonNegativityLoss(weight: 1.0);
        var prediction = new VectorN(new[] { 0.0, 1.5, 3.2, 10.0 });
        var target = new VectorN(new[] { 0.0, 0.0, 0.0, 0.0 });

        Assert.AreEqual(0.0, loss.Compute(prediction, target), 1e-12);

        // A negative entry contributes its square; gradient is 2·w·y there, 0 elsewhere.
        var withNegative = new VectorN(new[] { 1.0, -2.0, 3.0 });
        Assert.AreEqual(4.0, loss.Compute(withNegative, target), 1e-12);

        var grad = loss.Gradient(withNegative, target);
        Assert.AreEqual(0.0, grad[0], 1e-12);
        Assert.AreEqual(-4.0, grad[1], 1e-12);
        Assert.AreEqual(0.0, grad[2], 1e-12);
    }

    #endregion

    #region ConservationLoss

    [TestMethod]
    public void ConservationLoss_IsZero_WhenComponentsSumToTotal()
    {
        var loss = new ConservationLoss(weight: 1.0);
        var components = new VectorN(new[] { 3.0, 3.0, 4.0 });  // Σ = 10
        var target = new VectorN(new[] { 5.0, 5.0 });            // Σ = 10

        Assert.AreEqual(0.0, loss.Compute(components, target), 1e-12);

        // Mismatched sums: penalty = (Σtarget − Σpred)², gradient = 2(Σpred − Σtarget) per entry.
        var over = new VectorN(new[] { 4.0, 4.0, 4.0 });  // Σ = 12, target Σ = 10 → diff = 2
        Assert.AreEqual(4.0, loss.Compute(over, target), 1e-12);

        var grad = loss.Gradient(over, target);
        foreach (var g in grad.Values)
            Assert.AreEqual(4.0, g, 1e-12);   // 2 · 1 · (12 − 10)
    }

    #endregion

    #region CompositeLoss

    [TestMethod]
    public void CompositeLoss_Gradient_EqualsWeightedSumOfSubGradients()
    {
        var nonNeg = new NonNegativityLoss(weight: 1.0);
        var conservation = new ConservationLoss(weight: 1.0);

        const double w1 = 0.7, w2 = 2.5;
        var composite = new CompositeLoss()
            .Add(nonNeg, w1)
            .Add(conservation, w2);

        var prediction = new VectorN(new[] { -1.0, 2.0, 5.0 });
        var target = new VectorN(new[] { 2.0, 2.0, 2.0 });   // Σ = 6

        var compositeGrad = composite.Gradient(prediction, target);
        var nonNegGrad = nonNeg.Gradient(prediction, target);
        var consGrad = conservation.Gradient(prediction, target);

        for (int i = 0; i < prediction.Length; i++)
        {
            double expected = (w1 * nonNegGrad[i]) + (w2 * consGrad[i]);
            Assert.AreEqual(expected, compositeGrad[i], 1e-12, $"Component {i}");
        }

        // The composite value is likewise the weighted sum.
        double expectedValue = (w1 * nonNeg.Compute(prediction, target))
            + (w2 * conservation.Compute(prediction, target));
        Assert.AreEqual(expectedValue, composite.Compute(prediction, target), 1e-12);
    }

    #endregion

    #region ConstrainedTrainer

    [TestMethod]
    public void ConstrainedTrainer_WithConservation_ProducesSumWithinOnePercentOfTotal()
    {
        const int k = 3;
        const double total = 10.0;
        int samples = 40;
        var rng = new Random(17);

        var xValues = new double[samples, k];
        var targetValues = new double[samples, k];
        for (int i = 0; i < samples; i++)
        {
            double a = rng.NextDouble() + 0.1;
            double b = rng.NextDouble() + 0.1;
            double c = rng.NextDouble() + 0.1;
            double s = a + b + c;
            // Scale so each row's components sum exactly to the total.
            double[] comps = { total * a / s, total * b / s, total * c / s };
            for (int j = 0; j < k; j++)
            {
                xValues[i, j] = comps[j];
                targetValues[i, j] = comps[j];
            }
        }

        var model = new SequentialModel(
            new DenseLayer(k, 8, ActivationType.ReLU, seed: 1),
            new DenseLayer(8, k, ActivationType.Linear, seed: 2));

        var constraints = new CompositeLoss().Add(new ConservationLoss(weight: 10.0));
        var trainer = new ConstrainedTrainer(constraints)
        {
            Epochs = 400,
            LearningRate = 0.01,
            BatchSize = 16
        };

        trainer.Train(model, new Matrix(xValues), new Matrix(targetValues));

        var predictions = PredictRows(model, new Matrix(xValues), k);
        for (int i = 0; i < samples; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < k; j++) sum += predictions[i, j];
            Assert.IsTrue(Math.Abs(sum - total) < 0.01 * total,
                $"Sample {i}: component sum {sum:F4} should be within 1% of {total}.");
        }
    }

    #endregion

    #region SoftmaxConstraintHead

    [TestMethod]
    public void SoftmaxConstraintHead_ProducesValidPartition()
    {
        var head = new SoftmaxConstraintHead();
        var input = new[]
        {
            new VectorN(new[] { 1.0, 2.0, 3.0 }),
            new VectorN(new[] { 0.0, 0.0, 0.0 }),
            new VectorN(new[] { -5.0, 10.0, 2.0 }),
        };

        var output = head.Forward(input, training: true);

        foreach (var v in output)
        {
            double sum = 0.0;
            foreach (var x in v.Values)
            {
                Assert.IsTrue(x >= 0.0 && x <= 1.0, $"Output {x} must lie in [0, 1].");
                sum += x;
            }
            Assert.IsTrue(sum <= 1.0 + 1e-9, $"Output sum {sum} must not exceed 1.");
            Assert.AreEqual(1.0, sum, 1e-9, "Softmax outputs should sum to 1.");
        }

        // Backward pass is finite and shape-preserving.
        var grad = head.Backward(output);
        Assert.AreEqual(input.Length, grad.Length);
        foreach (var v in grad)
            foreach (var x in v.Values)
                Assert.IsFalse(double.IsNaN(x) || double.IsInfinity(x));
    }

    #endregion

    private static double[,] PredictRows(SequentialModel model, Matrix X, int outputSize)
    {
        var result = new double[X.rowLength, outputSize];
        for (int i = 0; i < X.rowLength; i++)
        {
            var output = model.ForwardSingle(new[] { new VectorN(X.RowSlice(i).Values) }, training: false);
            for (int j = 0; j < outputSize; j++)
                result[i, j] = output[j];
        }
        return result;
    }
}
