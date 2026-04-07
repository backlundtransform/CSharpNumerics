using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.ML.NeuralNetwork;
using CSharpNumerics.ML.NeuralNetwork.Layers;
using CSharpNumerics.ML.Sequence.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.SingleObjective;

namespace NumericTest;

[TestClass]
public class SequenceModelInfrastructureTests
{
    [TestMethod]
    public void VectorN_Concat_ShouldPreserveValueOrder()
    {
        var left = new VectorN(new[] { 1.0, 2.0 });
        var right = new VectorN(new[] { 3.0, 4.0, 5.0 });

        var concatenated = left.Concat(right);

        CollectionAssert.AreEqual(new[] { 1.0, 2.0, 3.0, 4.0, 5.0 }, concatenated.Values);
    }

    [TestMethod]
    public void Activations_ShouldExposeReusableFunctions()
    {
        var input = new VectorN(new[] { -1.0, 0.0, 1.0 });

        var relu = Activations.Apply(input, ActivationType.ReLU);
        var sigmoid = Activations.Apply(input, ActivationType.Sigmoid);
        var softmax = Activations.Softmax(input);
        var tanhDerivative = Activations.Derivative(Activations.Apply(input, ActivationType.Tanh), ActivationType.Tanh);

        CollectionAssert.AreEqual(new[] { 0.0, 0.0, 1.0 }, relu.Values);
        Assert.AreEqual(0.2689414213699951, sigmoid[0], 1e-12);
        Assert.AreEqual(0.7310585786300049, sigmoid[2], 1e-12);
        Assert.AreEqual(1.0, softmax.Values.Sum(), 1e-12);
        Assert.AreEqual(0.41997434161402614, tanhDerivative[0], 1e-12);
    }

    [TestMethod]
    public void DenseLayer_ApplyGradients_ShouldReduceTrainingLoss()
    {
        var weights = new Matrix(new double[,]
        {
            { 1.0, 0.0 },
            { 0.0, 1.0 }
        });
        var biases = new VectorN(new[] { 0.0, 0.0 });
        var layer = new DenseLayer(weights, biases, ActivationType.Linear);
        var input = new[] { new VectorN(new[] { 1.0, 2.0 }) };
        var target = new VectorN(new[] { 0.0, 0.0 });

        var before = layer.Forward(input)[0];
        var beforeLoss = SquaredError(before, target);

        var lossGradient = before - target;
        layer.Backward(new[] { lossGradient });
        layer.ApplyGradients(new GradientDescent(learningRate: 0.1), new GradientDescent(learningRate: 0.1), batchSize: 1);

        var after = layer.Forward(input)[0];
        var afterLoss = SquaredError(after, target);

        Assert.IsTrue(afterLoss < beforeLoss, $"Expected the updated layer to reduce loss, but {afterLoss} >= {beforeLoss}.");
    }

    [TestMethod]
    public void SequentialModel_ShouldComposeForwardBackwardAndParameterCounts()
    {
        var firstLayer = new DenseLayer(
            new Matrix(new double[,]
            {
                { 1.0, 0.0 },
                { 0.0, 1.0 }
            }),
            new VectorN(new[] { 0.0, 0.0 }),
            ActivationType.Linear);

        var secondLayer = new DenseLayer(
            new Matrix(new double[,]
            {
                { 2.0 },
                { 3.0 }
            }),
            new VectorN(new[] { 0.0 }),
            ActivationType.Linear);

        var model = new SequentialModel(firstLayer, secondLayer);
        var input = new[] { new VectorN(new[] { 4.0, 5.0 }) };

        var output = model.ForwardSingle(input);
        var inputGradient = model.BackwardSingle(new VectorN(new[] { 1.0 }));

        Assert.AreEqual(23.0, output[0], 1e-12);
        CollectionAssert.AreEqual(new[] { 2.0, 3.0 }, inputGradient.Values);
        Assert.AreEqual(9, model.ParameterCount);
    }

    [TestMethod]
    public void SequenceModelInterface_ShouldRemainCompatibleWithIModel()
    {
        ISequenceModel model = new StubSequenceModel
        {
            TimeSteps = 12,
            Features = 3
        };

        Assert.AreEqual(12, model.TimeSteps);
        Assert.AreEqual(3, model.Features);
        Assert.IsTrue(model is IModel);
    }

    private static double SquaredError(VectorN prediction, VectorN target)
    {
        double sum = 0.0;
        for (int i = 0; i < prediction.Length; i++)
        {
            double error = prediction[i] - target[i];
            sum += error * error;
        }

        return sum;
    }

    private sealed class StubSequenceModel : ISequenceModel
    {
        public int TimeSteps { get; set; }

        public int Features { get; set; }

        public void Fit(Matrix X, VectorN y)
        {
        }

        public VectorN Predict(Matrix X)
        {
            return new VectorN(X.rowLength);
        }

        public IModel Clone()
        {
            return new StubSequenceModel
            {
                TimeSteps = TimeSteps,
                Features = Features
            };
        }
    }
}