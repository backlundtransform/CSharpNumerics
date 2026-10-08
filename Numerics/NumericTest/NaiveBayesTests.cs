using CSharpNumerics.ML.Models.Classification;
using CSharpNumerics.Numerics.Objects;

namespace NumericTest;

/// <summary>
/// Tests for the Gaussian naive Bayes classifier, which had no coverage before.
/// </summary>
/// <remarks>
/// Added alongside the fix to <c>NumClasses</c>, which threw
/// <c>NotImplementedException</c> — the only such member left in the library.
/// </remarks>
[TestClass]
public class NaiveBayesTests
{
    /// <summary>
    /// Two well-separated Gaussian blobs: class 0 near the origin, class 1 far from it.
    /// </summary>
    private static (Matrix x, VectorN y) TwoBlobs()
    {
        var x = new double[,]
        {
            { 0.0, 0.0 }, { 0.2, -0.1 }, { -0.1, 0.2 }, { 0.1, 0.1 },
            { 5.0, 5.0 }, { 5.2, 4.9 }, { 4.9, 5.2 }, { 5.1, 5.1 }
        };

        var y = new double[] { 0, 0, 0, 0, 1, 1, 1, 1 };

        return (new Matrix(x), new VectorN(y));
    }

    [TestMethod]
    public void NumClasses_IsSetByFit()
    {
        var (x, y) = TwoBlobs();

        var model = new NaiveBayes();
        model.Fit(x, y);

        Assert.AreEqual(2, model.NumClasses);
    }

    [TestMethod]
    public void NumClasses_CountsTheHighestLabelPlusOne()
    {
        // Labels 0 and 2 with nothing at 1: NumClasses sizes a label-indexed array,
        // so it must be 3, matching the other classifiers.
        var x = new Matrix(new double[,] { { 0.0 }, { 1.0 }, { 8.0 }, { 9.0 } });
        var y = new VectorN(new double[] { 0, 0, 2, 2 });

        var model = new NaiveBayes();
        model.Fit(x, y);

        Assert.AreEqual(3, model.NumClasses);
    }

    [TestMethod]
    public void Predict_SeparatesTwoWellSeparatedBlobs()
    {
        var (x, y) = TwoBlobs();

        var model = new NaiveBayes();
        model.Fit(x, y);

        var predictions = model.Predict(x);

        for (var i = 0; i < y.Length; i++)
        {
            Assert.AreEqual(y[i], predictions[i], $"Sample {i}");
        }
    }

    [TestMethod]
    public void Predict_AssignsUnseenPointsToTheNearerClass()
    {
        var (x, y) = TwoBlobs();

        var model = new NaiveBayes();
        model.Fit(x, y);

        var unseen = new Matrix(new double[,] { { 0.3, 0.3 }, { 4.7, 5.3 } });
        var predictions = model.Predict(unseen);

        Assert.AreEqual(0.0, predictions[0]);
        Assert.AreEqual(1.0, predictions[1]);
    }

    [TestMethod]
    public void PredictBeforeFit_Throws()
    {
        var model = new NaiveBayes();

        Assert.ThrowsException<InvalidOperationException>(
            () => model.Predict(new Matrix(new double[,] { { 1.0, 2.0 } })));
    }

    [TestMethod]
    public void Clone_ReturnsAnUnfittedModel()
    {
        var (x, y) = TwoBlobs();

        var model = new NaiveBayes();
        model.Fit(x, y);

        // Clone follows the convention of the other models: hyperparameters only, no fitted
        // state, so grid search starts each fold from a clean estimator.
        var clone = (NaiveBayes)model.Clone();

        Assert.ThrowsException<InvalidOperationException>(() => clone.Predict(x));
    }
}
