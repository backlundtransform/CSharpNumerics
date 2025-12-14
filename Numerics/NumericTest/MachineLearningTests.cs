using CSharpNumerics.ML;
using CSharpNumerics.ML.Models.Classification;
using CSharpNumerics.ML.Models.Regression;
using CSharpNumerics.ML.Scalers;
using CSharpNumerics.Objects;
using Numerics.Objects;
using Xunit;
using Xunit.Sdk;

namespace NumericTest
{
    [TestClass]
    public class MachineLearningTests
    {
        [TestMethod]
        public void LinearRegression_RollingCV_PerfectLine_ShouldHaveNearZeroMSE()
        {
            // Arrange: y = 2x + 1
            double[,] Xdata =
            {
                { 0 },
                { 1 },
                { 2 },
                { 3 },
                { 4 },
                { 5 },
                { 6 },
                { 7 },
                { 8 },
                { 9 }
            };

            double[] ydata =
            {
                1, 3, 5, 7, 9, 11, 13, 15, 17, 19
            };

            var X = new Matrix(Xdata);
            var y = new VectorN(ydata);

            var model = new Linear(fitIntercept: true);

            var modelParams = new Dictionary<string, object>
            {
                { "LearningRate", 0.01 }
            };

            var pipeline = new Pipeline(
                model: model,
                modelParams: modelParams
            );

            var cv = new RollingCrossValidator(
                pipelines: new List<Pipeline> { pipeline },
                folds: 5
            );

            // Act
            var result = cv.Run(X, y);

            // Assert
            Assert.IsNotNull(result.BestPipeline);
            Assert.IsTrue(result.BestScore > -1e-10,
                $"Expected near-zero negative MSE, got {result.BestScore}");
        }


        [TestMethod]
        public void Pipeline_FitAndPredict_PerfectLinearData_ShouldPredictCorrectly()
        {
            // Arrange: y = 2x + 1
            double[,] Xdata =
            {
        { 0 },
        { 1 },
        { 2 },
        { 3 },
        { 4 }
    };

            double[] ydata = { 1, 3, 5, 7, 9 };

            var X = new Matrix(Xdata);
            var y = new VectorN(ydata);

            var model = new Linear(fitIntercept: true);

            var pipeline = new Pipeline(
                model: model,
                modelParams: new Dictionary<string, object>()
            );

            // Act
            pipeline.Fit(X, y);
            var preds = pipeline.Predict(X);

            // Assert
            Assert.AreEqual(y.Length, preds.Length);

            for (int i = 0; i < y.Length; i++)
            {
                Assert.IsTrue(
                    Math.Abs(preds[i] - y[i]) < 1e-10,
                    $"Prediction mismatch at {i}: {preds[i]} vs {y[i]}"
                );
            }
        }

        [TestMethod]
        public void LogisticRegression_RollingCV_PerfectSeparation_ShouldHaveHighScore()
        {
            
            double[,] Xdata =
            {
        { 0 },
        { 1 },
        { 2 },
        { 3 },
        { 4 },
        { 5 },
        { 6 },
        { 7 },
        { 8 },
        { 9 }
    };

          
            double[] ydata =
            {
        0, 0, 0, 0, 0,
        1, 1, 1, 1, 1
    };

            var X = new Matrix(Xdata);
            var y = new VectorN(ydata);

            var model = new Logistic(fitIntercept: true);

            var modelParams = new Dictionary<string, object>
    {
        { "LearningRate", 0.1 },
        { "MaxIterations", 2000 }
    };

            var pipeline = new Pipeline(
                model: model,
                modelParams: modelParams
            );

            var cv = new RollingCrossValidator(
                pipelines: new List<Pipeline> { pipeline },
                folds: 5
            );

         
            var result = cv.Run(X, y);

            Assert.IsNotNull(result.BestPipeline);

          
            Assert.IsTrue(result.BestScore > -0.05,
                $"Expected high classification score, got {result.BestScore}");
        }

        [TestMethod]
        public void StandardScaler_FitTransform_ShouldZeroMeanAndUnitVariance()
        {
            
            double[,] data =
            {
                { 1, 2 },
                { 3, 4 },
                { 5, 6 }
            };

            var X = new Matrix(data);
            var scaler = new StandardScaler();

        
            var Xscaled = scaler.FitTransform(X);

     
            int rows = Xscaled.rowLength;
            int cols = Xscaled.columnLength;

            for (int j = 0; j < cols; j++)
            {
                double mean = 0.0;
                for (int i = 0; i < rows; i++)
                    mean += Xscaled.values[i, j];
                mean /= rows;

                double variance = 0.0;
                for (int i = 0; i < rows; i++)
                    variance += Math.Pow(Xscaled.values[i, j] - mean, 2);
                variance /= rows;

                double std = Math.Sqrt(variance);

                Assert.IsTrue(Math.Abs(mean) < 1e-10, $"Mean col {j} = {mean}");
                Assert.IsTrue(Math.Abs(std - 1.0) < 1e-10, $"Std col {j} = {std}");
            }
        }
    }
}

