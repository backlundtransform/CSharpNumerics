using CSharpNumerics.ML;
using CSharpNumerics.ML.Models.Classification;
using CSharpNumerics.ML.Models.Regression;
using CSharpNumerics.ML.Scalers;
using CSharpNumerics.ML.Selector;
using CSharpNumerics.Objects;
using Numerics;
using Numerics.Objects;
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

            var model = new Linear();

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

            var model = new Linear();

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

            var model = new Logistic();

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
        [TestMethod]
        public void Fit_ComputeMeansAndStdDevs_Correctly()
        {
            // Arrange
            double[,] data = {
            { 1.0, 2.0, 3.0 },
            { 4.0, 5.0, 6.0 },
            { 7.0, 8.0, 9.0 }
        };
            StandardScaler scaler = new StandardScaler();

            // Act
            scaler.Fit(data);

            // Assert
            Assert.IsTrue(4.0 == scaler.Means[0]);
            Assert.IsTrue(5.0 == scaler.Means[1]);
            Assert.IsTrue(6.0 == scaler.Means[2]);
            Assert.IsTrue(2.4495 == Math.Round(scaler.StdDevs[0], 4));
            Assert.IsTrue(2.4495 == Math.Round(scaler.StdDevs[1], 4));
            Assert.IsTrue(2.4495 == Math.Round(scaler.StdDevs[2], 4));
        }
        [TestMethod]
        public void Transform_StandardizeData_Correctly()
        {
            // Arrange
            double[,] data = {
            { 1.0, 2.0, 3.0 },
            { 4.0, 5.0, 6.0 },
            { 7.0, 8.0, 9.0 }
        };
            StandardScaler scaler = new StandardScaler();
            scaler.Fit(data);

            // Act
            Matrix transformedDataMatrix = scaler.FitTransform(new Matrix(data));

            var transformedData = transformedDataMatrix.values;

            // Assert
            Assert.IsTrue(-1.2247 == Math.Round(transformedData[0, 0], 4));
            Assert.IsTrue(-1.2247 == Math.Round(transformedData[0, 1], 4));
            Assert.IsTrue(-1.2247 == Math.Round(transformedData[0, 2], 4));
            Assert.IsTrue(0.0 == Math.Round(transformedData[1, 0], 1));
            Assert.IsTrue(0.0 == Math.Round(transformedData[1, 1], 1));
            Assert.IsTrue(0.0 == Math.Round(transformedData[1, 2], 1));
            Assert.IsTrue(1.2247 == Math.Round(transformedData[2, 0], 4));
            Assert.IsTrue(1.2247 == Math.Round(transformedData[2, 1], 4));
            Assert.IsTrue(1.2247 == Math.Round(transformedData[2, 2], 4));
        }

        [TestMethod]
        public void FitTransform_StandardizeData_Correctly()
        {
            // Arrange
            double[,] data = {
            { 1.0, 2.0, 3.0 },
            { 4.0, 5.0, 6.0 },
            { 7.0, 8.0, 9.0 }
        };
            StandardScaler scaler = new StandardScaler();

            // Act
            Matrix transformedDataMatrix = scaler.FitTransform(new Matrix(data));

            var transformedData = transformedDataMatrix.values;

            // Assert
            Assert.IsTrue(-1.2247 == Math.Round(transformedData[0, 0], 4));
            Assert.IsTrue(-1.2247 == Math.Round(transformedData[0, 1], 4));
            Assert.IsTrue(-1.2247 == Math.Round(transformedData[0, 2], 4));
            Assert.IsTrue(0.0 == Math.Round(transformedData[1, 0], 1));
            Assert.IsTrue(0.0 == Math.Round(transformedData[1, 1], 1));
            Assert.IsTrue(0.0 == Math.Round(transformedData[1, 2], 1));
            Assert.IsTrue(1.2247 == Math.Round(transformedData[2, 0], 4));
            Assert.IsTrue(1.2247 == Math.Round(transformedData[2, 1], 4));
            Assert.IsTrue(1.2247 == Math.Round(transformedData[2, 2], 4));
        }
        [TestMethod]
        public void TestCrossValidatorPipeline()
        {

            int nSamples = 100;
            int nFeatures = 2;

            Matrix X = new Matrix(nSamples, nFeatures);
            VectorN y = new VectorN(nSamples);

            Random rnd = new Random(123);
            for (int i = 0; i < nSamples; i++)
            {
                double x1 = rnd.NextDouble() * 10;
                double x2 = rnd.NextDouble() * 10;
                X.values[i, 0] = x1;
                X.values[i, 1] = x2;
                y[i] = 3 * x1 + 2 * x2 + rnd.NextDouble() * 0.1;
            }

            var pipeline = new Pipeline(
                new Linear(),
                new Dictionary<string, object> { ["LearningRate"] = 0.01 },
                selector: new SelectKBest(),
                selectorParams: (Dictionary<string, object>)new Dictionary<string, object> { ["K"] = 1 }
            );

            var cv = new RollingCrossValidator([pipeline], 5);
            var results = cv.Run(X, y);


            // 5. Skriv ut resultat
            Console.WriteLine("Cross validation results:");


        }

        [TestMethod]
        public void TestRollingCrossValidator_DecisionTree_Multiclass()
        {
            // Arrange
            int nSamples = 120;
            int nFeatures = 2;

            Matrix X = new Matrix(nSamples, nFeatures);
            VectorN y = new VectorN(nSamples);

            Random rnd = new Random(123);

            for (int i = 0; i < nSamples; i++)
            {
                double x1 = rnd.NextDouble() * 10;
                double x2 = rnd.NextDouble() * 10;

                X.values[i, 0] = x1;
                X.values[i, 1] = x2;

                double s = x1 + x2;

                if (s < 7)
                    y[i] = 0;
                else if (s < 13)
                    y[i] = 1;
                else
                    y[i] = 2;
            }

            var tree = new DecisionTree();
            tree.SetHyperParameters(new()
            {
                ["MaxDepth"] = 5,
                ["MinSamplesSplit"] = 2
            });

            var pipeline = new Pipeline(
                model: tree,
                modelParams: (Dictionary<string, object>)new Dictionary<string, object>
                {
                    ["MaxDepth"] = 5,
                    ["MinSamplesSplit"] = 2
                },
                scaler: null,
                scalerParams: null,
                selector: null,
                selectorParams: null
            );

            var cv = new RollingCrossValidator(new List<Pipeline> { pipeline }, folds: 5);

            // Act
            var result = cv.Run(X, y);

            // Assert
            Assert.IsNotNull(result.BestPipeline);
            Assert.IsTrue(result.BestScore > 0.8, $"Accuracy too low: {result.BestScore}");

            Assert.IsNotNull(result.ConfusionMatrix);
            Assert.AreEqual(3, result.ConfusionMatrix.rowLength);
            Assert.AreEqual(3, result.ConfusionMatrix.columnLength);

            // Sanity check: diagonal dominance
            double diag = 0;
            double total = 0;

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    double v = result.ConfusionMatrix.values[i, j];
                    total += v;
                    if (i == j) diag += v;
                }
            }

            Assert.IsTrue(diag / total > 0.8);
        }



        [TestMethod]
        public void TestRollingCrossValidator_Classification()
        {
            // Arrange
            int nSamples = 120;
            int nFeatures = 2;

            Matrix X = new Matrix(nSamples, nFeatures);
            VectorN y = new VectorN(nSamples);

            Random rnd = new Random(123);

            for (int i = 0; i < nSamples; i++)
            {
                double x1 = rnd.NextDouble() * 10;
                double x2 = rnd.NextDouble() * 10;

                X.values[i, 0] = x1;
                X.values[i, 1] = x2;

                double s = x1 + x2;

                if (s < 7)
                    y[i] = 0;
                else if (s < 13)
                    y[i] = 1;
                else
                    y[i] = 2;
            }

            var tree = new DecisionTree();


            var pipeline1 = new Pipeline(
                model: tree,
                modelParams: (Dictionary<string, object>)new Dictionary<string, object>
                {
                    ["MaxDepth"] = 5,
                    ["MinSamplesSplit"] = 2
                }

            );

            var rf = new RandomForest();


            var pipeline2 = new Pipeline(
                rf,
                (Dictionary<string, object>)new Dictionary<string, object>
                {
                    ["NumTrees"] = 30,
                    ["MaxDepth"] = 6
                },
                null,
                null,
                null,
                null
            );




            var pipeline3 = new Pipeline(
                new Logistic(),
                (Dictionary<string, object>)new Dictionary<string, object>
                {
                    ["LearningRate"] = 0.1,
                    ["MaxIterations"] = 2000
                },
                null,
                null,
                null,
                null
            );

            var cv = new RollingCrossValidator(
    pipelines: new List<Pipeline> { pipeline1, pipeline2, pipeline3 },
    folds: 5
);

            // Act
            var result = cv.Run(X, y);

            // Assert
            Assert.IsNotNull(result.BestPipeline);
            Assert.IsTrue(result.BestScore > 0.8, $"Accuracy too low: {result.BestScore}");

            Assert.IsNotNull(result.ConfusionMatrix);
            Assert.AreEqual(3, result.ConfusionMatrix.rowLength);
            Assert.AreEqual(3, result.ConfusionMatrix.columnLength);

            // Sanity check: diagonal dominance
            double diag = 0;
            double total = 0;

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    double v = result.ConfusionMatrix.values[i, j];
                    total += v;
                    if (i == j) diag += v;
                }
            }

            Assert.IsTrue(diag / total > 0.8);
        }


        [TestMethod]
        public void TestRollingCrossValidator_Classification_Search()
        {
            // Arrange
            int nSamples = 120;
            int nFeatures = 2;

            Matrix X = new Matrix(nSamples, nFeatures);
            VectorN y = new VectorN(nSamples);

            Random rnd = new Random(123);

            for (int i = 0; i < nSamples; i++)
            {
                double x1 = rnd.NextDouble() * 10;
                double x2 = rnd.NextDouble() * 10;

                X.values[i, 0] = x1;
                X.values[i, 1] = x2;

                double s = x1 + x2;

                if (s < 7)
                    y[i] = 0;
                else if (s < 13)
                    y[i] = 1;
                else
                    y[i] = 2;
            }
            var pipelineGrid =
        new PipelineGrid()
            .AddModel<RandomForest>(g => g
                .Add("NumTrees", 50, 100, 200)
                .Add("MaxDepth", 5, 8, 10))

            .AddModel<Logistic>(g => g
                .Add("LearningRate", 0.05, 0.1)
                .Add("MaxIterations", 1000, 2000)
                .AddScaler<StandardScaler>(s => { })
                .AddSelector<SelectKBest>(s => s
                   .Add("K", 1, 2)))
            .AddModel<DecisionTree>(g => g
                .Add("MaxDepth", 3, 5, 8))
            .AddModel<KNearestNeighbors>(g => g
                .Add("K", 3, 5, 7))
            .AddModel<LinearSVC>(g => g
                .Add("C", 0.1, 1.0, 10.0)
                .Add("LearningRate", 0.001, 0.01)
                .Add("Epochs", 500, 1000));

            var cv = new RollingCrossValidator(pipelineGrid, folds: 5);

            // Act
            var result = cv.Run(X, y);

            // Assert
            Assert.IsNotNull(result.BestPipeline);
            Assert.IsTrue(result.BestScore > 0.8, $"Accuracy too low: {result.BestScore}");

            Assert.IsNotNull(result.ConfusionMatrix);
            Assert.AreEqual(3, result.ConfusionMatrix.rowLength);
            Assert.AreEqual(3, result.ConfusionMatrix.columnLength);

            // Sanity check: diagonal dominance
            double diag = 0;
            double total = 0;

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    double v = result.ConfusionMatrix.values[i, j];
                    total += v;
                    if (i == j) diag += v;
                }
            }

            Assert.IsTrue(diag / total > 0.8);
        }

        [TestMethod]
        public void TestRollingCrossValidator_Classification_SVM()
        {
            // Arrange
            int nSamples = 120;
            int nFeatures = 2;

            Matrix X = new Matrix(nSamples, nFeatures);
            VectorN y = new VectorN(nSamples);

            Random rnd = new Random(123);

            for (int i = 0; i < nSamples; i++)
            {
                double x1 = rnd.NextDouble() * 10;
                double x2 = rnd.NextDouble() * 10;

                X.values[i, 0] = x1;
                X.values[i, 1] = x2;

                double s = x1 + x2;

                if (s < 7)
                    y[i] = 0;
                else if (s < 13)
                    y[i] = 1;
                else
                    y[i] = 2;
            }
            var pipelineGrid =
        new PipelineGrid()

            .AddModel<LinearSVC>(g => g
                .Add("C", 0.1, 1.0, 10.0)
                .Add("LearningRate", 0.001, 0.01)
                .Add("Epochs", 500, 1000));

            var cv = new RollingCrossValidator(pipelineGrid, folds: 5);

            // Act
            var result = cv.Run(X, y);

            // Assert
            Assert.IsNotNull(result.BestPipeline);
            Assert.IsTrue(result.BestScore > 0.7, $"Accuracy too low: {result.BestScore}");

            Assert.IsNotNull(result.ConfusionMatrix);
            Assert.AreEqual(3, result.ConfusionMatrix.rowLength);
            Assert.AreEqual(3, result.ConfusionMatrix.columnLength);

            // Sanity check: diagonal dominance
            double diag = 0;
            double total = 0;

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    double v = result.ConfusionMatrix.values[i, j];
                    total += v;
                    if (i == j) diag += v;
                }
            }

            Assert.IsTrue(diag / total > 0.7);


            pipelineGrid =
                new PipelineGrid()
                .AddModel<KernelSVC>(g => g
                    .Add("C", 1.0)
                    .Add("Gamma", 0.1)
                    .Add("Kernel", KernelType.RBF));

            cv = new RollingCrossValidator(pipelineGrid, folds: 5);
            // Act
            result = cv.Run(X, y);
            // Assert
            Assert.IsNotNull(result.BestPipeline);
            Assert.IsTrue(result.BestScore > 0.8, $"Accuracy too low: {result.BestScore}");

        }



        [TestMethod]
        public void TestRollingCrossValidator_Regression_Search()
        {
            // Arrange: y = 2x + 1
            double[,] Xdata =
            {
        { 0 }, { 1 }, { 2 }, { 3 }, { 4 },
        { 5 }, { 6 }, { 7 }, { 8 }, { 9 }
    };

            double[] ydata = { 1, 3, 5, 7, 9, 11, 13, 15, 17, 19 };

            var X = new Matrix(Xdata);
            var y = new VectorN(ydata);

            var pipelineGrid =
                new PipelineGrid()
                    .AddModel<Linear>(g => g
                        .Add("FitIntercept", true)
                        .Add("LearningRate", 0.05, 0.1)
                        .AddScaler<StandardScaler>(s => { }))
                    .AddModel<Ridge>(g => g
                        .Add("Alpha", 0.001, 0.01, 0.1)
                        .AddScaler<StandardScaler>(s => { }));

            var cv = new RollingCrossValidator(pipelineGrid, folds: 3);

            // Act
            var result = cv.Run(X, y);

            // Assert
            Assert.IsNotNull(result.BestPipeline);


            Assert.IsTrue(result.CoefficientOfDetermination > 0.99, $"R² too low: {result.CoefficientOfDetermination}");
        }

    [TestMethod]
        public void Test_KernelSVR_SimpleRegression()
        {
            double[,] Xdata =
            {
        {0},{1},{2},{3},{4},{5},{6},{7},{8},{9}
    };

            double[] ydata =
            {
        1,3,5,7,9,11,13,15,17,19
    };

            var X = new Matrix(Xdata);
            var y = new VectorN(ydata);

            var grid = new PipelineGrid()
                .AddModel<KernelSVR>(g => g
                    .Add("Kernel", KernelType.RBF)
                    .Add("C", 1.0)
                    .Add("Gamma", 0.5)
                    .Add("Epsilon", 0.01)
                    .AddScaler<StandardScaler>(s => { }));

            var cv = new RollingCrossValidator(grid, folds: 3);
            var result = cv.Run(X, y);

            double r2 = result.CoefficientOfDetermination;
            Assert.IsTrue(r2 > 0.3);
        }
      
     [TestMethod]
        public void Test_MPL_SimpleRegression()
        {
            double[,] Xdata =
            {
        {0},{1},{2},{3},{4},{5},{6},{7},{8},{9}
    };

            double[] ydata =
            {
        1,3,5,7,9,11,13,15,17,19
    };

            var X = new Matrix(Xdata);
            var y = new VectorN(ydata);

            var grid = new PipelineGrid()
              .AddModel<MLPRegressor>(g => g
    .Add("HiddenLayers", new[] { 32, 16 }, new[] { 64, 32 })
    .Add("LearningRate", 0.001, 0.01)
    .Add("Epochs", 500, 1000)
    .Add("L2", 0.0, 0.001)
    .Add("ValidationSplit", 0.0)
    .Add("Patience", 100)
    .AddScaler<StandardScaler>(s => { }));

            var cv = new RollingCrossValidator(grid, folds: 3);
            var result = cv.Run(X, y);

            double r2 = result.CoefficientOfDetermination;
            Assert.IsTrue(r2 > 0.9);
        }
    
    }
}

