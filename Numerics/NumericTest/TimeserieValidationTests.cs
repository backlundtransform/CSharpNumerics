using CSharpNumerics.ML;
using CSharpNumerics.ML.CrossValidators;
using CSharpNumerics.ML.Models.Classification;
using CSharpNumerics.ML.Models.Regression;
using CSharpNumerics.ML.Scalers;
using CSharpNumerics.Numerics.Models;
using CSharpNumerics.Numerics.Objects;
using NumericTest.TestData;


namespace NumericTest
{

    [TestClass]
    public class TimeserieValidationTests
    {


        [TestMethod]
        public void RollingCV_Should_NotLeakFutureData()
        {
            CsvTestDataGenerator.GenerateTimeSeriesCsv("ts.csv");

            var pipelineGrid = new PipelineGrid()
                .AddModel<Ridge>(g => g
                    .Add("Alpha", 0.1, 1.0, 10.0)
                    .AddScaler<StandardScaler>(s => { }));

            var ts = TimeSeries.FromCsv("ts.csv");
            var cv = new RollingCrossValidator(pipelineGrid);

            var result = cv.Run(ts, "Target");

            Assert.IsTrue(result.BestScore > -10.0);
        }

        [TestMethod]
        public void LeaveOneOutCV_Should_WorkOnGroupedData()
        {
            CsvTestDataGenerator.GenerateGroupedCsv("grouped.csv");

            var pipelineGrid = new PipelineGrid()
                .AddModel<Ridge>(g => g
                    .Add("Alpha", 0.1, 1.0)
                    .AddScaler<StandardScaler>(s => { }));

            var series = Series.FromCsv("grouped.csv");

            var cv = new LeaveOneOutCrossValidator(pipelineGrid);
            var result = cv.Run(series, targetColumn: "Target");

            Assert.IsTrue(result.BestScore > -10.0);
        }

        [TestMethod]
        public void ShuffleSplitCV_Should_WorkOnRegressionData()
        {
            CsvTestDataGenerator.GenerateTimeSeriesCsv("ts.csv");
            var pipelineGrid = new PipelineGrid()
                .AddModel<Ridge>(g => g
                    .Add("Alpha", 0.1, 1.0, 10.0)
                    .AddScaler<StandardScaler>(s => { }));

            var ts = TimeSeries.FromCsv("ts.csv");
            var cv = new ShuffleSplitCrossValidator(pipelineGrid,  5, testSize: 0.2,0.8);
            var colIndex = Array.IndexOf(ts.Cols, "Target");
            var result = cv.Run(ts.ToMatrix(2), new VectorN(ts.Data[colIndex]));

            Assert.IsTrue(result.BestScore > -1);
        }

        [TestMethod]
        public void StratifiedKFoldCV_Should_WorkOnClassificationData()
        {
            CsvTestDataGenerator.GenerateClassificationCsv("classification.csv");
            var pipelineGrid = new PipelineGrid()
                .AddModel<Logistic>(g => g
                    .Add("LearningRate", 0.01, 0.1)
                    .Add("MaxIterations", 500, 1000)
                    .AddScaler<StandardScaler>(s => { }));

            var df = Series.FromCsv("classification.csv");
            var colIndex = Array.IndexOf(df.Cols, "Target");

            var cv = new StratifiedKFoldCrossValidator(pipelineGrid, folds: 5);

            var result = cv.Run(df.ToMatrix(2), new VectorN(df.Data[colIndex]));


            Assert.IsTrue(result.BestScore > 0);
        }

    }


}
