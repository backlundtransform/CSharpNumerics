using CSharpNumerics.ML;
using System;
using System.Linq;
using CSharpNumerics.ML.CrossValidators;
using CSharpNumerics.ML.Models.Classification;
using CSharpNumerics.ML.Models.Regression;
using CSharpNumerics.ML.Scalers;

using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Data;
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
        public void SeriesFromCsv_ColsMustAlignWithData()
        {
            // Regression test for the off-by-one that made every Series-based
            // cross-validation train on a leaked target and validate against a
            // feature: Cols was header.Skip(1), so IndexOf(Cols, name) pointed
            // one column left of the truth in Data.
            CsvTestDataGenerator.GenerateClassificationCsv("cols_alignment.csv");
            var df = Series.FromCsv("cols_alignment.csv");

            Assert.AreEqual(df.Data.Length, df.Cols.Length,
                "every data column must have a name");

            var target = df.Data[Array.IndexOf(df.Cols, "Target")];
            Assert.IsTrue(target.All(v => v == 0.0 || v == 1.0),
                "the column Cols calls 'Target' must hold the class labels, not a feature");
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

            var result = cv.Run(df.ToMatrix(colIndex), new VectorN(df.Data[colIndex]));


            Assert.IsTrue(result.BestScore > 0);
        }

    }


}
