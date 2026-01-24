using CSharpNumerics.ML.CrossValidators.Grouping;
using CSharpNumerics.ML.CrossValidators.Result;
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Objects;
using Numerics.Models;
using System;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;


namespace CSharpNumerics.ML.CrossValidators;

public class TimeSeriesCrossValidator
{
        public List<Pipeline> Pipelines { get; }
        public ITimeGrouping Grouping { get; }

        public TimeSeriesCrossValidator(List<Pipeline> pipelines, ITimeGrouping grouping)
        {
            Pipelines = pipelines;
            Grouping = grouping;
        }

        public RollingValidationResult Run(TimeSeries ts, string targetColumn)
        {
            
            var colIndex = Array.IndexOf(ts.Cols, targetColumn);
            if (colIndex < 0)
                throw new ArgumentException("Target column not found");

            var X = ts.ToMatrix(colIndex);
            var y = new VectorN(ts.Data[colIndex]);

      
            var groups = Grouping.GetGroups(ts.Time);
            var uniqueGroups = groups.Distinct().OrderBy(x => x).ToArray();

            var result = new RollingValidationResult();
            int totalPredictions = y.Length;
            Dictionary<Pipeline, (double[] pred, double[] actual)> cache = new(Pipelines.Count);

            foreach (var pipe in Pipelines)
            {
                var preds = ArrayPool<double>.Shared.Rent(totalPredictions);
                var actuals = ArrayPool<double>.Shared.Rent(totalPredictions);
                int predCount = 0;

                double score = 0;

                foreach (var testGroup in uniqueGroups)
                {
                    
                    var testIdx = groups
                        .Select((g, idx) => (g, idx))
                        .Where(t => t.g == testGroup)
                        .Select(t => t.idx)
                        .ToArray();

                    var trainIdx = Enumerable.Range(0, y.Length)
                        .Except(testIdx)
                        .ToArray();

                    var Xtrain = X.SliceRows(trainIdx);
                    var Xtest = X.SliceRows(testIdx);
                    var yTrain = y.Slice(trainIdx);
                    var yTest = y.Slice(testIdx);

                var cloned = new Pipeline(
                    pipe.Model, pipe.ModelParams,
                    pipe.Scaler, pipe.ScalerParams,
                    pipe.Selector, pipe.SelectorParams);

                    cloned.Fit(Xtrain, yTrain);
                    var pred = cloned.Predict(Xtest);

                    Array.Copy(pred.Values, 0, preds, predCount, pred.Length);
                    Array.Copy(yTest.Values, 0, actuals, predCount, yTest.Length);
                    predCount += pred.Length;

                    score += EvaluateScore(pipe.Model is IRegressionModel, pred, yTest) * yTest.Length;
                }

                score /= y.Length;
                result.Scores[pipe] = score;

                var actualPreds = new double[predCount];
                var actualActuals = new double[predCount];
                Array.Copy(preds, 0, actualPreds, 0, predCount);
                Array.Copy(actuals, 0, actualActuals, 0, predCount);
                cache[pipe] = (actualPreds, actualActuals);

                ArrayPool<double>.Shared.Return(preds);
                ArrayPool<double>.Shared.Return(actuals);
            }

            result.BestPipeline = result.Scores.OrderByDescending(x => x.Value).First().Key;
            result.BestScore = result.Scores[result.BestPipeline];

            var bestData = cache[result.BestPipeline];
            result.ActualValues = new VectorN(bestData.actual);
            result.PredictedValues = new VectorN(bestData.pred);

      
            var dataPairs = result.ActualValues.Values
                .Select((v, i) => (Pred: result.PredictedValues.Values[i], Actual: v));
            result.CoefficientOfDetermination = dataPairs.CoefficientOfDetermination(p => (p.Pred, p.Actual));

        if (!(result.BestPipeline.Model is IRegressionModel))
        {
            result.ConfusionMatrix =
                result.ActualValues.BuildConfusionMatrix(result.PredictedValues);
        }

        return result;
        }

        private double EvaluateScore(bool isRegression, VectorN pred, VectorN y)
        {
            if (isRegression)
            {
                double sum = 0;
                var predVals = pred.Values;
                var yVals = y.Values;
                for (int i = 0; i < pred.Length; i++)
                {
                    double diff = predVals[i] - yVals[i];
                    sum += diff * diff;
                }
                return -sum / pred.Length;
            }
            else
            {
                int correct = 0;
                var predVals = pred.Values;
                var yVals = y.Values;
                for (int i = 0; i < pred.Length; i++)
                    if (Math.Round(predVals[i]) == yVals[i]) correct++;
                return (double)correct / pred.Length;
            }
        }
}
