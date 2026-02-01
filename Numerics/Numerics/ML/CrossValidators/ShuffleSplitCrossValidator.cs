using CSharpNumerics.ML;
using CSharpNumerics.ML.CrossValidators.Interfaces;
using CSharpNumerics.ML.CrossValidators.Result;
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Objects;
using Numerics.Objects;
using System;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.CrossValidators
{
    public class ShuffleSplitCrossValidator : ICrossValidator
    {
        public List<Pipeline> Pipelines { get; }
        public int N_Splits { get; }
        public double TestSize { get; }
        public double TrainSize { get; }
        public int? RandomState { get; }

        public ShuffleSplitCrossValidator(
            List<Pipeline> pipelines,
            int n_splits = 5,
            double testSize = 0.2,
            double trainSize = 0.8,
            int? randomState = null)
        {
            Pipelines = pipelines;
            N_Splits = n_splits;
            TestSize = testSize;
            TrainSize = trainSize;
            RandomState = randomState;
        }

        public ShuffleSplitCrossValidator(
            PipelineGrid pipelineGrid,
            int n_splits = 5,
            double testSize = 0.2,
            double trainSize = 0.8,
            int? randomState = null)
        {
            Pipelines = pipelineGrid.Expand().ToList();
            N_Splits = n_splits;
            TestSize = testSize;
            TrainSize = trainSize;
            RandomState = randomState;
        }

        public CrossValidationResult Run(Matrix X, VectorN y)
        {
            var result = new CrossValidationResult();
            int n = y.Length;
            Random rng = RandomState.HasValue ? new Random(RandomState.Value) : new Random();
            Dictionary<Pipeline, (double[] pred, double[] actual)> cache = new(Pipelines.Count);

            foreach (var pipe in Pipelines)
            {
                var preds = ArrayPool<double>.Shared.Rent(n);
                var actuals = ArrayPool<double>.Shared.Rent(n);
                int predCount = 0;
                double totalScore = 0;

                for (int split = 0; split < N_Splits; split++)
                {
                    // 1️⃣ Slumpa index
                    var indices = Enumerable.Range(0, n).ToArray();
                    indices = indices.OrderBy(i => rng.Next()).ToArray();

                    int testCount = (int)(TestSize * n);
                    int trainCount = (int)(TrainSize * n);

                    var testIdx = indices.Take(testCount).ToArray();
                    var trainIdx = indices.Skip(testCount).Take(trainCount).ToArray();

                    var Xtrain = X.SubMatrix(trainIdx);
                    var Ytrain = new VectorN(trainIdx.Select(i => y.Values[i]).ToArray());
                    var Xval = X.SubMatrix(testIdx);
                    var Yval = new VectorN(testIdx.Select(i => y.Values[i]).ToArray());

                    var cloned = pipe.Clone();

                    cloned.Fit(Xtrain, Ytrain);
                    var foldPred = cloned.Predict(Xval);

                    Array.Copy(foldPred.Values, 0, preds, predCount, foldPred.Length);
                    Array.Copy(Yval.Values, 0, actuals, predCount, Yval.Length);
                    predCount += foldPred.Length;

                    double foldScore = ScoreOptimized(pipe.Model is IRegressionModel, foldPred, Yval);
                    totalScore += foldScore * Yval.Length;
                }

                result.Scores[pipe] = totalScore / (predCount);

                var actualPreds = new double[predCount];
                var actualActuals = new double[predCount];
                Array.Copy(preds, 0, actualPreds, 0, predCount);
                Array.Copy(actuals, 0, actualActuals, 0, predCount);
                cache[pipe] = (actualPreds, actualActuals);

                ArrayPool<double>.Shared.Return(preds);
                ArrayPool<double>.Shared.Return(actuals);
            }

            // Best pipeline
            result.BestPipeline = result.Scores.OrderByDescending(x => x.Value).First().Key;
            result.BestScore = result.Scores[result.BestPipeline];

            var bestData = cache[result.BestPipeline];
            result.ActualValues = new VectorN(bestData.actual);
            result.PredictedValues = new VectorN(bestData.pred);

            var data = result.ActualValues.Values
                .Select((v, i) => (Pred: result.PredictedValues.Values[i], Actual: v));

            result.CoefficientOfDetermination = data.CoefficientOfDetermination(p => (p.Pred, p.Actual));

            if (!(result.BestPipeline.Model is IRegressionModel))
            {
                result.ConfusionMatrix =
                    result.ActualValues.BuildConfusionMatrix(result.PredictedValues);
            }

            return result;
        }

        private double ScoreOptimized(bool isRegression, VectorN pred, VectorN y)
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
                {
                    if (Math.Round(predVals[i]) == yVals[i])
                        correct++;
                }
                return (double)correct / pred.Length;
            }
        }
    }
}
