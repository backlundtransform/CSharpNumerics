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

namespace CSharpNumerics.ML.CrossValidators;

public class StratifiedKFoldCrossValidator : ICrossValidator
{
    public List<Pipeline> Pipelines { get; }
    public int Folds { get; }

    public StratifiedKFoldCrossValidator(List<Pipeline> pipelines, int folds = 5)
    {
        Pipelines = pipelines;
        Folds = folds;
    }

    public StratifiedKFoldCrossValidator(PipelineGrid pipelineGrid, int folds = 5)
    {
        Pipelines = [.. pipelineGrid.Expand()];
        Folds = folds;
    }

    public CrossValidationResult Run(Matrix X, VectorN y)
    {
        var result = new CrossValidationResult();
        int n = y.Length;
        Dictionary<Pipeline, (double[] pred, double[] actual)> cache = new(Pipelines.Count);

  
        var classIndices = y.Values
            .Select((v, i) => (Value: v, Index: i))
            .GroupBy(p => p.Value)
            .ToDictionary(g => g.Key, g => g.Select(p => p.Index).ToList());

        foreach (var pipe in Pipelines)
        {

            if (!(pipe.Model is IClassificationModel))
                throw new InvalidOperationException(
                    $"Pipeline {pipe.Model.GetType().Name} is not a classification model. " +
                    "StratifiedKFoldCrossValidator can only be used with classification models."
                );

            var preds = ArrayPool<double>.Shared.Rent(n);
            var actuals = ArrayPool<double>.Shared.Rent(n);
            int predCount = 0;
            double totalScore = 0;

            var foldIndices = new List<List<int>>();
            for (int f = 0; f < Folds; f++)
                foldIndices.Add(new List<int>());

            foreach (var indices in classIndices.Values)
            {
                int i = 0;
                foreach (var idx in indices)
                {
                    foldIndices[i % Folds].Add(idx);
                    i++;
                }
            }

           
            for (int fold = 0; fold < Folds; fold++)
            {
                var testIdx = foldIndices[fold].OrderBy(x => x).ToArray();
                var trainIdx = Enumerable.Range(0, n).Except(testIdx).ToArray();

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

            result.Scores[pipe] = totalScore / n;

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
