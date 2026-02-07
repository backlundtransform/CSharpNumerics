using CSharpNumerics.ML;
using CSharpNumerics.ML.CrossValidators.Grouping;
using CSharpNumerics.ML.CrossValidators.Interfaces;
using CSharpNumerics.ML.CrossValidators.Result;
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Models;
using CSharpNumerics.Objects;
using Numerics.Models;
using Numerics.Objects;
using System;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;
namespace CSharpNumerics.ML.CrossValidators;

public class RollingCrossValidator: ICrossValidator, ITimeSeriesCrossValidator, ISeriesCrossValidator
{
    public List<Pipeline> Pipelines { get; }
   
    public RollingCrossValidator(List<Pipeline> pipelines)
    {
        Pipelines = pipelines;
  
    }

    public RollingCrossValidator(PipelineGrid pipelineGrid)
    {
        Pipelines = [.. pipelineGrid.Expand()];
      
    }

    public CrossValidationResult Run(
    TimeSeries ts,
    string targetColumn,
    ITimeGrouping grouping = null)
    {
        var colIndex = Array.IndexOf(ts.Cols, targetColumn);
        if (colIndex < 0)
            throw new ArgumentException("Target column not found");

        var X = ts.ToMatrix(colIndex);
        var y = new VectorN(ts.Data[colIndex]);

        var groups = grouping != null
            ? grouping.GetGroups(ts.Time)
            : Enumerable.Range(0, y.Length).ToArray();

        return RunInternal(X, y, groups);
    }

    public CrossValidationResult Run(Series serie, string targetColumn)
    {

        var colIndex = Array.IndexOf(serie.Cols, targetColumn);
        if (colIndex < 0)
            throw new ArgumentException("Target column not found");

        var X = serie.ToMatrix(colIndex);
        var y = new VectorN(serie.Data[colIndex]);


        var groups = serie.Groups ?? Enumerable.Range(0, serie.Data.Length).ToArray();

        return RunInternal(X, y, groups);


    }

    public CrossValidationResult Run(Matrix X, VectorN y)
    {
        var groups = Enumerable.Range(0, y.Length).ToArray();
        return RunInternal(X, y, groups);
    }

    private CrossValidationResult RunInternal(
    Matrix X,
    VectorN y,
    int[] groups)
    {

        var result = new CrossValidationResult();
        int totalPredictions = y.Length;

        Dictionary<Pipeline, (double[] pred, double[] actual)> cache = new(Pipelines.Count);

        foreach (var pipe in Pipelines)
        {
            var preds = ArrayPool<double>.Shared.Rent(totalPredictions);
            var actuals = ArrayPool<double>.Shared.Rent(totalPredictions);
            int predCount = 0;

            double score = EvaluatePipelineRolling(pipe, X, y, groups, preds, actuals, ref predCount);

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

    private double EvaluatePipelineRolling(
       Pipeline pipe,
       Matrix X,
       VectorN y,
       int[] groups,
       double[] allPred,
       double[] allActual,
       ref int predCount)
    {
        var uniqueGroups = groups
            .Distinct()
            .OrderBy(g => g)
            .ToArray();

        if (uniqueGroups.Length < 2)
            throw new InvalidOperationException(
                "Rolling CV requires at least 2 groups");

        double totalScore = 0;
        int totalItems = 0;

        for (int i = 1; i < uniqueGroups.Length; i++)
        {
            var trainGroups = uniqueGroups.Take(i).ToHashSet();
            int valGroup = uniqueGroups[i];

            var trainIdx = groups
                .Select((g, idx) => (g, idx))
                .Where(t => trainGroups.Contains(t.g))
                .Select(t => t.idx)
                .ToArray();

            var valIdx = groups
                .Select((g, idx) => (g, idx))
                .Where(t => t.g == valGroup)
                .Select(t => t.idx)
                .ToArray();

            if (valIdx.Length == 0)
                continue;

            var Xtrain = X.SubMatrix(trainIdx);
            var yTrain = y.SubVector(trainIdx);
            var Xval = X.SubMatrix(valIdx);
            var yVal = y.SubVector(valIdx);

            var cloned = pipe.Clone();

            cloned.Fit(Xtrain, yTrain);
            var pred = cloned.Predict(Xval);

            Array.Copy(pred.Values, 0, allPred, predCount, pred.Length);
            Array.Copy(yVal.Values, 0, allActual, predCount, yVal.Length);
            predCount += pred.Length;

            double foldScore =
                ScoreOptimized(pipe.Model is IRegressionModel, pred, yVal);

            totalScore += foldScore * yVal.Length;
            totalItems += yVal.Length;
        }

        return totalItems > 0
            ? totalScore / totalItems
            : double.NaN;
    }

    private double ScoreOptimized(bool isRegression, VectorN pred, VectorN y)
    {
        if (isRegression)
        {
            // MSE - optimized loop
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
            // Accuracy
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


