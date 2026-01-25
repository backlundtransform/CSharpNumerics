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

public class KFoldCrossValidator: ICrossValidator
{
    public List<Pipeline> Pipelines { get; }
    public int Folds { get; }

    public KFoldCrossValidator(List<Pipeline> pipelines, int folds = 5)
    {
        Pipelines = pipelines;
        Folds = folds;
    }

    public KFoldCrossValidator(PipelineGrid pipelineGrid, int folds = 5)
    {
        Pipelines = [.. pipelineGrid.Expand()];
        Folds = folds;
    }

    public CrossValidationResult Run(Matrix X, VectorN y)
    {
        var result = new CrossValidationResult();
        int totalPredictions = y.Length;

        Dictionary<Pipeline, (double[] pred, double[] actual)> cache = new(Pipelines.Count);

        foreach (var pipe in Pipelines)
        {
            var preds = ArrayPool<double>.Shared.Rent(totalPredictions);
            var actuals = ArrayPool<double>.Shared.Rent(totalPredictions);
            int predCount = 0;

            double score = EvaluatePipelineRolling(pipe, X, y, preds, actuals, ref predCount);

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
      double[] allPred,
      double[] allActual,
      ref int predCount)
    {
        int n = y.Length;
        int foldSize = n / Folds;

        double totalScore = 0;
        int totalItems = 0;

        for (int fold = 0; fold < Folds; fold++)
        {
            int start = fold * foldSize;
            int end = (fold == Folds - 1) ? n : start + foldSize;

            var (Xtrain, Ytrain, Xval, Yval) = SplitOptimized(X, y, start, end);
            if (Yval.Length == 0)
                continue;

            var cloned = new Pipeline(
                pipe.Model, pipe.ModelParams,
                pipe.Scaler, pipe.ScalerParams,
                pipe.Selector, pipe.SelectorParams);

            cloned.Fit(Xtrain, Ytrain);
            var pred = cloned.Predict(Xval);


            Array.Copy(pred.Values, 0, allPred, predCount, pred.Length);
            Array.Copy(Yval.Values, 0, allActual, predCount, Yval.Length);
            predCount += pred.Length;

            double foldScore = ScoreOptimized(pipe.Model is IRegressionModel, pred, Yval);

            totalScore += foldScore * (end - start);
            totalItems += (end - start);
        }

        return totalScore / totalItems;
    }

    private (Matrix, VectorN, Matrix, VectorN) SplitOptimized(Matrix X, VectorN y, int start, int end)
    {
        int n = X.rowLength;
        int valSize = end - start;
        int trainSize = n - valSize;

        double[,] Xt = new double[trainSize, X.columnLength];
        double[,] Xv = new double[valSize, X.columnLength];
        double[] yt = new double[trainSize];
        double[] yv = new double[valSize];

        int ti = 0, vi = 0;

        for (int i = 0; i < start; i++)
        {
            for (int j = 0; j < X.columnLength; j++)
                Xt[ti, j] = X.values[i, j];
            yt[ti] = y.Values[i];
            ti++;
        }

        for (int i = start; i < end; i++)
        {
            for (int j = 0; j < X.columnLength; j++)
                Xv[vi, j] = X.values[i, j];
            yv[vi] = y.Values[i];
            vi++;
        }

        for (int i = end; i < n; i++)
        {
            for (int j = 0; j < X.columnLength; j++)
                Xt[ti, j] = X.values[i, j];
            yt[ti] = y.Values[i];
            ti++;
        }

        return (new Matrix(Xt), new VectorN(yt), new Matrix(Xv), new VectorN(yv));
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