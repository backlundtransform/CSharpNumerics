using CSharpNumerics.ML;
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Objects;
using Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;

public class RollingCrossValidator
{
    public List<Pipeline> Pipelines { get; }
    public int Folds { get; }

    public RollingCrossValidator(List<Pipeline> pipelines, int folds = 5)
    {
        Pipelines = pipelines;
        Folds = folds;
    }

    public RollingCrossValidator(PipelineGrid pipelineGrid, int folds = 5)
    {
        Pipelines = [.. pipelineGrid.Expand()];
        Folds = folds;
    }

    public RollingValidationResult Run(Matrix X, VectorN y)
    {
        var result = new RollingValidationResult();

        Dictionary<Pipeline, (List<double> pred, List<double> actual)> cache
            = new();

        foreach (var pipe in Pipelines)
        {
            var preds = new List<double>();
            var actuals = new List<double>();

            double score = EvaluatePipelineRolling(pipe, X, y, preds, actuals);

            result.Scores[pipe] = score;
            cache[pipe] = (preds, actuals);
        }

        result.BestPipeline = result.Scores.OrderByDescending(x => x.Value).First().Key;
        result.BestScore = result.Scores[result.BestPipeline];

        var bestData = cache[result.BestPipeline];
        result.ActualValues = new VectorN(bestData.actual.ToArray());
        result.PredictedValues = new VectorN(bestData.pred.ToArray());

        if (!(result.BestPipeline.Model is IRegressionModel))
        {
            result.ConfusionMatrix =
                BuildConfusionMatrix(result.ActualValues, result.PredictedValues);
        }

        return result;
    }

    private double EvaluatePipelineRolling(
      Pipeline pipe,
      Matrix X,
      VectorN y,
      List<double> allPred,
      List<double> allActual)
    {
        int n = y.Length;
        int foldSize = n / Folds;

        double totalScore = 0;
        int totalItems = 0;

        for (int fold = 0; fold < Folds; fold++)
        {
            int start = fold * foldSize;
            int end = (fold == Folds - 1) ? n : start + foldSize;

            var (Xtrain, Ytrain, Xval, Yval) = Split(X, y, start, end);
            if (Yval.Length == 0)
                continue;

            var cloned = new Pipeline(
                pipe.Model, pipe.ModelParams,
                pipe.Scaler, pipe.ScalerParams,
                pipe.Selector, pipe.SelectorParams);

            cloned.Fit(Xtrain, Ytrain);
            var pred = cloned.Predict(Xval);

         
            for (int i = 0; i < pred.Length; i++)
            {
                allPred.Add(pred[i]);
                allActual.Add(Yval[i]);
            }

            double foldScore = Score(cloned, pred, Yval);

            totalScore += foldScore * (end - start);
            totalItems += (end - start);
        }

        return totalScore / totalItems;
    }
    private Matrix BuildConfusionMatrix(VectorN yTrue, VectorN yPred)
    {
    
        int numClasses = (int)Math.Max(
            yTrue.Values.Max(),
            yPred.Values.Max()
        ) + 1;

        var cm = new Matrix(numClasses, numClasses);

        for (int i = 0; i < yTrue.Length; i++)
        {
            int actual = (int)yTrue[i];
            int predicted = (int)Math.Round(yPred[i]);

            if (predicted < 0 || predicted >= numClasses)
                continue;

            cm.values[actual, predicted]++;
        }

        return cm;
    }
    private (Matrix, VectorN, Matrix, VectorN) Split(Matrix X, VectorN y, int start, int end)
    {
        int n = X.rowLength;
        int valSize = end - start;
        int trainSize = n - valSize;

        double[,] Xt = new double[trainSize, X.columnLength];
        double[,] Xv = new double[valSize, X.columnLength];

        double[] yt = new double[trainSize];
        double[] yv = new double[valSize];

        int ti = 0, vi = 0;

        for (int i = 0; i < n; i++)
        {
            if (i >= start && i < end)
            {
                for (int j = 0; j < X.columnLength; j++)
                    Xv[vi, j] = X.values[i, j];
                yv[vi] = y[i];
                vi++;
            }
            else
            {
                for (int j = 0; j < X.columnLength; j++)
                    Xt[ti, j] = X.values[i, j];
                yt[ti] = y[i];
                ti++;
            }
        }

        return (new Matrix(Xt), new VectorN(yt), new Matrix(Xv), new VectorN(yv));
    }

    private double Score(Pipeline pipe, VectorN pred, VectorN y)
    {
        if (pipe.Model is IRegressionModel)
        {
            // MSE
            double sum = 0;
            for (int i = 0; i < pred.Length; i++)
                sum += Math.Pow(pred[i] - y[i], 2);
            return -sum / pred.Length; 
        }
        else
        {
            // Accuracy
            int correct = 0;
            for (int i = 0; i < pred.Length; i++)
                if (Math.Round(pred[i]) == y[i])
                    correct++;
            return (double)correct / pred.Length;
        }
    }
}

public class RollingValidationResult
{
    public Dictionary<Pipeline, double> Scores { get; } = new();
    public Pipeline BestPipeline { get; set; }
    public double BestScore { get; set; }

    public Matrix ConfusionMatrix { get; set; } = new Matrix();

    public VectorN ActualValues { get; set; } = new VectorN();

    public VectorN PredictedValues { get; set; } = new VectorN();
}

