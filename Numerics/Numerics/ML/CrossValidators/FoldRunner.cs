using CSharpNumerics.ML.CrossValidators.Result;
using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.Evaluation;
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.CrossValidators;

/// <summary>
/// Runs pipelines over a fixed assignment of rows to folds, where every row is tested exactly
/// once, and scores the pooled out-of-fold predictions. Shared by the grouped cross-validators.
/// </summary>
internal static class FoldRunner
{
    public static CrossValidationResult Run(List<Pipeline> pipelines, Matrix X, VectorN y, int[] foldOfRow, int folds,
        VectorN? sampleWeights, ClassificationScoring scoring)
    {
        int n = y.Length;
        if (X.rowLength != n) throw new ArgumentException("X.Rows must match y.Length.");
        if (foldOfRow.Length != n) throw new ArgumentException("Every row needs a fold.");
        if (pipelines == null || pipelines.Count == 0) throw new ArgumentException("At least one pipeline is required.");
        if (sampleWeights is { Values: not null } w && w.Length != n)
            throw new ArgumentException($"Expected {n} sample weights, got {w.Length}.", nameof(sampleWeights));

        var testRows = Enumerable.Range(0, folds).Select(f => new List<int>()).ToArray();
        for (int i = 0; i < n; i++) testRows[foldOfRow[i]].Add(i);

        var result = new CrossValidationResult();
        var outputs = new Dictionary<Pipeline, (double[] Predicted, Matrix? Probabilities)>();
        int numClasses = (int)y.Values.Max() + 1;

        foreach (var pipe in pipelines)
        {
            bool regression = pipe.Model is IRegressionModel;
            bool probabilistic = pipe.Model is IProbabilisticClassifier;
            if (!regression && scoring != ClassificationScoring.Accuracy && !probabilistic)
                throw new InvalidOperationException(
                    $"Scoring by {scoring} needs class probabilities, but {pipe.Model.GetType().Name} does not implement IProbabilisticClassifier.");

            var predicted = new double[n];
            Matrix? probabilities = probabilistic ? new Matrix(n, numClasses) : null;

            for (int f = 0; f < folds; f++)
            {
                int[] test = testRows[f].ToArray();
                if (test.Length == 0) continue;
                var inTest = new bool[n];
                foreach (int i in test) inTest[i] = true;
                int[] train = Enumerable.Range(0, n).Where(i => !inTest[i]).ToArray();

                var model = pipe.Clone();
                model.Fit(X.SubMatrix(train), Subset(y, train),
                    sampleWeights is { Values: not null } weights ? Subset(weights, train) : (VectorN?)null);

                var Xtest = X.SubMatrix(test);
                VectorN foldPredicted = model.Predict(Xtest);
                for (int k = 0; k < test.Length; k++) predicted[test[k]] = foldPredicted[k];

                if (probabilities is { } all)
                {
                    // A fold without some class trains a model with fewer columns; the rest stay 0.
                    Matrix foldProbabilities = model.PredictProbabilities(Xtest);
                    int columns = Math.Min(foldProbabilities.columnLength, numClasses);
                    for (int k = 0; k < test.Length; k++)
                        for (int c = 0; c < columns; c++)
                            all.values[test[k], c] = foldProbabilities.values[k, c];
                }
            }

            result.Scores[pipe] = Score(regression, scoring, y, predicted, probabilities, sampleWeights);
            outputs[pipe] = (predicted, probabilities);
        }

        result.BestPipeline = result.Scores.OrderByDescending(s => s.Value).First().Key;
        result.BestScore = result.Scores[result.BestPipeline];
        var best = outputs[result.BestPipeline];
        result.ActualValues = new VectorN((double[])y.Values.Clone());
        result.PredictedValues = new VectorN(best.Predicted);
        result.PredictedProbabilities = best.Probabilities;
        result.RowIndices = Enumerable.Range(0, n).ToArray();

        var pairs = result.ActualValues.Values.Select((v, i) => (Pred: best.Predicted[i], Actual: v));
        result.CoefficientOfDetermination = pairs.CoefficientOfDetermination(p => (p.Pred, p.Actual));
        if (result.BestPipeline.Model is not IRegressionModel)
            result.ConfusionMatrix = result.ActualValues.BuildConfusionMatrix(result.PredictedValues);

        return result;
    }

    private static double Score(bool regression, ClassificationScoring scoring, VectorN y, double[] predicted,
        Matrix? probabilities, VectorN? sampleWeights)
    {
        int n = y.Length;
        double[] w = sampleWeights is { Values: not null } given ? given.Values : Enumerable.Repeat(1.0, n).ToArray();
        double total = w.Sum();

        if (regression)
        {
            double sum = 0;
            for (int i = 0; i < n; i++) sum += w[i] * (predicted[i] - y[i]) * (predicted[i] - y[i]);
            return -sum / total;
        }

        if (scoring == ClassificationScoring.Accuracy)
        {
            double correct = 0;
            for (int i = 0; i < n; i++)
                if (Math.Round(predicted[i]) == y[i]) correct += w[i];
            return correct / total;
        }

        if (probabilities is not { } p || p.columnLength != 2)
            throw new InvalidOperationException($"Scoring by {scoring} is for binary classification (labels 0 and 1).");

        VectorN positive = ClassificationMetrics.ProbabilityOf(p, 1);
        return scoring switch
        {
            ClassificationScoring.RocAuc => ClassificationMetrics.RocAuc(y, positive, sampleWeights),
            ClassificationScoring.AveragePrecision => ClassificationMetrics.AveragePrecision(y, positive, sampleWeights),
            ClassificationScoring.NegativeLogLoss => -ClassificationMetrics.LogLoss(y, positive, sampleWeights),
            ClassificationScoring.NegativeBrierScore => -ClassificationMetrics.BrierScore(y, positive, sampleWeights),
            _ => throw new ArgumentOutOfRangeException(nameof(scoring)),
        };
    }

    internal static VectorN Subset(VectorN v, int[] rows)
    {
        var values = new double[rows.Length];
        for (int k = 0; k < rows.Length; k++) values[k] = v[rows[k]];
        return new VectorN(values);
    }

    /// <summary>Group id per row → index of each distinct group, in order of first appearance.</summary>
    internal static (int[] GroupOfRow, List<long> Keys) IndexGroups(IReadOnlyList<long> groups, int rows)
    {
        if (groups == null) throw new ArgumentNullException(nameof(groups));
        if (groups.Count != rows)
            throw new ArgumentException($"Expected one group per row ({rows}), got {groups.Count}.", nameof(groups));

        var index = new Dictionary<long, int>();
        var keys = new List<long>();
        var groupOfRow = new int[rows];
        for (int i = 0; i < rows; i++)
        {
            if (!index.TryGetValue(groups[i], out int g))
            {
                g = keys.Count;
                index[groups[i]] = g;
                keys.Add(groups[i]);
            }
            groupOfRow[i] = g;
        }
        return (groupOfRow, keys);
    }
}
