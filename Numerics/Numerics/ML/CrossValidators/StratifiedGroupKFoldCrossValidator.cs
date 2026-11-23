using CSharpNumerics.ML.CrossValidators.Interfaces;
using CSharpNumerics.ML.CrossValidators.Result;
using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.CrossValidators;

/// <summary>
/// Group k-fold for classification that also keeps the class proportions about the same in every
/// fold: all rows of a group stay together, and groups are placed so that each fold gets its
/// share of every class.
/// </summary>
/// <remarks>
/// Greedy, as scikit-learn's StratifiedGroupKFold without shuffling: groups with the most uneven
/// class mix are placed first, each in the fold where it makes the per-class shares across folds
/// most even (smallest mean standard deviation); ties go to the smaller fold. Deterministic.
/// </remarks>
public class StratifiedGroupKFoldCrossValidator : ICrossValidator
{
    public List<Pipeline> Pipelines { get; }
    public int Folds { get; }

    /// <summary>Group id of each row of the data passed to <see cref="Run"/>.</summary>
    public IReadOnlyList<long> Groups { get; }

    /// <summary>Optional weight per row, used for training (the models must support it) and scoring.</summary>
    public VectorN? SampleWeights { get; }

    public ClassificationScoring Scoring { get; }

    public StratifiedGroupKFoldCrossValidator(List<Pipeline> pipelines, IReadOnlyList<long> groups, int folds = 5,
        VectorN? sampleWeights = null, ClassificationScoring scoring = ClassificationScoring.Accuracy)
    {
        if (folds < 2) throw new ArgumentOutOfRangeException(nameof(folds), "At least 2 folds are required.");
        Pipelines = pipelines;
        Groups = groups ?? throw new ArgumentNullException(nameof(groups));
        Folds = folds;
        SampleWeights = sampleWeights;
        Scoring = scoring;
    }

    public StratifiedGroupKFoldCrossValidator(PipelineGrid pipelineGrid, IReadOnlyList<long> groups, int folds = 5,
        VectorN? sampleWeights = null, ClassificationScoring scoring = ClassificationScoring.Accuracy)
        : this([.. pipelineGrid.Expand()], groups, folds, sampleWeights, scoring)
    {
    }

    public CrossValidationResult Run(Matrix X, VectorN y)
    {
        foreach (var pipe in Pipelines)
            if (pipe.Model is not IClassificationModel)
                throw new InvalidOperationException(
                    $"{pipe.Model.GetType().Name} is not a classification model; StratifiedGroupKFoldCrossValidator is for classification.");

        return FoldRunner.Run(Pipelines, X, y, FoldOfEachRow(y), Folds, SampleWeights, Scoring);
    }

    /// <summary>The fold (0 … Folds − 1) each row is tested in, for class labels <paramref name="y"/>.</summary>
    public int[] FoldOfEachRow(VectorN y)
    {
        var (groupOfRow, keys) = FoldRunner.IndexGroups(Groups, y.Length);
        if (keys.Count < Folds)
            throw new InvalidOperationException($"{keys.Count} groups cannot fill {Folds} folds; use fewer folds.");

        var classes = y.Values.Distinct().OrderBy(v => v).ToList();
        var classIndex = classes.Select((c, i) => (c, i)).ToDictionary(p => p.c, p => p.i);
        int numClasses = classes.Count;

        var groupCounts = new double[keys.Count, numClasses];
        var classTotals = new double[numClasses];
        var groupSize = new int[keys.Count];
        for (int i = 0; i < y.Length; i++)
        {
            int c = classIndex[y[i]];
            groupCounts[groupOfRow[i], c]++;
            classTotals[c]++;
            groupSize[groupOfRow[i]]++;
        }

        // Most uneven class mix first (standard deviation of the group's class counts).
        double Spread(int g)
        {
            double mean = 0;
            for (int c = 0; c < numClasses; c++) mean += groupCounts[g, c];
            mean /= numClasses;
            double sum = 0;
            for (int c = 0; c < numClasses; c++) sum += (groupCounts[g, c] - mean) * (groupCounts[g, c] - mean);
            return Math.Sqrt(sum / numClasses);
        }

        var foldCounts = new double[Folds, numClasses];
        var foldSize = new int[Folds];
        var foldOfGroup = new int[keys.Count];
        foreach (int g in Enumerable.Range(0, keys.Count).OrderByDescending(Spread).ThenBy(g => keys[g]))
        {
            int best = -1;
            double bestEvenness = double.MaxValue;
            for (int f = 0; f < Folds; f++)
            {
                for (int c = 0; c < numClasses; c++) foldCounts[f, c] += groupCounts[g, c];
                double evenness = MeanStdOfShares(foldCounts, classTotals);
                for (int c = 0; c < numClasses; c++) foldCounts[f, c] -= groupCounts[g, c];

                bool better = best < 0
                    || evenness < bestEvenness - 1e-12
                    || (Math.Abs(evenness - bestEvenness) <= 1e-12 && foldSize[f] < foldSize[best]);
                if (better)
                {
                    best = f;
                    bestEvenness = evenness;
                }
            }

            foldOfGroup[g] = best;
            foldSize[best] += groupSize[g];
            for (int c = 0; c < numClasses; c++) foldCounts[best, c] += groupCounts[g, c];
        }

        return groupOfRow.Select(g => foldOfGroup[g]).ToArray();
    }

    /// <summary>Mean over classes of the standard deviation, across folds, of each fold's share of that class.</summary>
    private double MeanStdOfShares(double[,] foldCounts, double[] classTotals)
    {
        double total = 0;
        for (int c = 0; c < classTotals.Length; c++)
        {
            double mean = 0;
            for (int f = 0; f < Folds; f++) mean += foldCounts[f, c] / classTotals[c];
            mean /= Folds;
            double sum = 0;
            for (int f = 0; f < Folds; f++)
            {
                double d = foldCounts[f, c] / classTotals[c] - mean;
                sum += d * d;
            }
            total += Math.Sqrt(sum / Folds);
        }
        return total / classTotals.Length;
    }
}
