using CSharpNumerics.ML.CrossValidators.Interfaces;
using CSharpNumerics.ML.CrossValidators.Result;
using CSharpNumerics.ML.Enums;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.CrossValidators;

/// <summary>
/// K-fold cross-validation where all rows of a group stay in the same fold, so a model is never
/// tested on a group it saw in training (for example several signals from the same star).
/// </summary>
/// <remarks>
/// Groups are dealt out largest first, each to the fold with the fewest rows so far, which keeps
/// the folds about equally large. The assignment is deterministic. Works for classification and
/// regression; for classes of very different sizes see <see cref="StratifiedGroupKFoldCrossValidator"/>.
/// </remarks>
public class GroupKFoldCrossValidator : ICrossValidator
{
    public List<Pipeline> Pipelines { get; }
    public int Folds { get; }

    /// <summary>Group id of each row of the data passed to <see cref="Run"/>.</summary>
    public IReadOnlyList<long> Groups { get; }

    /// <summary>Optional weight per row, used for training (the models must support it) and scoring.</summary>
    public VectorN? SampleWeights { get; }

    public ClassificationScoring Scoring { get; }

    public GroupKFoldCrossValidator(List<Pipeline> pipelines, IReadOnlyList<long> groups, int folds = 5,
        VectorN? sampleWeights = null, ClassificationScoring scoring = ClassificationScoring.Accuracy)
    {
        if (folds < 2) throw new ArgumentOutOfRangeException(nameof(folds), "At least 2 folds are required.");
        Pipelines = pipelines;
        Groups = groups ?? throw new ArgumentNullException(nameof(groups));
        Folds = folds;
        SampleWeights = sampleWeights;
        Scoring = scoring;
    }

    public GroupKFoldCrossValidator(PipelineGrid pipelineGrid, IReadOnlyList<long> groups, int folds = 5,
        VectorN? sampleWeights = null, ClassificationScoring scoring = ClassificationScoring.Accuracy)
        : this([.. pipelineGrid.Expand()], groups, folds, sampleWeights, scoring)
    {
    }

    public CrossValidationResult Run(Matrix X, VectorN y) =>
        FoldRunner.Run(Pipelines, X, y, FoldOfEachRow(y.Length), Folds, SampleWeights, Scoring);

    /// <summary>The fold (0 … Folds − 1) each of <paramref name="rows"/> rows is tested in.</summary>
    public int[] FoldOfEachRow(int rows)
    {
        var (groupOfRow, keys) = FoldRunner.IndexGroups(Groups, rows);
        if (keys.Count < Folds)
            throw new InvalidOperationException($"{keys.Count} groups cannot fill {Folds} folds; use fewer folds.");

        var sizes = new int[keys.Count];
        foreach (int g in groupOfRow) sizes[g]++;

        var foldOfGroup = new int[keys.Count];
        var foldSize = new int[Folds];
        foreach (int g in Enumerable.Range(0, keys.Count).OrderByDescending(g => sizes[g]).ThenBy(g => keys[g]))
        {
            int fold = 0;
            for (int f = 1; f < Folds; f++)
                if (foldSize[f] < foldSize[fold]) fold = f;
            foldOfGroup[g] = fold;
            foldSize[fold] += sizes[g];
        }

        return groupOfRow.Select(g => foldOfGroup[g]).ToArray();
    }
}
