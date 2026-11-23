using CSharpNumerics.ML.CrossValidators;
using CSharpNumerics.ML.CrossValidators.Interfaces;
using CSharpNumerics.ML.Enums;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.Experiment;

/// <summary>
/// Describes a cross-validation strategy without binding to specific pipelines.
/// Used by <see cref="SupervisedExperimentBuilder"/> to create cross-validators
/// once the pipeline list has been built from the grid.
///
/// <example>
/// <code>
/// CrossValidatorConfig.KFold(folds: 5)
/// CrossValidatorConfig.StratifiedKFold(folds: 10)
/// CrossValidatorConfig.ShuffleSplit(nSplits: 20, testSize: 0.25)
/// CrossValidatorConfig.LeaveOneOut()
/// CrossValidatorConfig.StratifiedGroupKFold(starIds, folds: 5, scoring: ClassificationScoring.RocAuc)
/// CrossValidatorConfig.MonteCarlo(iterations: 200, testSize: 0.2, seed: 42)
/// </code>
/// </example>
/// </summary>
public class CrossValidatorConfig
{
    /// <summary>Human-readable name for this CV strategy (used as key in result dictionaries).</summary>
    public string Name { get; }

    private readonly Func<List<Pipeline>, ICrossValidator> _factory;

    private CrossValidatorConfig(string name, Func<List<Pipeline>, ICrossValidator> factory)
    {
        Name = name;
        _factory = factory;
    }

    /// <summary>Create an <see cref="ICrossValidator"/> bound to the given pipeline list.</summary>
    internal ICrossValidator Create(List<Pipeline> pipelines) => _factory(pipelines);

    // ═══════════════════════════════════════════════════════════════
    //  Static factory methods
    // ═══════════════════════════════════════════════════════════════

    /// <summary>Standard K-fold cross-validation.</summary>
    public static CrossValidatorConfig KFold(int folds = 5)
        => new("KFold", p => new KFoldCrossValidator(p, folds));

    /// <summary>Stratified K-fold (preserves class distribution per fold). Classification only.</summary>
    public static CrossValidatorConfig StratifiedKFold(int folds = 5)
        => new("StratifiedKFold", p => new StratifiedKFoldCrossValidator(p, folds));

    /// <summary>
    /// Group K-fold: all rows with the same group id (one per row of the experiment's data) stay
    /// in the same fold. Optional sample weights are used for training and scoring.
    /// </summary>
    public static CrossValidatorConfig GroupKFold(IReadOnlyList<long> groups, int folds = 5,
        VectorN? sampleWeights = null, ClassificationScoring scoring = ClassificationScoring.Accuracy)
        => new("GroupKFold", p => new GroupKFoldCrossValidator(p, groups, folds, sampleWeights, scoring));

    /// <summary>
    /// Stratified group K-fold: groups stay together and each fold keeps about the overall class
    /// proportions. Classification only.
    /// </summary>
    public static CrossValidatorConfig StratifiedGroupKFold(IReadOnlyList<long> groups, int folds = 5,
        VectorN? sampleWeights = null, ClassificationScoring scoring = ClassificationScoring.Accuracy)
        => new("StratifiedGroupKFold", p => new StratifiedGroupKFoldCrossValidator(p, groups, folds, sampleWeights, scoring));

    /// <summary>Random shuffle-split cross-validation.</summary>
    public static CrossValidatorConfig ShuffleSplit(
        int nSplits = 5,
        double testSize = 0.2,
        double trainSize = 0.8,
        int? randomState = null)
        => new("ShuffleSplit", p => new ShuffleSplitCrossValidator(p, nSplits, testSize, trainSize, randomState));

    /// <summary>Leave-one-out cross-validation.</summary>
    public static CrossValidatorConfig LeaveOneOut()
        => new("LeaveOneOut", p => new LeaveOneOutCrossValidator(p));

    /// <summary>Monte Carlo cross-validation with full score distributions.</summary>
    public static CrossValidatorConfig MonteCarlo(
        int iterations = 100,
        double testSize = 0.2,
        double trainSize = 0.8,
        int? seed = null,
        double confidenceLevel = 0.95)
        => new("MonteCarlo", p => new MonteCarloCrossValidator(
            p, iterations, testSize, trainSize, seed, confidenceLevel));

    /// <summary>
    /// Create a custom cross-validator configuration from any factory function.
    /// </summary>
    /// <param name="name">Name for the CV strategy (used as dictionary key).</param>
    /// <param name="factory">Factory that receives the expanded pipeline list and returns a ready-to-run CV.</param>
    public static CrossValidatorConfig Custom(string name, Func<List<Pipeline>, ICrossValidator> factory)
        => new(name, factory);
}
