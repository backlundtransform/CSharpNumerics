namespace CSharpNumerics.ML.Enums;

/// <summary>
/// How a cross-validator scores a classification pipeline. Higher is always better, so losses are
/// negated. All but <see cref="Accuracy"/> need a binary <c>IProbabilisticClassifier</c> and use
/// the probability of class 1.
/// </summary>
public enum ClassificationScoring
{
    /// <summary>Share of rows classified correctly.</summary>
    Accuracy,

    /// <summary>Area under the ROC curve.</summary>
    RocAuc,

    /// <summary>Average precision (area under the precision–recall curve).</summary>
    AveragePrecision,

    /// <summary>Minus the log loss (cross-entropy).</summary>
    NegativeLogLoss,

    /// <summary>Minus the Brier score.</summary>
    NegativeBrierScore,
}
