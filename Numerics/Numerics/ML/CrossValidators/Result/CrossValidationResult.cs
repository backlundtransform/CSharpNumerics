
using CSharpNumerics.Numerics.Objects;
using System.Collections.Generic;


namespace CSharpNumerics.ML.CrossValidators.Result;

public class CrossValidationResult {
    public Dictionary<Pipeline, double> Scores { get; } = new();
    public Pipeline BestPipeline { get; set; }
    public double BestScore { get; set; } 
    public double CoefficientOfDetermination { get; set; } = 0; 
    public Matrix ConfusionMatrix { get; set; } = new Matrix(); 
    public VectorN ActualValues { get; set; } = new VectorN(); 
    public VectorN PredictedValues { get; set; } = new VectorN(); 

    /// <summary>
    /// Out-of-fold class probabilities of the best pipeline, row by row in the same order as
    /// <see cref="ActualValues"/>; null when the validator or model does not provide them.
    /// </summary>
    public Matrix? PredictedProbabilities { get; set; }

    /// <summary>
    /// Index in the input data of each row of <see cref="ActualValues"/>; null when the validator
    /// does not record it.
    /// </summary>
    public int[] RowIndices { get; set; }
}
