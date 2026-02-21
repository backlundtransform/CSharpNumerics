
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
}
