using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.SingleObjective;
using CSharpNumerics.Numerics.Optimization.Strategies;
using System;
using System.Collections.Generic;


namespace CSharpNumerics.ML.Models.Classification;


public class Logistic : IProbabilisticClassifier, ISampleWeightedModel, IHasHyperparameters
{
    private VectorN _weights;
    private bool _fitted = false;
    public bool FitIntercept { get; private set; } = true;
    public int NumClasses { get; private set; } = 2;
    public double LearningRate { get; private set; } = 0.01;
    public int MaxIterations { get; private set; } = 1000;
    public double RegularizationStrength { get; private set; } = 0.0;
    public double Tolerance { get; private set; } = 1e-7;

    public void SetHyperParameters(Dictionary<string, object> parameters)
    {
        if (parameters == null) return;

        if (parameters.ContainsKey("LearningRate"))
            LearningRate = Convert.ToDouble(parameters["LearningRate"]);

        if (parameters.ContainsKey("MaxIterations"))
            MaxIterations = Convert.ToInt32(parameters["MaxIterations"]);

        if (parameters.ContainsKey("FitIntercept"))
            FitIntercept = Convert.ToBoolean(parameters["FitIntercept"]);

        if (parameters.ContainsKey("RegularizationStrength"))
            RegularizationStrength = Convert.ToDouble(parameters["RegularizationStrength"]);

        if (parameters.ContainsKey("Tolerance"))
            Tolerance = Convert.ToDouble(parameters["Tolerance"]);
    }

    public void Fit(Matrix X, VectorN y) => FitWeighted(X, y, null);

    /// <summary>
    /// Trains with a weight per row: the gradient is the weighted mean of the per-row gradients.
    /// </summary>
    public void Fit(Matrix X, VectorN y, VectorN sampleWeights) => FitWeighted(X, y, sampleWeights);

    private void FitWeighted(Matrix X, VectorN y, VectorN? sampleWeights)
    {
        if (X.rowLength != y.Length)
            throw new ArgumentException("X.Rows must match y.Length.");

        double[] weights = ClassWeights.Validate(sampleWeights, y.Length);
        double weightTotal = 0;
        foreach (double w in weights) weightTotal += w;

        Matrix Xdesign = FitIntercept ? AddInterceptColumn(X) : X;

        int nSamples = Xdesign.rowLength;
        int nFeatures = Xdesign.columnLength;
        _weights = new VectorN(nFeatures);

        var optimizer = new GradientDescent(LearningRate, l2: RegularizationStrength / nSamples);
        var convergence = new MaxIterationsOrTolerance(MaxIterations, Tolerance);

        for (int iter = 0; ; iter++)
        {
            VectorN predictions = Sigmoid(Xdesign * _weights);
            VectorN errors = predictions - y;
            for (int i = 0; i < errors.Length; i++)
                errors[i] *= weights[i];

            VectorN gradient = (Xdesign.Transpose() * errors) / weightTotal;

            if (convergence.HasConverged(iter, 0, gradient.Norm()))
                break;

            _weights = optimizer.Step(_weights, gradient);
        }

        _fitted = true;
    }

    public VectorN Predict(Matrix X)
    {
        VectorN probabilities = PositiveProbability(X);

        var labels = new double[probabilities.Length];
        for (int i = 0; i < probabilities.Length; i++)
            labels[i] = probabilities[i] >= 0.5 ? 1.0 : 0.0;

        return new VectorN(labels);
    }

    /// <summary>Column 0: probability of class 0; column 1: probability of class 1.</summary>
    public Matrix PredictProbabilities(Matrix X)
    {
        VectorN p = PositiveProbability(X);
        var result = new Matrix(p.Length, 2);
        for (int i = 0; i < p.Length; i++)
        {
            result.values[i, 0] = 1.0 - p[i];
            result.values[i, 1] = p[i];
        }
        return result;
    }

    private VectorN PositiveProbability(Matrix X)
    {
        if (!_fitted)
            throw new InvalidOperationException("Model has not been fitted.");

        Matrix Xdesign = FitIntercept ? AddInterceptColumn(X) : X;
        return Sigmoid(Xdesign * _weights);
    }

    private VectorN Sigmoid(VectorN z)
    {
        VectorN result = new VectorN(z.Length);
        for (int i = 0; i < z.Length; i++)
        {
            if (z[i] >= 0)
            {
                double ex = Math.Exp(-z[i]);
                result[i] = 1.0 / (1.0 + ex);
            }
            else
            {
                double ex = Math.Exp(z[i]);
                result[i] = ex / (1.0 + ex);
            }
        }
        return result;
    }

    private Matrix AddInterceptColumn(Matrix X)
    {
        Matrix M = new Matrix(X.rowLength, X.columnLength + 1);

        for (int i = 0; i < X.rowLength; i++)
        {
            M.values[i, 0] = 1.0;
            for (int j = 0; j < X.columnLength; j++)
                M.values[i, j + 1] = X.values[i, j];
        }
        return M;
    }

    public IModel Clone()
    {
        var clone = new Logistic();
        clone.SetHyperParameters(new Dictionary<string, object>
        {
            ["LearningRate"] = LearningRate,
            ["MaxIterations"] = MaxIterations,
            ["FitIntercept"] = FitIntercept,
            ["RegularizationStrength"] = RegularizationStrength,
            ["Tolerance"] = Tolerance
        });
        return clone;
    }
}
