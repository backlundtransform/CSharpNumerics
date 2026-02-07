using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Objects;
using Numerics.Objects;
using System;
using System.Collections.Generic;


namespace CSharpNumerics.ML.Models.Classification;


public class Logistic : IClassificationModel, IHasHyperparameters
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

    public void Fit(Matrix X, VectorN y)
    {
        if (X.rowLength != y.Length)
            throw new ArgumentException("X.Rows must match y.Length.");

        Matrix Xdesign = FitIntercept ? AddInterceptColumn(X) : X;

        int nSamples = Xdesign.rowLength;
        int nFeatures = Xdesign.columnLength;
        _weights = new VectorN(nFeatures);

        for (int iter = 0; iter < MaxIterations; iter++)
        {
            VectorN predictions = Sigmoid(Xdesign * _weights);
            VectorN errors = predictions - y;

            VectorN gradient = (Xdesign.Transpose() * errors) / nSamples;

            if (RegularizationStrength > 0.0)
                gradient += (RegularizationStrength / nSamples) * _weights;

            _weights -= LearningRate * gradient;

            if (gradient.Norm() < Tolerance)
                break;
        }

        _fitted = true;
    }

    public VectorN Predict(Matrix X)
    {
        if (!_fitted)
            throw new InvalidOperationException("Model has not been fitted.");

        Matrix Xdesign = FitIntercept ? AddInterceptColumn(X) : X;
        VectorN probabilities = Sigmoid(Xdesign * _weights);

        var labels = new double[probabilities.Length];
        for (int i = 0; i < probabilities.Length; i++)
            labels[i] = probabilities[i] >= 0.5 ? 1.0 : 0.0;

        return new VectorN(labels);
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
