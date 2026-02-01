using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Objects;
using Numerics.Objects;

namespace CSharpNumerics.ML.Models.Regression;

using System;
using System.Collections.Generic;

public class Linear : IRegressionModel, IHasHyperparameters
{
    private VectorN _weights;  
    private bool _fitted = false;

    public bool FitIntercept { get; private set; } = true;
    public double LearningRate { get; private set; } = 0.01;
 

    public void SetHyperParameters(Dictionary<string, object> parameters)
    {
        if (parameters == null) return;

        if (parameters.TryGetValue("LearningRate", out object value))
            LearningRate = Convert.ToDouble(value);
        if (parameters.TryGetValue("FitIntercept", out object value1))
            FitIntercept = Convert.ToBoolean(value1);
    }

    public void Fit(Matrix X, VectorN y)
    {
        if (X.rowLength != y.Length)
            throw new ArgumentException("X.Rows must match y.Length.");

    
        Matrix Xdesign = FitIntercept ? AddInterceptColumn(X) : X;

  
        Matrix Xt = Xdesign.Transpose();
        Matrix XtX = Xt * Xdesign;
        Matrix XtX_inv = XtX.Inverse();   
        VectorN XtY = Xt * y;            

        _weights = XtX_inv * XtY;       
        _fitted = true;
    }

    public VectorN Predict(Matrix X)
    {
        if (!_fitted)
            throw new InvalidOperationException("Model has not been fitted.");

        Matrix Xdesign = FitIntercept ? AddInterceptColumn(X) : X;

        return Xdesign * _weights; 
    }

    private Matrix AddInterceptColumn(Matrix X)
    {
        Matrix M = new Matrix(X.rowLength, X.columnLength + 1);

        for (int i = 0; i < X.rowLength; i++)
        {
            M.values[i, 0] = 1.0; // bias
            for (int j = 0; j < X.columnLength; j++)
                M.values[i, j + 1] = X.values[i, j];
        }
        return M;
    }

    public IModel Clone()
    {
        var clone = new Linear();
        clone.SetHyperParameters(new Dictionary<string, object>
        {
            ["LearningRate"] = LearningRate,
            ["FitIntercept"] = FitIntercept
        });
        return clone;
    }
}
