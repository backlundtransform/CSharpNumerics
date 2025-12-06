using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Objects;
using Numerics.Objects;

namespace CSharpNumerics.ML.Models.Regression
{
    using System;

    public class Linear : IRegressionModel
    {
        private VectorN _weights;   // β-koefficienterna
        private bool _fitted = false;

        public bool FitIntercept { get; }

        public Linear(bool fitIntercept = true)
        {
            FitIntercept = fitIntercept;
        }

        public void Fit(Matrix X, VectorN y)
        {
            if (X.rowLength != y.Length)
                throw new ArgumentException("X.Rows must match y.Length.");

            // 1. Lägg till en kolumn med 1:or om vi vill ha intercept
            Matrix Xdesign = FitIntercept ? AddInterceptColumn(X) : X;

            // 2. Normalekvation: β = (XᵀX)^(-1) Xᵀ y
            Matrix Xt = Xdesign.Transpose();
            Matrix XtX = Xt * Xdesign;
            Matrix XtX_inv = XtX.Inverse();   // Krav: du har en Inverse() i Matrix
            VectorN XtY = Xt * y;             // Xt (MxN) * y (N) = M

            _weights = XtX_inv * XtY;         // β
            _fitted = true;
        }

        public VectorN Predict(Matrix X)
        {
            if (!_fitted)
                throw new InvalidOperationException("Model has not been fitted.");

            Matrix Xdesign = FitIntercept ? AddInterceptColumn(X) : X;

            return Xdesign * _weights; // ger VectorN
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
    }
}
 