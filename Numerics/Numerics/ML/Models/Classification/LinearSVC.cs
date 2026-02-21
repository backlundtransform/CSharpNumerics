
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace CSharpNumerics.ML.Models.Classification
{
    public class LinearSVC : IClassificationModel, IHasHyperparameters
    {
        public double C { get; set; } = 1.0;
        public int Epochs { get; set; } = 1000;
        public double LearningRate { get; set; } = 0.01;

        public int NumClasses { get; private set; }

        private VectorN _weights;
        private double _bias;

        public void Fit(Matrix X, VectorN y)
        {
            NumClasses = (int)y.Values.Max() + 1;
            int nSamples = X.rowLength;
            int nFeatures = X.columnLength;

            _weights = new VectorN(new double[nFeatures]);
            _bias = 0;

            for (int epoch = 0; epoch < Epochs; epoch++)
            {
                for (int i = 0; i < nSamples; i++)
                {
                    var xi = X.RowSlice(i);
                    double yi = y[i] == 0 ? -1 : 1;

                    double condition = yi * (_weights.Dot(xi) + _bias);

                    if (condition >= 1)
                    {
                        _weights -= LearningRate * _weights;
                    }
                    else
                    {
                        _weights -= LearningRate * (_weights - C * yi * xi);
                        _bias += LearningRate * C * yi;
                    }
                }
            }
        }

        public VectorN Predict(Matrix X)
        {
            int n = X.rowLength;
            double[] preds = new double[n];

            for (int i = 0; i < n; i++)
            {
                double score = _weights.Dot(X.RowSlice(i)) + _bias;
                preds[i] = score >= 0 ? 1 : 0;
            }

            return new VectorN(preds);
        }


        public void SetHyperParameters(Dictionary<string, object> p)
        {
            if (p.TryGetValue("C", out var c)) C = Convert.ToDouble(c);
            if (p.TryGetValue("LearningRate", out var lr)) LearningRate = Convert.ToDouble(lr);
            if (p.TryGetValue("Epochs", out var it)) Epochs = Convert.ToInt32(it);
        }

        public IModel Clone()
        {
            return new LinearSVC
            {
                C = C,
                Epochs = Epochs,
                LearningRate = LearningRate
            };
        }
    }
}
