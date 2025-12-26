using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Objects;
using Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace CSharpNumerics.ML.Models.Classification
{
    public class NaiveBayes : IClassificationModel
    {
        private Dictionary<int, VectorN> _means;
        private Dictionary<int, VectorN> _variances;
        private Dictionary<int, double> _priors;
        private int _numFeatures;
        private bool _fitted;
       
        public int NumClasses => throw new NotImplementedException();

        public void Fit(Matrix X, VectorN y)
        {
            _numFeatures = X.columnLength;
            _means = new();
            _variances = new();
            _priors = new();

            var classes = y.Values.Distinct().Select(v => (int)v).ToArray();

            foreach (var cls in classes)
            {
                var idx = Enumerable.Range(0, y.Length)
                                    .Where(i => (int)y[i] == cls)
                                    .ToArray();

                _priors[cls] = (double)idx.Length / y.Length;

                var mean = new VectorN(_numFeatures);
                var var = new VectorN(_numFeatures);

                for (int j = 0; j < _numFeatures; j++)
                {
                    var vals = idx.Select(i => X.values[i, j]).ToArray();
                    mean[j] = vals.Average();
                    var[j] = vals.Select(v => Math.Pow(v - mean[j], 2)).Average() + 1e-9;
                }

                _means[cls] = mean;
                _variances[cls] = var;
            }

            _fitted = true;
        }

        public VectorN Predict(Matrix X)
        {
            if (!_fitted)
                throw new InvalidOperationException("Model not fitted");

            VectorN preds = new VectorN(X.rowLength);

            for (int i = 0; i < X.rowLength; i++)
            {
                double bestScore = double.NegativeInfinity;
                int bestClass = -1;

                foreach (var cls in _priors.Keys)
                {
                    double logProb = Math.Log(_priors[cls]);

                    for (int j = 0; j < _numFeatures; j++)
                    {
                        double mean = _means[cls][j];
                        double var = _variances[cls][j];
                        double x = X.values[i, j];

                        logProb += -0.5 * Math.Log(2 * Math.PI * var)
                                   - Math.Pow(x - mean, 2) / (2 * var);
                    }

                    if (logProb > bestScore)
                    {
                        bestScore = logProb;
                        bestClass = cls;
                    }
                }

                preds[i] = bestClass;
            }

            return preds;
        }
    }

}
