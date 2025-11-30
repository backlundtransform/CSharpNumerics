using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CSharpNumerics.Regression
{
    namespace CSharpNumerics.Regression
    {
        internal class Linear
        {
            // Model parameters
            private double[] _weights;  // w0 ... wn
            private double _bias;       // b

            // Hyperparameters
            public double LearningRate { get; set; } = 0.01;
            public int Epochs { get; set; } = 1000;
            public double L2Regularization { get; set; } = 0.0;  // lambda

            public Linear(int featureCount)
            {
                _weights = new double[featureCount];
                _bias = 0.0;
            }

            public void Fit(double[][] X, double[] y)
            {
                int n = X.Length;            // number of samples
                int d = _weights.Length;     // number of features

                for (int epoch = 0; epoch < Epochs; epoch++)
                {
                    double[] gradW = new double[d];
                    double gradB = 0.0;

                    // Compute gradients
                    for (int i = 0; i < n; i++)
                    {
                        double pred = PredictSingle(X[i]);
                        double error = pred - y[i];

                        // accumulate gradients
                        for (int j = 0; j < d; j++)
                            gradW[j] += error * X[i][j];

                        gradB += error;
                    }

                    // Apply L2 regularization on weights
                    for (int j = 0; j < d; j++)
                        gradW[j] = gradW[j] / n + L2Regularization * _weights[j];

                    gradB /= n;

                    // Update rule
                    for (int j = 0; j < d; j++)
                        _weights[j] -= LearningRate * gradW[j];

                    _bias -= LearningRate * gradB;
                }
            }

            public double Predict(double[] x)
            {
                return PredictSingle(x);
            }

            public double[] Predict(double[][] X)
            {
                var preds = new double[X.Length];
                for (int i = 0; i < X.Length; i++)
                    preds[i] = PredictSingle(X[i]);

                return preds;
            }

            private double PredictSingle(double[] x)
            {
                double sum = _bias;
                for (int j = 0; j < _weights.Length; j++)
                    sum += _weights[j] * x[j];

                return sum;
            }

            public double MeanSquaredError(double[][] X, double[] y)
            {
                double loss = 0;
                int n = X.Length;
                for (int i = 0; i < n; i++)
                {
                    double e = PredictSingle(X[i]) - y[i];
                    loss += e * e;
                }
                return loss / n;
            }
        }
    }
