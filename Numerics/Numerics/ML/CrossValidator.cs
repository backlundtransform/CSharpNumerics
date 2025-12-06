using CSharpNumerics.Objects;
using Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML
{
    public class CrossValidator
    {
        public List<Pipeline> Pipelines { get; }
        public Matrix X { get; }
        public VectorN Y { get; }

        public CrossValidator(List<Pipeline> pipelines, Matrix features, VectorN labels)
        {
            Pipelines = pipelines;
            X = features;
            Y = labels;
        }

        public Dictionary<Pipeline, double> Run(int folds = 5)
        {
            var results = new Dictionary<Pipeline, double>();

            foreach (var pipe in Pipelines)
            {
                double score = EvaluatePipeline(pipe, folds);
                results[pipe] = score;
            }

            return results;
        }

        private double EvaluatePipeline(Pipeline pipe, int folds)
        {
            int n = Y.Length;
            int foldSize = n / folds;
            double totalScore = 0;

            for (int f = 0; f < folds; f++)
            {
                int start = f * foldSize;
                int end = (f == folds - 1) ? n : start + foldSize; 

                var (Xtrain, Ytrain, Xval, Yval) = SplitFold(X, Y, start, end);

                pipe.Fit(Xtrain, Ytrain);

                var pred = pipe.Predict(Xval);

                totalScore += Accuracy(pred, Yval);
            }

            return totalScore / folds;
        }

        private double Accuracy(VectorN pred, VectorN actual)
        {
            int correct = 0;
            for (int i = 0; i < pred.Length; i++)
                if (Math.Round(pred[i]) == actual[i])
                    correct++;
            return (double)correct / pred.Length;
        }

        private (Matrix Xtrain, VectorN Ytrain, Matrix Xval, VectorN Yval) SplitFold(Matrix X, VectorN Y, int start, int end)
        {
            int n = X.rowLength;
            int valSize = end - start;
            int trainSize = n - valSize;

            var XtrainValues = new double[trainSize, X.columnLength];
            var XvalValues = new double[valSize, X.columnLength];

            var YtrainValues = new double[trainSize];
            var YvalValues = new double[valSize];

            int trainIdx = 0;
            int valIdx = 0;

            for (int i = 0; i < n; i++)
            {
                if (i >= start && i < end)
                {
                    for (int j = 0; j < X.columnLength; j++)
                        XvalValues[valIdx, j] = X.values[i, j];
                    YvalValues[valIdx] = Y[i];
                    valIdx++;
                }
                else
                {
                    for (int j = 0; j < X.columnLength; j++)
                        XtrainValues[trainIdx, j] = X.values[i, j];
                    YtrainValues[trainIdx] = Y[i];
                    trainIdx++;
                }
            }

            return (new Matrix(XtrainValues), new VectorN(YtrainValues),
                    new Matrix(XvalValues), new VectorN(YvalValues));
        }
    }
}
