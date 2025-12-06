
using Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Text;

namespace CSharpNumerics.ML
{
    public class CrossValidator(List<Pipeline> pipelines, Matrix Features, double[] Labels)
    {
        public List<Pipeline> Pipelines { get; } = pipelines;
        public Matrix X { get; } = Features;
        public double[] Y { get; } = Labels;

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
                int end = start + foldSize;

                var (Xtrain, Ytrain, Xval, Yval) = SplitFold(X, Y, start, end);

                pipe.Fit(Xtrain, Ytrain);

                var pred = pipe.Predict(Xval);

                totalScore += Accuracy(pred, Yval);
            }

            return totalScore / folds;
        }

        private double Accuracy(double[] pred, double[] actual)
        {
            int correct = 0;

            for (int i = 0; i < pred.Length; i++)
                if (Math.Round(pred[i]) == actual[i])
                    correct++;

            return (double)correct / pred.Length;
        }

        private (Matrix, double[], Matrix, double[]) SplitFold(Matrix X, double[] Y, int start, int end)
        {
            // Här kan du själv implementera slicing av Matrix
       
            throw new NotImplementedException();
        }
    }
}
