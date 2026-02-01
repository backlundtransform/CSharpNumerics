using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.ML.Selector.Interfaces;
using CSharpNumerics.Objects;
using Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Text;

namespace CSharpNumerics.ML.Selector
{
    public class SelectKBest: ISelector, IHasHyperparameters
    {
        public int K { get; set; }
        private int[] selectedIndices;

        public void SetHyperParameters(Dictionary<string, object> parameters)
        {
            if (parameters == null) return;

            if (parameters.ContainsKey("K"))
                K = Convert.ToInt32(parameters["K"]);
        }

        /// <summary>
        /// Fit selector on X and y. Computes a score per column and selects top K.
        /// </summary>
        public void Fit(Matrix X, VectorN y)
        {
            if (X.columnLength < K)
                throw new ArgumentException("K cannot be greater than number of features.");

            int nFeatures = X.columnLength;
            double[] scores = new double[nFeatures];

         
            for (int j = 0; j < nFeatures; j++)
            {
                VectorN column = X.ColumnSlice(j);
                scores[j] = Math.Abs(Correlation(column, y));
            }

            selectedIndices = ArgSortDescending(scores, K);
        }

        public Matrix Transform(Matrix X)
        {
            if (selectedIndices == null)
                throw new InvalidOperationException("SelectKBest must be fitted before calling Transform.");

            var newX = new double[X.rowLength, K];

            for (int i = 0; i < X.rowLength; i++)
            {
                for (int k = 0; k < K; k++)
                {
                    int colIndex = selectedIndices[k];
                    newX[i, k] = X.values[i, colIndex];
                }
            }

            return new Matrix(newX);
        }

        public Matrix FitTransform(Matrix X, VectorN y)
        {
            Fit(X, y);
            return Transform(X);
        }

        public ISelector Clone()
        {
            return new SelectKBest { K = K };
        }

      
        private double Correlation(VectorN a, VectorN b)
        {
            if (a.Length != b.Length)
                throw new ArgumentException("Vectors must have same length.");

            int n = a.Length;
            double meanA = 0, meanB = 0;

            for (int i = 0; i < n; i++)
            {
                meanA += a[i];
                meanB += b[i];
            }

            meanA /= n;
            meanB /= n;

            double num = 0;
            double denomA = 0;
            double denomB = 0;

            for (int i = 0; i < n; i++)
            {
                double da = a[i] - meanA;
                double db = b[i] - meanB;
                num += da * db;
                denomA += da * da;
                denomB += db * db;
            }

            double denom = Math.Sqrt(denomA * denomB);

            if (denom == 0) return 0;

            return num / denom;
        }

      
        private int[] ArgSortDescending(double[] scores, int topK)
        {
            var indices = new List<int>();
            for (int i = 0; i < scores.Length; i++)
                indices.Add(i);

            indices.Sort((i, j) => scores[j].CompareTo(scores[i]));

            return indices.GetRange(0, topK).ToArray();
        }
    }
}
