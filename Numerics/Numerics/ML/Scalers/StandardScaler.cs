using CSharpNumerics.ML.Scalers.Interfaces;
using Numerics.Objects;
using System;

namespace CSharpNumerics.ML.Scalers
{
    public class StandardScaler : IScaler
    {
        public double[] Means;
        public double[] StdDevs;

        public void Fit(double[,] data)
        {
            int rows = data.GetLength(0);
            int cols = data.GetLength(1);

            Means = new double[cols];
            StdDevs = new double[cols];

            for (int j = 0; j < cols; j++)
            {
                double sum = 0.0;
                for (int i = 0; i < rows; i++)
                {
                    sum += data[i, j];
                }
                Means[j] = sum / rows;

                double sumSquares = 0.0;
                for (int i = 0; i < rows; i++)
                {
                    sumSquares += Math.Pow(data[i, j] - Means[j], 2);
                }
                StdDevs[j] = Math.Sqrt(sumSquares / rows);
            }
        }

        public Matrix Transform(Matrix X)
        {
            int rows = X.rowLength;
            int cols = X.columnLength;
            double[,] transformedData = new double[rows, cols];

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    transformedData[i, j] = (X.values[i, j] - Means[j]) / (StdDevs[j] + 1e-12);
                }
            }

            return new Matrix(transformedData);
        }

        public Matrix FitTransform(Matrix X)
        {
            Fit(X.values);
            return Transform(X);
        }

        public IScaler Clone()
        {
            return new StandardScaler
            {
                Means = Means == null ? null : (double[])Means.Clone(),
                StdDevs = StdDevs == null ? null : (double[])StdDevs.Clone()
            };
        }

    }


}
