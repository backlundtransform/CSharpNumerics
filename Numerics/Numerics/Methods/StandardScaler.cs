using System;

namespace CSharpNumerics.Methods
{
    public class StandardScaler
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

        public double[,] Transform(double[,] data)
        {
            int rows = data.GetLength(0);
            int cols = data.GetLength(1);
            double[,] transformedData = new double[rows, cols];

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    transformedData[i, j] = (data[i, j] - Means[j]) / StdDevs[j];
                }
            }

            return transformedData;
        }

        public double[,] FitTransform(double[,] data)
        {
            Fit(data);
            return Transform(data);
        }
    }
}
