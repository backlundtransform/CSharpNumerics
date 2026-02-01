using CSharpNumerics.ML.Scalers.Interfaces;
using Numerics.Objects;


namespace CSharpNumerics.ML.Scalers;

public class MinMaxScaler : IScaler
{
    private double[] min;
    private double[] max;

    public Matrix FitTransform(Matrix X)
    {
        int rows = X.rowLength;
        int cols = X.columnLength;

        min = new double[cols];
        max = new double[cols];

        double[,] result = new double[rows, cols];

        for (int j = 0; j < cols; j++)
        {
            double mn = double.PositiveInfinity;
            double mx = double.NegativeInfinity;

            for (int i = 0; i < rows; i++)
            {
                double v = X.values[i, j];
                if (v < mn) mn = v;
                if (v > mx) mx = v;
            }

            min[j] = mn;
            max[j] = mx;
        }

        return Transform(X);
    }

    public Matrix Transform(Matrix X)
    {
        int rows = X.rowLength;
        int cols = X.columnLength;

        double[,] result = new double[rows, cols];

        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                result[i, j] = (X.values[i, j] - min[j]) /
                               (max[j] - min[j] + 1e-9);

        return new Matrix(result);
    }

    public IScaler Clone()
    {
        var clone = new MinMaxScaler();
        if (min != null) clone.min = (double[])min.Clone();
        if (max != null) clone.max = (double[])max.Clone();
        return clone;
    }
}
