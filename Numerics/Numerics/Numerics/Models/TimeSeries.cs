using Numerics.Models;
using Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Text;

namespace CSharpNumerics.Numerics.Models;

public class TimeSeries(DateTime[] time, double[][] data)
{
    public DateTime[] Time { get; } = time;
    public double[][] Data { get; } = data;

    public static TimeSeries FromTimeSerie(List<TimeSerie> serie)
    {
        var time = new DateTime[serie.Count];
        var values = new double[serie.Count];

        for (int i = 0; i < serie.Count; i++)
        {
            time[i] = serie[i].TimeStamp;
            values[i] = serie[i].Value;
        }

        return new TimeSeries(
            time,
            new[] { values }
        );
    }

    public Matrix ToMatrix()
    {
        int rows = Time.Length;
        int cols = Data.Length;

        var values = new double[rows, cols];

        for (int c = 0; c < cols; c++)
        {
            var column = Data[c];
             
            for (int r = 0; r < rows; r++)
            {
                values[r, c] = column[r];
            }
        }

        return new Matrix(values);
    }
}
