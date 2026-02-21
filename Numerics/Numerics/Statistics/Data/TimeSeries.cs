using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace CSharpNumerics.Statistics.Data;

public class TimeSeries
{
    public DateTime[] Time { get; }
    public double[][] Data { get; }
    public string[] Cols { get; }

    public int RowCount => Time.Length;
    public int ColumnCount => Data.Length;

    public TimeSeries(DateTime[] time, double[][] data, string[] cols)
    {
        if (data.Length != cols.Length)
            throw new ArgumentException("Data columns must match Cols");

        if (data.Any(c => c.Length != time.Length))
            throw new ArgumentException("All columns must match Time length");

        Time = time;
        Data = data;
        Cols = cols;
    }


    public Matrix ToMatrix(int? excludeCol = null)
    {
        int rows = Time.Length;
        int cols = Data.Length;

        var values = new double[rows, cols];

        int cNew = 0;
        for (int c = 0; c < Data.Length; c++)
        {
            if (c == excludeCol) continue;

            var column = Data[c];
            for (int r = 0; r < rows; r++)
                values[r, cNew] = column[r];

            cNew++;
        }


        return new Matrix(values);
    }
    public static TimeSeries FromTimeSerie(
    List<TimeSerie> serie,
    string columnName = "Value")
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
            [values],
            [columnName]
        );
    }

    public static TimeSeries FromCsv(string path, char sep = ',')
    {
        var lines = File.ReadAllLines(path);

        var header = lines[0].Split(sep);
        if (header.Length < 2)
            throw new InvalidOperationException("CSV must have at least one feature column");

        var cols = header.Skip(1).ToArray();
        int featureCount = cols.Length;
        int rowCount = lines.Length - 1;

        var time = new DateTime[rowCount];
        var data = new double[featureCount][];

        for (int c = 0; c < featureCount; c++)
            data[c] = new double[rowCount];

        for (int r = 0; r < rowCount; r++)
        {
            var parts = lines[r + 1].Split(sep);

            time[r] = DateTime.Parse(parts[0]);

            for (int c = 0; c < featureCount; c++)
            {
                data[c][r] = double.Parse(
                    parts[c + 1],
                    CultureInfo.InvariantCulture
                );
            }
        }

        return new TimeSeries(time, data, cols);
    }
}
