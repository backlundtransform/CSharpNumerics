using Numerics.Models;
using Numerics.Objects;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
namespace CSharpNumerics.Numerics.Models;



public class Series(int[] index, double[][] data, string[] cols, int[] groups = null)
{

 
    public int[] Index { get; } = index;
    public double[][] Data { get; } = data;
    public int[] Groups { get; } = groups;
    public string[] Cols { get; } = cols;

    public static Series FromSerie(List<Serie> serie)
    {
        var index = new int[serie.Count];
        var values = new double[serie.Count];

        for (int i = 0; i < serie.Count; i++)
        {
            index[i] = (int)serie[i].Index;
            values[i] = serie[i].Value;
        }

        return new Series(
            index,
            new[] { values },               
            new[] { "Value" }

        );
    }

    public Matrix ToMatrix(int? excludeCol = null)
    {
        int rows = Index.Length;
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

    public static Series FromCsv(
string path,
char sep = ',',
string targetColumn = null,
string groupColumn = null)
    {
        var lines = File.ReadAllLines(path);
        var header = lines[0].Split(sep);
     
        if (header.Length < 2)
            throw new InvalidOperationException("CSV must have at least one feature column");

       
        int targetIdx = targetColumn != null
            ? Array.IndexOf(header, targetColumn)
            : -1;

        int groupIdx = groupColumn != null
            ? Array.IndexOf(header, groupColumn)
            : -1;

        var featureIdx = header
            .Select((h, i) => i)
            .Where(i => i != targetIdx && i != groupIdx)
            .ToArray();

        int rows = lines.Length - 1;
        int cols = featureIdx.Length;

        var data = new double[cols][];
        for (int c = 0; c < cols; c++)
            data[c] = new double[rows];

        var groups = groupIdx >= 0 ? new int[rows] : null;
        var groupMap = new Dictionary<string, int>();
        int nextGroupId = 0;

        for (int r = 0; r < rows; r++)
        {
            var parts = lines[r + 1].Split(sep);

            for (int c = 0; c < cols; c++)
                data[c][r] = double.Parse(parts[featureIdx[c]]);

            if (groupIdx >= 0)
            {
                var key = parts[groupIdx];
                if (!groupMap.TryGetValue(key, out var id))
                    groupMap[key] = id = nextGroupId++;

                groups[r] = id;
            }
        }

        var index = Enumerable.Range(0, rows).ToArray();
        return new Series(index, data, header.Skip(1).ToArray(), groups);
    }
}
