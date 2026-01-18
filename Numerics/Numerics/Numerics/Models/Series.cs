using Numerics.Models;
using Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Text;

namespace CSharpNumerics.Numerics.Models
{
    public class Series(int[] index, double[][] data)
    {
        public int[] Index { get; } = index;
        public double[][] Data { get; } = data;

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
                new[] { values }
            );
        }

        public Matrix ToMatrix()
        {
            int rows = Index.Length;
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
}
