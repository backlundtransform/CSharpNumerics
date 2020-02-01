using System;

namespace Numerics.Objects
{
    public struct Matrix
    {

        public double[,] values;
        public double[,] identity;
        public int rowLength;
        public int columnLength;

        public Matrix(double[,] matrix)
        {
            rowLength = matrix.GetLength(0);
            columnLength = matrix.GetLength(1);
            identity = new double[rowLength, columnLength];

            for (var i = 0; i < rowLength; i++)
            {

                for (var j = 0; j < columnLength; j++)
                {

                    if (i == j)
                    {
                        identity[i, j] = 1;
                    }

                    if (i != j)
                    {
                        identity[i, j] = 0;

                    }

                }
            }

            values = matrix;
        }

        public Matrix Inverse()
        {
            var determinant =Determinant();
            if (determinant == 0) {
                throw new Exception("This matrix is not invertible");
            }
            var adj = Adjugate();
            return adj / determinant;
        }


        public Matrix Transpose()
        {

            var matrix = new double[columnLength, rowLength];

            for (var i = 0; i < columnLength; i++)
            {
                for (var j = 0; j < rowLength; j++)
                {
                    matrix[i, j] = values[j, i];
                }

            }

            return new Matrix(matrix);
        }


        public Matrix Pascal()
        {
            var pascal = new double[rowLength,columnLength];

            for (var i = 0; i < rowLength; i++)
            {

                for (var j = 0; j < columnLength; j++)
                {
                    pascal[i, j] = Pascal(i + 1, j + 1);

                }
            }
            return new Matrix(pascal);

        }

        private int Pascal(int row, int column)
        {
            if (column == 0 || column > row)
            {
                return 0;
            }

            if (row == 1 && column == 1)
            {
                return 1;
            }

            return (Pascal(row - 1, column - 1) + Pascal(row - 1, column));

        }

        public double Determinant()
        {
            if (rowLength == columnLength)
            {
                return Determinant(values, rowLength);
            }
            return 0;

        }

        private double Determinant(double[,] matrix, int n)
        {
            var determinant = 0.0;
            if (n == 1)
            {
                return matrix[0, 0];
            }

            var temp = new double[rowLength, columnLength];

            var sign = 1;

            for (var i = 0; i < n; i++)
            {

                temp = GetCofactor(matrix, temp, 0, i, n);
                determinant += sign * matrix[0, i] * Determinant(temp, n - 1);
                sign = -sign;
            }
            return determinant;
        }

        public Matrix Adjugate()
        {

            if (columnLength != rowLength)
            {
                throw new Exception("Is not a NxN matrix");
            }
            var adj = new double[rowLength, columnLength];
            if (rowLength == 1)
            {
                adj[0, 0] = 1;
                return new Matrix(adj);
            }
            var sign = 1;
            var temp = new double[rowLength, columnLength];

            for (var i = 0; i < rowLength; i++)
            {
                for (var j = 0; j < columnLength; j++)
                {

                    temp = GetCofactor(values, temp, i, j, rowLength);

                    sign = ((i + j) % 2 == 0) ? 1 : -1;

                    adj[j, i] = (sign) * (Determinant(temp, rowLength - 1));
                }
            }
            return new Matrix(adj);
        }

        public static Matrix operator *(Matrix a, Matrix b) => a.GetMultiplicationResult(b);
        public static Vector operator *(Matrix a, Vector b) => a.GetMultiplicationResult(b);
        public static Matrix operator *(double a, Matrix b) => b.GetResult(b, 0, a);
        public static Matrix operator /(Matrix a, double b) => a.GetResult(a, 0, 1/b);
        public static Matrix operator -(Matrix a, Matrix b) => a.GetResult(b, -1);
        public static Matrix operator +(Matrix a, Matrix b) => a.GetResult(b, 1);


        private double[,] GetCofactor(double[,] matrix, double[,] temp, int rowIndex, int columnIndex, int length)
        {
            var i = 0;
            var j = 0;

            for (var row = 0; row < length; row++)
            {
                for (var col = 0; col < length; col++)
                {
                    if (row != rowIndex && col != columnIndex)
                    {
                        temp[i, j++] = matrix[row, col];

                        if (j == length - 1)
                        {
                            j = 0;
                            i++;
                        }
                    }
                }
            }

            return temp;
        }

        private Matrix GetResult(Matrix b, int sign, double multiplier = 1)
        {
            if (rowLength != b.rowLength && columnLength != b.columnLength)
            {
                throw new Exception("The matrixes does not have the same size");
            }
            var result = new double[rowLength, columnLength];

            for (var j = 0; j < columnLength; j++)
            {
                for (var i = 0; i < rowLength; i++)
                {
                    result[i, j] = multiplier * values[i, j] + sign * b.values[i, j];
                }
            }

            return new Matrix(result);

        }


        private Matrix GetMultiplicationResult(Matrix b)
        {
            if (columnLength != b.rowLength)
            {
                throw new Exception("The column length of the first matrix is not same as row length of the second matrix");
            }
            var result = new double[rowLength, b.columnLength];


            for (var i = 0; i < rowLength; i++)
            {
                for (var j = 0; j < b.columnLength; j++)
                {
                    result[i, j] = 0;
                    for (var k = 0; k < b.rowLength; k++)
                    {
                        result[i, j] += values[i, k] * b.values[k, j];

                    }

                }
            }

            return new Matrix(result);

        }

        private Vector GetMultiplicationResult(Vector b)
        {
            if (columnLength != 3)
            {
                throw new Exception("The column length of the matrix should be 3");
            }
               return new Vector(values[0, 0] * b.x + values[0, 1] * b.y + values[0, 2] * b.z,
                values[1, 0] * b.x + values[1, 1] * b.y + values[1, 2] * b.z,
                values[2, 0] * b.x + values[2, 1] * b.y + values[2, 2] * b.z); 

        }
    }
}
