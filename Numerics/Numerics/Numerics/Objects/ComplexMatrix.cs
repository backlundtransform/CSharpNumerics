using System;

namespace CSharpNumerics.Numerics.Objects
{
    public struct ComplexMatrix
    {

        public ComplexNumber[,] values;
        public int rowLength;
        public int columnLength;

        public ComplexMatrix(ComplexNumber[,] matrix)
        {
            rowLength = matrix.GetLength(0);
            columnLength = matrix.GetLength(1);
            values = matrix;
        }

        public ComplexMatrix(int rows, int cols)
        {
            rowLength = rows;
            columnLength = cols;
            values = new ComplexNumber[rows, cols];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    values[i, j] = new ComplexNumber(0, 0);
        }

        public ComplexMatrix Transpose()
        {
            var matrix = new ComplexNumber[columnLength, rowLength];

            for (var i = 0; i < columnLength; i++)
                for (var j = 0; j < rowLength; j++)
                    matrix[i, j] = values[j, i];

            return new ComplexMatrix(matrix);
        }

        public ComplexMatrix ConjugateTranspose()
        {
            var matrix = new ComplexNumber[columnLength, rowLength];

            for (var i = 0; i < columnLength; i++)
                for (var j = 0; j < rowLength; j++)
                    matrix[i, j] = values[j, i].GetConjugate();

            return new ComplexMatrix(matrix);
        }

        public ComplexNumber Determinant()
        {
            if (rowLength == columnLength)
                return Determinant(values, rowLength);
            return new ComplexNumber(0, 0);
        }

        private ComplexNumber Determinant(ComplexNumber[,] matrix, int n)
        {
            var determinant = new ComplexNumber(0, 0);
            if (n == 1)
                return matrix[0, 0];

            var temp = new ComplexNumber[rowLength, columnLength];
            for (int r = 0; r < rowLength; r++)
                for (int c = 0; c < columnLength; c++)
                    temp[r, c] = new ComplexNumber(0, 0);

            var sign = 1;

            for (var i = 0; i < n; i++)
            {
                temp = GetCofactor(matrix, temp, 0, i, n);
                var signComplex = new ComplexNumber(sign, 0);
                determinant = determinant + signComplex * matrix[0, i] * Determinant(temp, n - 1);
                sign = -sign;
            }
            return determinant;
        }

        public ComplexMatrix Adjugate()
        {
            if (columnLength != rowLength)
                throw new Exception("Is not a NxN matrix");

            var adj = new ComplexNumber[rowLength, columnLength];
            if (rowLength == 1)
            {
                adj[0, 0] = new ComplexNumber(1, 0);
                return new ComplexMatrix(adj);
            }

            var temp = new ComplexNumber[rowLength, columnLength];
            for (int r = 0; r < rowLength; r++)
                for (int c = 0; c < columnLength; c++)
                    temp[r, c] = new ComplexNumber(0, 0);

            for (var i = 0; i < rowLength; i++)
            {
                for (var j = 0; j < columnLength; j++)
                {
                    temp = GetCofactor(values, temp, i, j, rowLength);

                    var sign = ((i + j) % 2 == 0) ? 1 : -1;

                    adj[j, i] = new ComplexNumber(sign, 0) * Determinant(temp, rowLength - 1);
                }
            }
            return new ComplexMatrix(adj);
        }

        public ComplexMatrix Inverse()
        {
            var determinant = Determinant();
            if (determinant.GetMagnitude() == 0)
                throw new Exception("This matrix is not invertible");
            var adj = Adjugate();
            return adj / determinant;
        }


        public static ComplexMatrix operator +(ComplexMatrix a, ComplexMatrix b)
        {
            if (a.rowLength != b.rowLength || a.columnLength != b.columnLength)
                throw new Exception("The matrices do not have the same size");
            var result = new ComplexNumber[a.rowLength, a.columnLength];
            for (var i = 0; i < a.rowLength; i++)
                for (var j = 0; j < a.columnLength; j++)
                    result[i, j] = a.values[i, j] + b.values[i, j];
            return new ComplexMatrix(result);
        }

        public static ComplexMatrix operator -(ComplexMatrix a, ComplexMatrix b)
        {
            if (a.rowLength != b.rowLength || a.columnLength != b.columnLength)
                throw new Exception("The matrices do not have the same size");
            var result = new ComplexNumber[a.rowLength, a.columnLength];
            for (var i = 0; i < a.rowLength; i++)
                for (var j = 0; j < a.columnLength; j++)
                    result[i, j] = a.values[i, j] - b.values[i, j];
            return new ComplexMatrix(result);
        }

        public static ComplexMatrix operator *(ComplexMatrix a, ComplexMatrix b)
        {
            if (a.columnLength != b.rowLength)
                throw new Exception("The column length of the first matrix is not same as row length of the second matrix");
            var result = new ComplexNumber[a.rowLength, b.columnLength];

            for (var i = 0; i < a.rowLength; i++)
            {
                for (var j = 0; j < b.columnLength; j++)
                {
                    result[i, j] = new ComplexNumber(0, 0);
                    for (var k = 0; k < b.rowLength; k++)
                        result[i, j] = result[i, j] + a.values[i, k] * b.values[k, j];
                }
            }

            return new ComplexMatrix(result);
        }

        public static ComplexMatrix operator *(ComplexNumber a, ComplexMatrix b)
        {
            var result = new ComplexNumber[b.rowLength, b.columnLength];
            for (var i = 0; i < b.rowLength; i++)
                for (var j = 0; j < b.columnLength; j++)
                    result[i, j] = a * b.values[i, j];
            return new ComplexMatrix(result);
        }

        public static ComplexMatrix operator *(double a, ComplexMatrix b) =>
            new ComplexNumber(a, 0) * b;

        public static ComplexMatrix operator /(ComplexMatrix a, ComplexNumber b)
        {
            var result = new ComplexNumber[a.rowLength, a.columnLength];
            for (var i = 0; i < a.rowLength; i++)
                for (var j = 0; j < a.columnLength; j++)
                    result[i, j] = a.values[i, j] / b;
            return new ComplexMatrix(result);
        }

        public static ComplexMatrix operator /(ComplexMatrix a, double b) =>
            a / new ComplexNumber(b, 0);

        public static ComplexVector operator *(ComplexMatrix a, ComplexVector b)
        {
            if (a.columnLength == 2)
                return new ComplexVector(
                    a.values[0, 0] * b.x + a.values[0, 1] * b.y,
                    a.values[1, 0] * b.x + a.values[1, 1] * b.y,
                    new ComplexNumber(0, 0));

            if (a.columnLength != 3)
                throw new Exception("The column length of the matrix should be 3");

            return new ComplexVector(
                a.values[0, 0] * b.x + a.values[0, 1] * b.y + a.values[0, 2] * b.z,
                a.values[1, 0] * b.x + a.values[1, 1] * b.y + a.values[1, 2] * b.z,
                a.values[2, 0] * b.x + a.values[2, 1] * b.y + a.values[2, 2] * b.z);
        }

        public static ComplexVectorN operator *(ComplexMatrix a, ComplexVectorN b)
        {
            if (a.columnLength != b.Length)
                throw new Exception("Matrix column length must match vector length.");

            var result = new ComplexNumber[a.rowLength];

            for (int i = 0; i < a.rowLength; i++)
            {
                result[i] = new ComplexNumber(0, 0);
                for (int j = 0; j < a.columnLength; j++)
                    result[i] = result[i] + a.values[i, j] * b[j];
            }

            return new ComplexVectorN(result);
        }

        public static implicit operator ComplexMatrix(Matrix m)
        {
            var result = new ComplexNumber[m.rowLength, m.columnLength];
            for (int i = 0; i < m.rowLength; i++)
                for (int j = 0; j < m.columnLength; j++)
                    result[i, j] = new ComplexNumber(m.values[i, j], 0);
            return new ComplexMatrix(result);
        }


        private ComplexNumber[,] GetCofactor(ComplexNumber[,] matrix, ComplexNumber[,] temp, int rowIndex, int columnIndex, int length)
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
    }
}
