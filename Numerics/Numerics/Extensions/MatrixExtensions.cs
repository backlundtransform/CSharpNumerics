namespace Numerics.Objects
{
    public static class MatrixExtensions
    {
        public static Matrix Transpose(this Matrix matrix)
        {
            var totalLength = matrix.Values.Length;
            var rowLength = matrix.Values.GetLength(0);
            var columnLength = totalLength / rowLength;
          
            var values = new double[columnLength, rowLength];

            for (var i = 0; i < columnLength; i++)
            {
                
                for (var j = 0; j < rowLength  ; j++)
                {
                    values[i, j] = matrix.Values[j, i];
                }

            }

            return new Matrix(values);
        }
    }
 
}
