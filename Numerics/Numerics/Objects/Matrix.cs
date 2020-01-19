namespace Numerics.Objects
{
    public struct Matrix
    {

        public double[,] Values;
        public double[,] Identity;

        public Matrix(double[,] values)
        {

            var totalLength = values.Length;
            var rowLength = values.GetLength(0);
            var columnLength = totalLength / rowLength;
            var identity = new double[rowLength, columnLength];

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

            Values = values;
            Identity = identity;
        }
    }
}
