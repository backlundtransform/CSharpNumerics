using Numerics.Objects;
using System.Collections.Generic;
using System.Drawing;

namespace System
{
   public static class AlgebraExtensions
    {
        public static double Factorial(this int number)
        {
            if (number <= 1)
            {
                return 1;
            }

            return number * (number - 1).Factorial();
        }
        public static double Factorial(this double number)
        {
            return ((int)number).Factorial();
        }

        public static double Dot(this List<double> firstList, List<double> secondList)
        {
            var value = 0.0;

            for (var j = 0; j < firstList.Count; j++)
            {
                value += firstList[j] * secondList[j];

            }
            return value;

        }

        public static double NewtonRaphson(this Func<double,double> func, double xZero =1.0)
        {
            var value = xZero;

            for (var j = 0; j < 100; j++)
            {
                var y = func(value);
                var yPrime = func.Derivate(value);

                value -= y / yPrime;

            }
            return value;

        }
        public static Matrix WithBiasColumn(this Matrix X)
        {
            var result = new double[X.rowLength, X.columnLength + 1];

            for (int i = 0; i < X.rowLength; i++)
            {
               result[i, 0] = 1.0; 
                for (int j = 0; j < X.columnLength; j++)
                     result[i, j + 1] = X.values[i, j];
            }

                return new Matrix(result);
        }
    }
}
