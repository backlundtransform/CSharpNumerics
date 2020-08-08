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
            var value = 0.0;

            for (var j = 0; j < 100; j++)
            {
                var y = func(xZero);
                var yPrime = func.Derivate(xZero);
                if (Math.Abs(yPrime) <0.0000000000001) {
                    break;
                }
                value = xZero - y / yPrime;

                xZero = value;


            }
            return value;

        }

    }
}
