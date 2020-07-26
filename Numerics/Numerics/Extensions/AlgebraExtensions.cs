using System.Collections.Generic;

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


    }
}
