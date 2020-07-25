using System;
using System.Collections.Generic;
using System.Text;

namespace CSharpNumerics.Extensions
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

    }
}
