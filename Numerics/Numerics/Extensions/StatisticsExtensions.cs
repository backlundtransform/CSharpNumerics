using System.Collections.Generic;
namespace System.Linq
{

    public static class StatisticsExtensions
    {
   
        public static double Median<T>(this IEnumerable<T> enumerable, Func<T, double> func)
        {

            var count = enumerable.Count();

            if (count == 0)
            {
                return 0;
            }

            var list = enumerable.OrderBy(func).Select(func);
            var midpoint = count / 2;

            if (count % 2 == 0)
            {
                return (Convert.ToDouble(list.ElementAt(midpoint - 1)) + Convert.ToDouble(list.ElementAt(midpoint))) / 2.0;
            }

            return Convert.ToDouble(list.ElementAt(midpoint));


        }

    }
}
