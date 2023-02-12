using System.Collections.Generic;


namespace System
{
    public static class ArithmeticExtensions
    {
        public static bool IsPrime(this int n)
        {
            if (n <= 1)
            {
                return false;
            }

            for (int i = 2; i < n; i++)
            {
                if (n % i == 0)
                {
                    return false;
                }
            }

            return true;
        }

        public static List<int> GetPrimeFactors(this int n)
        {
            var factors = new List<int>();
            for (int i = 2; i <= n / i; i++)
            {
                while (n % i == 0)
                {
                    factors.Add(i);
                    n /= i;
                }
            }
            if (n > 1)
            {
                factors.Add(n);
            }
            return factors;
        }
    }
}
