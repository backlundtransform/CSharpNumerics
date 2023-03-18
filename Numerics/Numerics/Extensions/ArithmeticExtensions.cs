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

        public static bool IsHappy(this int n)
        {
            var visited = new HashSet<int>();
            while (n != 1 && !visited.Contains(n))
            {
                visited.Add(n);
                int sum = 0;
                while (n != 0)
                {
                    int digit = n % 10;
                    sum += digit * digit;
                    n /= 10;
                }
                n = sum;
            }
            return n == 1;
        }

        public static bool IsPerfectNumber(this int number)
        {
            var sum = 1;

            for (var i = 2; i <= Math.Sqrt(number); i++)
            {
                if (number % i == 0)
                {
                    sum += i;
                    sum += number / i;
                }
            }

            return (sum == number && number != 1);

        }

        public static int GetDecimalPlaces(this double n)
        {
            n = Math.Abs(n);
            n -= (int)n;
            var decimalPlaces = 0;
            while (n > 0)
            {
                decimalPlaces++;
                n *= 10;
                n -= (int)n;
            }
            return decimalPlaces;
        }
    }
}
