using Numerics.Objects;
using System.Collections.Generic;
using System.Linq;

namespace System
{
    public static class TransformExtensions
    {

        public static List<ComplexNumber> Fouriertransform(this List<ComplexNumber> numbers)
        {
            return numbers.Fouriertransform(-1);
        }


        public static List<ComplexNumber> Fouriertransform(this List<double> numbers)
        {
            return numbers.Select(p=>new ComplexNumber(p,0)).Fouriertransform(-1);
        }


        public static List<ComplexNumber> InverseFouriertransform(this List<ComplexNumber> numbers)
        {
     
            return numbers.Fouriertransform(1).Select(p=>new ComplexNumber(p.realPart/numbers.Count(),p.imaginaryPart)).ToList();

        }

        public static List<ComplexNumber> Fouriertransform(this IEnumerable<ComplexNumber> input, int sign)
        {
            var numbers = new List<ComplexNumber>(input);

            var bits = (int)Math.Log(numbers.Count(), 2);
            for (var j = 1; j < numbers.Count(); j++)
            {
                var swapPos = BitReverse(j, bits);
                if (swapPos <= j)
                {
                    continue;
                }
                var temp = numbers[j];
                numbers[j] = numbers[swapPos];
                numbers[swapPos] = temp;
            }
            for (var n = 2; n <= numbers.Count(); n <<= 1)
            {
                for (var i = 0; i < numbers.Count(); i += n)
                {
                    for (var k = 0; k < n / 2; k++)
                    {

                        var evenIndex = i + k;
                        var oddIndex = i + k + (n / 2);

                        if (evenIndex >= numbers.Count()|| oddIndex >= numbers.Count()) {
                            continue;
                        }

                            var even = numbers[evenIndex];
                            var odd = numbers[oddIndex];
                            var exp = new ComplexNumber(0, sign * 2 * Math.PI * k / n);
                            exp = exp.Exponential();
                            numbers[evenIndex] = (even + exp * odd);
                            numbers[oddIndex] = (even - exp * odd);
                  
                 
                  

                    }
                }
            }
            return numbers;


        }

        public static ComplexNumber Exponential(this ComplexNumber complex)
        {
            var realPart = Math.Pow(Math.E, complex.realPart) * Math.Cos(complex.imaginaryPart);
            var  imaginaryPart = Math.Sin(complex.imaginaryPart);
            return new ComplexNumber(realPart, imaginaryPart);
        }


        public static ComplexNumber Pow(this ComplexNumber complex, int power)
        {
            var mag = complex.GetMagnitude();
            var arg = complex.GetArgument();
            var realPart = Math.Pow(mag, power) * Math.Cos(power * arg);
            var imaginaryPart = Math.Pow(mag, power) * Math.Sin(power * arg);
            return new ComplexNumber(realPart, imaginaryPart);
        }

        private static int BitReverse(int n, int bits)
        {
            int reversedN = n;
            int count = bits - 1;

            n >>= 1;
            while (n > 0)
            {
                reversedN = (reversedN << 1) | (n & 1);
                count--;
                n >>= 1;
            }

            return ((reversedN << count) & ((1 << bits) - 1));
        }

    }
}
