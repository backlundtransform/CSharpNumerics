using Numerics.Models;
using Numerics.Objects;
using System.Collections.Generic;
using System.Linq;

namespace System
{
    public static class TransformExtensions
    {
        private const int defaultNvalue= 14;
 

        public static List<Serie> ToFrequencyResolution(this List<ComplexNumber> numbers, double samplingFrequency, int pow=1)
        {
            return numbers.Select((p, index)=>new Serie() { Value=Math.Pow(p.GetMagnitude(),pow) / numbers.Count(), Index = index* samplingFrequency/numbers.Count() }).ToList();
        }


        public static List<ComplexNumber> FastFourierTransform(this Func<double,double> func,double minValue, double maxValue, double stepSize)
        {
         return func.GetSeries(minValue, maxValue, stepSize).Select(p => new ComplexNumber(p.Value, 0)).FastFourierTransform(-1);
        }

        public static List<ComplexNumber> DiscreteFourierTransform(this Func<double, double> func, double minValue, double maxValue, double stepSize)
        {
            return func.GetSeries(minValue, maxValue, stepSize).Select(p => new ComplexNumber(p.Value, 0)).FastFourierTransform(-1);
        }


        public static List<ComplexNumber> FastFourierTransform(this List<ComplexNumber> numbers)
        {
            return numbers.FastFourierTransform(-1);
        }

        public static List<ComplexNumber> DiscreteFourierTransform(this List<ComplexNumber> numbers)
        {
            return  numbers.DiscreteFourierTransform(-1);
        }

        public static List<ComplexNumber> FastFourierTransform(this List<double> numbers)
        {
            return numbers.Select(p=>new ComplexNumber(p,0)).FastFourierTransform(-1);
        }

        public static List<ComplexNumber> DiscreteFourierTransform(this List<double> numbers)
        {
            return numbers.Select(p => new ComplexNumber(p, 0)).DiscreteFourierTransform(-1);
        }


        public static List<ComplexNumber> InverseFastFourierTransform(this Func<double, double> func, double minValue, double maxValue, double stepSize)
        {
            return func.GetSeries(minValue, maxValue, stepSize).Select(p => new ComplexNumber(p.Value, 0)).FastFourierTransform(1).Select(p => new ComplexNumber(p.realPart / stepSize, p.imaginaryPart)).ToList();
        }

        public static List<ComplexNumber> InverseDiscreteFourierTransform(this Func<double, double> func, double minValue, double maxValue, double stepSize)
        {
            return func.GetSeries(minValue, maxValue, stepSize).Select(p => new ComplexNumber(p.Value, 0)).DiscreteFourierTransform(1).Select(p => new ComplexNumber(p.realPart / stepSize, p.imaginaryPart)).ToList();
        }


        public static List<ComplexNumber> InverseFastFourierTransform(this List<ComplexNumber> numbers)
        {
            return numbers.FastFourierTransform(1).Select(p => new ComplexNumber(p.realPart / numbers.Count(), p.imaginaryPart)).ToList();
        }

        public static List<ComplexNumber> InverseDiscreteFourierTransform(this List<ComplexNumber> numbers)
        {
            return numbers.DiscreteFourierTransform(1).Select(p => new ComplexNumber(p.realPart / numbers.Count(), p.imaginaryPart)).ToList();
        }

        public static List<ComplexNumber> InverseFastFourierTransform(this List<double> numbers)
        {
            return numbers.Select(p => new ComplexNumber(p, 0)).FastFourierTransform(1).Select(p => new ComplexNumber(p.realPart / numbers.Count(), p.imaginaryPart)).ToList();
        }

        public static List<ComplexNumber> InverseDiscreteFourierTransform(this List<double> numbers)
        {
            return numbers.Select(p => new ComplexNumber(p, 0)).DiscreteFourierTransform(1).Select(p => new ComplexNumber(p.realPart / numbers.Count(), p.imaginaryPart)).ToList();
        }


        public static List<ComplexNumber> DiscreteFourierTransform(this IEnumerable<ComplexNumber> input, int sign)
        {

            var numbers = new List<ComplexNumber>(input);
            var result= new List<ComplexNumber>();
            for (var k = 0; k < input.Count(); k++)
            {
                result.Add(new ComplexNumber(0, 0));
                for (var n = 0; n < input.Count(); n++)
                {
                    var exp = new ComplexNumber(0, sign * 2 * Math.PI * k*n / input.Count());
                    exp = exp.Exponential();

                    exp *= numbers[n];
                    result[k] += exp;
                }
            }
            return result;
        }


        public static List<ComplexNumber> FastFourierTransform(this IEnumerable<ComplexNumber> input, int sign)
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

        public static double LaplaceTransform(this Func<double, double> func, double s)
        {
            Func<double, double> funcLaplace=(double u) => Math.Pow(u, s - 1) * func(-Math.Log(u));

           return funcLaplace.Integrate(0,1);
        }


        public static double InverseLaplaceTransform(this Func<double, double> func, double t)
        {
            var numbers = defaultNvalue;
            var result = 0.0;
            var sign = -1;

            for (var i = 0; i < numbers; i++)
            {
                var v = 0.0;
                sign = -sign;

                for (int k = (i + 2) / 2; k <= Math.Min(i + 1, numbers / 2); k++)
                {
                    v += (Math.Pow(k, 7)) * (2 * k).Factorial()
                     / ((2 * k - i - 1).Factorial() * (7 - k).Factorial()
                    * (k - 1).Factorial() * (i + 1 - k).Factorial() * k.Factorial());
                }

                result += v * sign * func(i * Math.Log(2) / t);

            }
            return Math.Log(2.0) / t * result / 2;

        }
   

        public static List<double> LowPassFilter(this List<double> input, List<double> output,double alpha = 0.25)
        {
            if (output == null)
            {
                return input.ToList();
            }
            for (int i = 0; i < input.Count(); i++)
            {
                output[i] = output[i] + alpha * (input[i] - output[i]);
            }
            return output;
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
            int reversed = n;
            int count = bits - 1;

            n >>= 1;
            while (n > 0)
            {
                reversed = (reversed << 1) | (n & 1);
                count--;
                n >>= 1;
            }

            return ((reversed << count) & ((1 << bits) - 1));
        }

    }
}
