using Numerics.Objects;
using System.Collections.Generic;
using System.Linq;

namespace System
{
    public static class DifferentialEquationExtensions
    {

        public static double RungeKutta(this Func<(double t, double y), double> func, double min, double max, double stepSize, double yInitial)
        {

            return func.RungeKutta(min, max, stepSize, yInitial, new Matrix(new double[,] { { 0, 0, 0 }, { 0.5, 0, 0 }, { 0, 0.5, 0 }, { 0, 0, 1 } }), new double[] { 1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0 }, new double[] { 0.0, 0.5, 0.5, 1 });

        }


        public static double RungeKutta(this Func<(double t, double y), double> func, double min, double max, double stepSize, double yInitial, Matrix rungeKuttaMatrix  , double[] weights, double[] nodes)
        {

            var y = yInitial;

            for (var t = min + stepSize; t <= max; t += stepSize)
            {
                var kList = new List<double>();


                for (var i = 0; i < nodes.Length; i++)
                {
                    var value = 0.0;

                    var columnLength = rungeKuttaMatrix.values.GetLength(1);

                   for (var j = 0; j < columnLength; j++)
                    {
                        value += rungeKuttaMatrix.values[i,j]*(j>= kList.Count?0: kList[j]);

                    }
                    kList.Add(func((t + nodes[i] * stepSize, y + stepSize * value)));

                }


                y += stepSize * kList.Sum(p => p * weights[kList.IndexOf(p)]);

            }
            return y;

        }

        public static double TrapezoidalRule(this Func<(double t, double y), double> func, double min, double max, double stepSize, double yInitial)
        {

            var y = yInitial;


            for (var i = min; i <=max; i += stepSize)
            {
                y += stepSize * (func((i, y)) + func((i + stepSize, y + stepSize))) / 2;
         
            }
            return y;
        }

        public static Vector LinearSystemSolver(this Matrix matrix, Vector vector)
        {
            var values = matrix.Inverse() * vector;


            return values;
        }

        public static List<double> LinearSystemSolver(this Matrix matrix, List<double> vector)
        {
            var values = matrix.Inverse() * vector;


            return values;
        }

 

        public static List<double> EigenValues(this Matrix matrix)
        {

            var largest = matrix.LargestEigenValue();

            var smallest = matrix.SmallestEigenValue();

            var results = new List<double>();

            for (var i = smallest; i <= largest; i ++)
            {
           
                if (Math.Round((matrix - i * new Matrix(matrix.identity)).Determinant()) == 0.00)
                {

                    results.Add(i);

                }

            }

            return results;

        }


        public static double SmallestEigenValue(this Matrix matrix)
        {

            var largest = matrix.LargestEigenValue();

            var result = (matrix - largest * matrix.Inverse()).LargestEigenValue();

            return result-largest;


        }


        public static double LargestEigenValue(this Matrix matrix)
        {

            var results = new List<double>();
            for (var i = 0; i < matrix.rowLength; i++)
            {
                var result =0.0;
                for (var j = 0; j < matrix.columnLength; j++)
                {
                    result +=matrix.values[i, j];
                }
                results.Add(result);
            }
    

            for (var i = results.Min(); i <= results.Max(); i +=0.1)
            {
                var lambda = Math.Round(i, 2);
                if (Math.Round((matrix - lambda * new Matrix(matrix.identity)).Determinant(), 2) == 0.00) {

                    return lambda;           

                };

            }
            return 0.0;
        }



    }
}
