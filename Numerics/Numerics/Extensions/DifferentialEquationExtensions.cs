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

    }
}
