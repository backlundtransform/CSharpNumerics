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


        public static Vector GaussElimination(this Matrix matrix, Vector vector)
        {
            var values = matrix.GaussElimination(new List<double>() { vector.x,vector.y, vector.z });

            return new Vector(values[0], values[1], values[2]);
        }


        public static List<double>  GaussElimination(this Matrix matrix, List<double> vector)
        {
            var n = vector.Count;

            for (var p = 0; p < n; p++)
            {

                var max = p;
                for (var i = p + 1; i < n; i++)
                {
                    if (Math.Abs(matrix.values[i,p]) > Math.Abs(matrix.values[max,p]))
                    {
                        max = i;
                    }
                }
      

                for (var c = 0; c < matrix.columnLength; c++) {
                    var temp = matrix.values[p, c];

                    matrix.values[p, c] = matrix.values[max, c];

                    matrix.values[max, c] = temp;
                }

               var t = vector[p];
                vector[p] = vector[max];
                vector[max] = t;

         
                for (var i = p + 1; i < n; i++)
                {
                    var alpha = matrix.values[i,p] / matrix.values[p,p];
                    vector[i] -= alpha * vector[p];
                    for (var j = p; j < n; j++)
                    {
                        matrix.values[i,j] -= alpha * matrix.values[p,j];
                    }
                }
            }

            var x = new double[n];
            for (var i = n - 1; i >= 0; i--)
            {
                var sum = 0.0;
                for (var j = i + 1; j < n; j++)
                {
                    sum += matrix.values[i,j] * x[j];
                }
               if( matrix.values[i, i] != 0)
                {

                    if((vector[i] - sum) != 0)
                    {
                        x[i] = (vector[i] - sum) / matrix.values[i, i];
                    }

                    if ((vector[i] - sum) == 0)
                    {
                        x[i] = matrix.values[i, i];
                    }

                }
              
            }
            return x.ToList();

     
        }


        public static Vector EigenVector(this Matrix matrix, double egienValue)
        {

            for (var i = 0; i < matrix.rowLength; i++)
            {

                for (var j = 0; j < matrix.columnLength; j++)
                {

                    if (i == j)
                    {
                        matrix.values[i, j] -= egienValue-0.1;
                    } 

                }
            }
            var vector = new Vector(1, 1, 1);

            for (var i = 0; i < 5; i++)
            {
                vector = matrix.Inverse() * vector;

            }
            var min = new List<double>() { Math.Abs(vector.x), Math.Abs(vector.y), Math.Abs(vector.z) }.Where(p => p != 0).Min(p => p);

            return new Vector(Math.Abs(Math.Round(vector.x/min)), Math.Abs(Math.Round(vector.y/min)), Math.Abs(Math.Round(vector.z/min)));
        }


        public static List<Func<double,double>> OdeSolver(this Matrix matrix, double tZero)
        {
            var eigenValues = matrix.EigenValues();
            var eigenVectors = new Dictionary<double, List<double>>();

            foreach (var value in eigenValues) {
                var vector = matrix.EigenVector(value);
                eigenVectors.Add(value, new List<double>() { vector.x, vector.y, vector.z });
            }

            var eigenMatrix = matrix;

            for (var i = 0; i < matrix.rowLength; i++)
            {

                for (var j = 0; j < matrix.columnLength; j++)
                {
                
                    eigenMatrix.values[i, j] = eigenVectors[eigenValues[j]][i]; 

                }
            }

            var constants = eigenMatrix.GaussElimination(new List<double>() { tZero, tZero, tZero });

            double fx(double t) => GetFunc(0, constants, eigenValues, eigenVectors, t);
            double fy(double t) => GetFunc(1, constants, eigenValues, eigenVectors, t);
            double fz(double t) => GetFunc(2, constants, eigenValues, eigenVectors, t);

            return new List<Func<double, double>>() { fx, fy, fz };
        }

        private static double GetFunc(int index, List<double> constants, List<double> eigenValues, Dictionary<double, List<double>> eigenVectors, double t)
        {
            return Math.Exp(eigenValues[0] * t) * eigenVectors[eigenValues[0]][index] * constants[0] +
                                       Math.Exp(eigenValues[1] * t) * eigenVectors[eigenValues[1]][index] * constants[1] +
                                           Math.Exp(eigenValues[2] * t) * eigenVectors[eigenValues[2]][index] * constants[2];
        }

        public static List<double> EigenValues(this Matrix matrix)
        {


            var results = new List<double>();

            for (var i = -100.0; i <= 100; i ++)
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
