using CSharpNumerics.Objects;
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

        public static double EulerMetod(this Func<(double t, double y), double> func, double min, double max, double stepSize, double yInitial)
        {

            var y = yInitial;

            for (var t = min + stepSize; t <= max; t += stepSize)
            {
              
                y += stepSize * func((t,y));

            }
            return y;

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

        public static VectorN LinearSystemSolver(this Matrix matrix, VectorN vector)
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


        public static List<double> EigenVector(this Matrix matrix, double egienValue)
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
            var vector = Enumerable.Range(1, matrix.rowLength).Select(x => 1 + 0.00001 * x).ToList();

            for (var i = 0; i < 5; i++)
            {
                vector = matrix.Inverse() * vector;

            }
            var min = vector.Where(p => p != 0).Min(p => Math.Abs(p));

            return vector.Select(c => Math.Abs(Math.Round(c / min))).ToList();
        }

        public static List<double> DominantEigenVector(this Matrix matrix)
        {

            var vector = Enumerable.Range(1, matrix.rowLength).Select(x => 1 + 0.00001 * x).ToList();

            for (var i = 0; i < 5; i++)
            {
                vector = matrix * vector;

            }
            var min = vector.Where(p => p != 0).Min(p => Math.Abs(p));

            return vector.Select(c => Math.Abs(Math.Round(c / min))).ToList();
                
         }
        public static List<Func<double, double>> OdeSolver(this Matrix matrix, double tZero)
        {
            var list = Enumerable.Range(1, matrix.rowLength).Select(x => tZero).ToList();
          

            return matrix.OdeSolver(list);
        }

        public static List<Func<double,double>> OdeSolver(this Matrix matrix, List<double> tZeros)
        {
            var eigenValues = matrix.EigenValues();
            var eigenVectors = new Dictionary<double, List<double>>();

            eigenValues.Reverse();

            foreach (var value in eigenValues) {
                var vector = matrix.EigenVector(value);
                eigenVectors.Add(value, vector);
            }

            var eigenMatrix = matrix;

            for (var i = 0; i < matrix.rowLength; i++)
            {

                for (var j = 0; j < matrix.columnLength; j++)
                {
                    if (eigenValues.Count > i)
                    {
                        eigenMatrix.values[i, j] = eigenVectors[eigenValues[j]][i];
                    }
    
                }
            }

            var constants = eigenMatrix.GaussElimination(tZeros);

            double fx(double t) => GetFunc(0, constants, eigenValues, eigenVectors, t);
            double fy(double t) => GetFunc(1, constants, eigenValues, eigenVectors, t);
            double fz(double t) => GetFunc(2, constants, eigenValues, eigenVectors, t);

            return new List<Func<double, double>>() { fx, fy, fz };
        }

        private static double GetFunc(int index, List<double> constants, List<double> eigenValues, Dictionary<double, List<double>> eigenVectors, double t)
        {
            return Math.Exp(eigenValues[0] * t) * eigenVectors[eigenValues[0]][index] * constants[0] +
                                       Math.Exp(eigenValues[1] * t) * eigenVectors[eigenValues[1]][index] * constants[1] +
                                          (constants.Count>2? Math.Exp(eigenValues[2] * t) * eigenVectors[eigenValues[2]][index] * constants[2]:0.0);
        }

        public static List<double> EigenValues(this Matrix matrix)
        {
            var results = new List<double>();
            var dominantEigenvector = matrix.DominantEigenVector();

            results.Add((matrix * dominantEigenvector).Dot(dominantEigenvector) / dominantEigenvector.Dot(dominantEigenvector));

            for (var s = 0; s < matrix.rowLength-1; s++)
            {
     
                var MatrixB = new Matrix(matrix.identity);

                var MatrixA = new Matrix(matrix.values);


                for (var j = 0; j < MatrixB.rowLength; j++)
                {
                    for (var i = 0; i < MatrixB.columnLength; i++)
                    {

                        MatrixB.values[j, i] = MatrixA.values[j, i] - dominantEigenvector[j] * MatrixA.values[s, i] / dominantEigenvector[s];
                    }

                }
                dominantEigenvector = MatrixB.DominantEigenVector();

                results.Add((MatrixB * dominantEigenvector).Dot(dominantEigenvector) / dominantEigenvector.Dot(dominantEigenvector));
            }

         

            return results;

        }

        public static double DominantEigenValue(this Matrix matrix) {

            var result = matrix.DominantEigenVector();

            var eigenvalue = (matrix * result).Dot(result) / result.Dot(result);

            return eigenvalue;
        
        }

       

    }
}
