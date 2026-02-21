using CSharpNumerics.Numerics.Objects;
using System.Collections.Generic;
using System.Linq;

namespace System
{
    public static class DifferentialEquationExtensions
    {
        /// <summary>
        /// Solves a first-order ODE y' = f(t, y) using the classical 4th order Runge-Kutta (RK4) method.
        /// </summary>
        /// <param name="func">Right-hand side function f(t, y).</param>
        /// <param name="min">Initial time.</param>
        /// <param name="max">Final time.</param>
        /// <param name="stepSize">Time step size.</param>
        /// <param name="yInitial">Initial value y(min).</param>
        /// <returns>Approximate y(max) computed by RK4.</returns>
        public static double RungeKutta(this Func<(double t, double y), double> func, double min, double max, double stepSize, double yInitial)
        {
            return func.RungeKutta(
                min,
                max,
                stepSize,
                yInitial,
                new Matrix(new double[,] { { 0, 0, 0 }, { 0.5, 0, 0 }, { 0, 0.5, 0 }, { 0, 0, 1 } }),
                new double[] { 1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0 },
                new double[] { 0.0, 0.5, 0.5, 1 });
        }

        /// <summary>
        /// Solves a first-order ODE y' = f(t, y) using the explicit Euler method.
        /// Update: y_{n+1} = y_n + h * f(t_n, y_n).
        /// </summary>
        /// <param name="func">Right-hand side function f(t, y).</param>
        /// <param name="min">Initial time.</param>
        /// <param name="max">Final time.</param>
        /// <param name="stepSize">Time step size.</param>
        /// <param name="yInitial">Initial value y(min).</param>
        /// <returns>Approximate y(max) computed by Euler's method.</returns>
        public static double EulerMetod(this Func<(double t, double y), double> func, double min, double max, double stepSize, double yInitial)
        {
            var y = yInitial;

            for (var t = min + stepSize; t <= max; t += stepSize)
            {
                y += stepSize * func((t, y));
            }
            return y;
        }

        /// <summary>
        /// Solves a first-order ODE y' = f(t, y) using a general explicit Runge-Kutta method defined by a Butcher tableau.
        /// </summary>
        /// <param name="func">Right-hand side function f(t, y).</param>
        /// <param name="min">Initial time.</param>
        /// <param name="max">Final time.</param>
        /// <param name="stepSize">Time step size.</param>
        /// <param name="yInitial">Initial value y(min).</param>
        /// <param name="rungeKuttaMatrix">Runge-Kutta A-matrix (stage coefficients).</param>
        /// <param name="weights">Runge-Kutta b-vector (stage weights).</param>
        /// <param name="nodes">Runge-Kutta c-vector (stage nodes).</param>
        /// <returns>Approximate y(max) computed by the specified Runge-Kutta scheme.</returns>
        public static double RungeKutta(
            this Func<(double t, double y), double> func,
            double min,
            double max,
            double stepSize,
            double yInitial,
            Matrix rungeKuttaMatrix,
            double[] weights,
            double[] nodes)
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
                        value += rungeKuttaMatrix.values[i, j] * (j >= kList.Count ? 0 : kList[j]);
                    }
                    kList.Add(func((t + nodes[i] * stepSize, y + stepSize * value)));
                }

                y += stepSize * kList.Sum(p => p * weights[kList.IndexOf(p)]);
            }
            return y;
        }

        /// <summary>
        /// Solves a first-order ODE y' = f(t, y) using a trapezoidal-rule style update.
        /// </summary>
        /// <param name="func">Right-hand side function f(t, y).</param>
        /// <param name="min">Initial time.</param>
        /// <param name="max">Final time.</param>
        /// <param name="stepSize">Time step size.</param>
        /// <param name="yInitial">Initial value y(min).</param>
        /// <returns>Approximate y(max) computed by the trapezoidal update.</returns>
        public static double TrapezoidalRule(this Func<(double t, double y), double> func, double min, double max, double stepSize, double yInitial)
        {
            var y = yInitial;

            for (var i = min; i <= max; i += stepSize)
            {
                y += stepSize * (func((i, y)) + func((i + stepSize, y + stepSize))) / 2;
            }
            return y;
        }

        #region Vector ODE Solvers

        /// <summary>
        /// Solves a 3D vector ODE y' = f(t, y) using the classical 4th order Runge-Kutta method.
        /// </summary>
        /// <param name="func">Right-hand side function f(t, y) returning a Vector derivative.</param>
        /// <param name="tMin">Initial time.</param>
        /// <param name="tMax">Final time.</param>
        /// <param name="stepSize">Time step size.</param>
        /// <param name="y0">Initial value y(tMin).</param>
        /// <returns>Approximate y(tMax) computed by RK4.</returns>
        public static Vector RungeKutta(
            this Func<(double t, Vector y), Vector> func,
            double tMin, double tMax, double stepSize, Vector y0)
        {
            var y = new Vector(y0.x, y0.y, y0.z);

            for (var t = tMin; t < tMax; t += stepSize)
            {
                var h = Math.Min(stepSize, tMax - t);

                var k1 = func((t, y));
                var k2 = func((t + 0.5 * h, y + 0.5 * h * k1));
                var k3 = func((t + 0.5 * h, y + 0.5 * h * k2));
                var k4 = func((t + h, y + h * k3));

                y = y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
            }

            return y;
        }

        /// <summary>
        /// Solves a 3D vector ODE y' = f(t, y) using the explicit Euler method.
        /// </summary>
        /// <param name="func">Right-hand side function f(t, y) returning a Vector derivative.</param>
        /// <param name="tMin">Initial time.</param>
        /// <param name="tMax">Final time.</param>
        /// <param name="stepSize">Time step size.</param>
        /// <param name="y0">Initial value y(tMin).</param>
        /// <returns>Approximate y(tMax) computed by Euler's method.</returns>
        public static Vector EulerMethod(
            this Func<(double t, Vector y), Vector> func,
            double tMin, double tMax, double stepSize, Vector y0)
        {
            var y = new Vector(y0.x, y0.y, y0.z);

            for (var t = tMin; t < tMax; t += stepSize)
            {
                var h = Math.Min(stepSize, tMax - t);
                y = y + h * func((t, y));
            }

            return y;
        }

        /// <summary>
        /// Solves a 3D vector ODE and returns the full trajectory as a list of (time, state) pairs.
        /// Uses the classical 4th order Runge-Kutta method.
        /// </summary>
        /// <param name="func">Right-hand side function f(t, y).</param>
        /// <param name="tMin">Initial time.</param>
        /// <param name="tMax">Final time.</param>
        /// <param name="stepSize">Time step size.</param>
        /// <param name="y0">Initial value y(tMin).</param>
        /// <returns>Full trajectory as a list of (t, y) pairs.</returns>
        public static List<(double t, Vector y)> RungeKuttaTrajectory(
            this Func<(double t, Vector y), Vector> func,
            double tMin, double tMax, double stepSize, Vector y0)
        {
            var trajectory = new List<(double t, Vector y)>();
            var y = new Vector(y0.x, y0.y, y0.z);
            trajectory.Add((tMin, y));

            for (var t = tMin; t < tMax; t += stepSize)
            {
                var h = Math.Min(stepSize, tMax - t);

                var k1 = func((t, y));
                var k2 = func((t + 0.5 * h, y + 0.5 * h * k1));
                var k3 = func((t + 0.5 * h, y + 0.5 * h * k2));
                var k4 = func((t + h, y + h * k3));

                y = y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
                trajectory.Add((t + h, y));
            }

            return trajectory;
        }

        #endregion

        #region System ODE Solvers (VectorN)

        /// <summary>
        /// Solves a system of ODEs y' = f(t, y) using the classical 4th order Runge-Kutta method.
        /// Supports arbitrary dimensions — use for dynamics by packing state as [x, y, z, vx, vy, vz].
        /// </summary>
        /// <param name="func">Right-hand side function f(t, y) returning the derivative vector.</param>
        /// <param name="tMin">Initial time.</param>
        /// <param name="tMax">Final time.</param>
        /// <param name="stepSize">Time step size.</param>
        /// <param name="y0">Initial state vector y(tMin).</param>
        /// <returns>Approximate y(tMax) computed by RK4.</returns>
        public static VectorN RungeKutta(
            this Func<(double t, VectorN y), VectorN> func,
            double tMin, double tMax, double stepSize, VectorN y0)
        {
            var y = new VectorN(y0.Values);

            for (var t = tMin; t < tMax; t += stepSize)
            {
                var h = Math.Min(stepSize, tMax - t);

                var k1 = func((t, y));
                var k2 = func((t + 0.5 * h, y + 0.5 * h * k1));
                var k3 = func((t + 0.5 * h, y + 0.5 * h * k2));
                var k4 = func((t + h, y + h * k3));

                y = y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
            }

            return y;
        }

        /// <summary>
        /// Solves a system of ODEs y' = f(t, y) using the explicit Euler method.
        /// </summary>
        /// <param name="func">Right-hand side function f(t, y) returning the derivative vector.</param>
        /// <param name="tMin">Initial time.</param>
        /// <param name="tMax">Final time.</param>
        /// <param name="stepSize">Time step size.</param>
        /// <param name="y0">Initial state vector y(tMin).</param>
        /// <returns>Approximate y(tMax) computed by Euler's method.</returns>
        public static VectorN EulerMethod(
            this Func<(double t, VectorN y), VectorN> func,
            double tMin, double tMax, double stepSize, VectorN y0)
        {
            var y = new VectorN(y0.Values);

            for (var t = tMin; t < tMax; t += stepSize)
            {
                var h = Math.Min(stepSize, tMax - t);
                y = y + h * func((t, y));
            }

            return y;
        }

        /// <summary>
        /// Solves a system of ODEs and returns the full trajectory.
        /// Uses the classical 4th order Runge-Kutta method.
        /// </summary>
        /// <param name="func">Right-hand side function f(t, y).</param>
        /// <param name="tMin">Initial time.</param>
        /// <param name="tMax">Final time.</param>
        /// <param name="stepSize">Time step size.</param>
        /// <param name="y0">Initial state vector y(tMin).</param>
        /// <returns>Full trajectory as a list of (t, y) pairs.</returns>
        public static List<(double t, VectorN y)> RungeKuttaTrajectory(
            this Func<(double t, VectorN y), VectorN> func,
            double tMin, double tMax, double stepSize, VectorN y0)
        {
            var trajectory = new List<(double t, VectorN y)>();
            var y = new VectorN(y0.Values);
            trajectory.Add((tMin, new VectorN(y.Values)));

            for (var t = tMin; t < tMax; t += stepSize)
            {
                var h = Math.Min(stepSize, tMax - t);

                var k1 = func((t, y));
                var k2 = func((t + 0.5 * h, y + 0.5 * h * k1));
                var k3 = func((t + 0.5 * h, y + 0.5 * h * k2));
                var k4 = func((t + h, y + h * k3));

                y = y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
                trajectory.Add((t + h, new VectorN(y.Values)));
            }

            return trajectory;
        }

        /// <summary>
        /// Convenience overload: solves a system of ODEs using double[] via VectorN RK4.
        /// </summary>
        public static double[] RungeKutta(
            this Func<(double t, double[] y), double[]> func,
            double tMin, double tMax, double stepSize, double[] y0)
        {
            Func<(double t, VectorN y), VectorN> wrapped = v =>
                new VectorN(func((v.t, v.y.Values)));

            return wrapped.RungeKutta(tMin, tMax, stepSize, new VectorN(y0)).Values;
        }

        /// <summary>
        /// Convenience overload: solves a system of ODEs using double[] via VectorN Euler.
        /// </summary>
        public static double[] EulerMethod(
            this Func<(double t, double[] y), double[]> func,
            double tMin, double tMax, double stepSize, double[] y0)
        {
            Func<(double t, VectorN y), VectorN> wrapped = v =>
                new VectorN(func((v.t, v.y.Values)));

            return wrapped.EulerMethod(tMin, tMax, stepSize, new VectorN(y0)).Values;
        }

        /// <summary>
        /// Convenience overload: solves a system of ODEs and returns trajectory using double[] via VectorN RK4.
        /// </summary>
        public static List<(double t, double[] y)> RungeKuttaTrajectory(
            this Func<(double t, double[] y), double[]> func,
            double tMin, double tMax, double stepSize, double[] y0)
        {
            Func<(double t, VectorN y), VectorN> wrapped = v =>
                new VectorN(func((v.t, v.y.Values)));

            return wrapped.RungeKuttaTrajectory(tMin, tMax, stepSize, new VectorN(y0))
                .Select(p => (p.t, p.y.Values))
                .ToList();
        }

        #endregion

        #region Semi-Implicit Euler (Symplectic)

        /// <summary>
        /// Solves a second-order ODE system using semi-implicit (symplectic) Euler.
        /// State y = [positions, velocities] (must have even length).
        /// func returns the full derivative [velocities, accelerations].
        /// Update order: v += a·dt, then x += v_new·dt (uses updated velocity for stability).
        /// </summary>
        public static VectorN SemiImplicitEuler(
            this Func<(double t, VectorN y), VectorN> func,
            double tMin, double tMax, double stepSize, VectorN y0)
        {
            var y = new VectorN(y0.Values);
            int half = y.Length / 2;

            for (var t = tMin; t < tMax; t += stepSize)
            {
                var h = Math.Min(stepSize, tMax - t);
                var dy = func((t, y));

                for (int i = half; i < y.Length; i++)
                    y[i] += h * dy[i];

                for (int i = 0; i < half; i++)
                    y[i] += h * y[half + i];
            }

            return y;
        }

        /// <summary>
        /// Semi-implicit Euler with full trajectory output.
        /// </summary>
        public static List<(double t, VectorN y)> SemiImplicitEulerTrajectory(
            this Func<(double t, VectorN y), VectorN> func,
            double tMin, double tMax, double stepSize, VectorN y0)
        {
            var trajectory = new List<(double t, VectorN y)>();
            var y = new VectorN(y0.Values);
            int half = y.Length / 2;
            trajectory.Add((tMin, new VectorN(y.Values)));

            for (var t = tMin; t < tMax; t += stepSize)
            {
                var h = Math.Min(stepSize, tMax - t);
                var dy = func((t, y));

                for (int i = half; i < y.Length; i++)
                    y[i] += h * dy[i];

                for (int i = 0; i < half; i++)
                    y[i] += h * y[half + i];

                trajectory.Add((t + h, new VectorN(y.Values)));
            }

            return trajectory;
        }

        /// <summary>
        /// Convenience overload: semi-implicit Euler with double[].
        /// </summary>
        public static double[] SemiImplicitEuler(
            this Func<(double t, double[] y), double[]> func,
            double tMin, double tMax, double stepSize, double[] y0)
        {
            Func<(double t, VectorN y), VectorN> wrapped = v =>
                new VectorN(func((v.t, v.y.Values)));

            return wrapped.SemiImplicitEuler(tMin, tMax, stepSize, new VectorN(y0)).Values;
        }

        /// <summary>
        /// Convenience overload: semi-implicit Euler trajectory with double[].
        /// </summary>
        public static List<(double t, double[] y)> SemiImplicitEulerTrajectory(
            this Func<(double t, double[] y), double[]> func,
            double tMin, double tMax, double stepSize, double[] y0)
        {
            Func<(double t, VectorN y), VectorN> wrapped = v =>
                new VectorN(func((v.t, v.y.Values)));

            return wrapped.SemiImplicitEulerTrajectory(tMin, tMax, stepSize, new VectorN(y0))
                .Select(p => (p.t, p.y.Values))
                .ToList();
        }

        #endregion

        #region Velocity Verlet

        /// <summary>
        /// Solves a second-order ODE system using velocity Verlet integration.
        /// State y = [positions, velocities] (must have even length).
        /// accelerationFunc takes (t, state) and returns accelerations only (half the state length).
        /// Provides O(dt²) accuracy and excellent energy conservation.
        /// </summary>
        public static VectorN VelocityVerlet(
            this Func<(double t, VectorN y), VectorN> accelerationFunc,
            double tMin, double tMax, double stepSize, VectorN y0)
        {
            var y = new VectorN(y0.Values);
            int half = y.Length / 2;
            var a = accelerationFunc((tMin, y));

            for (var t = tMin; t < tMax; t += stepSize)
            {
                var h = Math.Min(stepSize, tMax - t);

                for (int i = 0; i < half; i++)
                    y[i] += h * y[half + i] + 0.5 * h * h * a[i];

                var aNew = accelerationFunc((t + h, y));

                for (int i = 0; i < half; i++)
                    y[half + i] += 0.5 * h * (a[i] + aNew[i]);

                a = aNew;
            }

            return y;
        }

        /// <summary>
        /// Velocity Verlet with full trajectory output.
        /// </summary>
        public static List<(double t, VectorN y)> VelocityVerletTrajectory(
            this Func<(double t, VectorN y), VectorN> accelerationFunc,
            double tMin, double tMax, double stepSize, VectorN y0)
        {
            var trajectory = new List<(double t, VectorN y)>();
            var y = new VectorN(y0.Values);
            int half = y.Length / 2;
            trajectory.Add((tMin, new VectorN(y.Values)));
            var a = accelerationFunc((tMin, y));

            for (var t = tMin; t < tMax; t += stepSize)
            {
                var h = Math.Min(stepSize, tMax - t);

                for (int i = 0; i < half; i++)
                    y[i] += h * y[half + i] + 0.5 * h * h * a[i];

                var aNew = accelerationFunc((t + h, y));

                for (int i = 0; i < half; i++)
                    y[half + i] += 0.5 * h * (a[i] + aNew[i]);

                a = aNew;
                trajectory.Add((t + h, new VectorN(y.Values)));
            }

            return trajectory;
        }

        /// <summary>
        /// Convenience overload: velocity Verlet with double[].
        /// </summary>
        public static double[] VelocityVerlet(
            this Func<(double t, double[] y), double[]> accelerationFunc,
            double tMin, double tMax, double stepSize, double[] y0)
        {
            Func<(double t, VectorN y), VectorN> wrapped = v =>
                new VectorN(accelerationFunc((v.t, v.y.Values)));

            return wrapped.VelocityVerlet(tMin, tMax, stepSize, new VectorN(y0)).Values;
        }

        /// <summary>
        /// Convenience overload: velocity Verlet trajectory with double[].
        /// </summary>
        public static List<(double t, double[] y)> VelocityVerletTrajectory(
            this Func<(double t, double[] y), double[]> accelerationFunc,
            double tMin, double tMax, double stepSize, double[] y0)
        {
            Func<(double t, VectorN y), VectorN> wrapped = v =>
                new VectorN(accelerationFunc((v.t, v.y.Values)));

            return wrapped.VelocityVerletTrajectory(tMin, tMax, stepSize, new VectorN(y0))
                .Select(p => (p.t, p.y.Values))
                .ToList();
        }

        #endregion

        /// <summary>
        /// Solves a linear system A x = b by computing x = A^{-1} b.
        /// </summary>
        /// <param name="matrix">Coefficient matrix A.</param>
        /// <param name="vector">Right-hand side vector b.</param>
        /// <returns>Solution vector x.</returns>
        public static Vector LinearSystemSolver(this Matrix matrix, Vector vector)
        {
            var values = matrix.Inverse() * vector;

            return values;
        }

        /// <summary>
        /// Solves a linear system A x = b by computing x = A^{-1} b.
        /// </summary>
        /// <param name="matrix">Coefficient matrix A.</param>
        /// <param name="vector">Right-hand side vector b.</param>
        /// <returns>Solution vector x.</returns>
        public static VectorN LinearSystemSolver(this Matrix matrix, VectorN vector)
        {
            var values = matrix.Inverse() * vector;

            return values;
        }

        /// <summary>
        /// Solves a linear system A x = b by computing x = A^{-1} b.
        /// </summary>
        /// <param name="matrix">Coefficient matrix A.</param>
        /// <param name="vector">Right-hand side vector b.</param>
        /// <returns>Solution vector x as a list.</returns>
        public static List<double> LinearSystemSolver(this Matrix matrix, List<double> vector)
        {
            var values = matrix.Inverse() * vector;

            return values;
        }

        /// <summary>
        /// Solves a 3x3 linear system using Gaussian elimination.
        /// </summary>
        /// <param name="matrix">Coefficient matrix.</param>
        /// <param name="vector">Right-hand side vector.</param>
        /// <returns>Solution as a <see cref="Vector"/>.</returns>
        public static Vector GaussElimination(this Matrix matrix, Vector vector)
        {
            var values = matrix.GaussElimination(new List<double>() { vector.x, vector.y, vector.z });

            return new Vector(values[0], values[1], values[2]);
        }

        /// <summary>
        /// Solves a linear system using Gaussian elimination with partial pivoting.
        /// </summary>
        /// <param name="matrix">Coefficient matrix (will be modified in-place).</param>
        /// <param name="vector">Right-hand side vector (will be modified in-place).</param>
        /// <returns>Solution vector as a list.</returns>
        public static List<double> GaussElimination(this Matrix matrix, List<double> vector)
        {
            var n = vector.Count;

            for (var p = 0; p < n; p++)
            {
                var max = p;
                for (var i = p + 1; i < n; i++)
                {
                    if (Math.Abs(matrix.values[i, p]) > Math.Abs(matrix.values[max, p]))
                    {
                        max = i;
                    }
                }

                for (var c = 0; c < matrix.columnLength; c++)
                {
                    var temp = matrix.values[p, c];

                    matrix.values[p, c] = matrix.values[max, c];

                    matrix.values[max, c] = temp;
                }

                var t = vector[p];
                vector[p] = vector[max];
                vector[max] = t;

                for (var i = p + 1; i < n; i++)
                {
                    var alpha = matrix.values[i, p] / matrix.values[p, p];
                    vector[i] -= alpha * vector[p];
                    for (var j = p; j < n; j++)
                    {
                        matrix.values[i, j] -= alpha * matrix.values[p, j];
                    }
                }
            }

            var x = new double[n];
            for (var i = n - 1; i >= 0; i--)
            {
                var sum = 0.0;
                for (var j = i + 1; j < n; j++)
                {
                    sum += matrix.values[i, j] * x[j];
                }
                if (matrix.values[i, i] != 0)
                {
                    if ((vector[i] - sum) != 0)
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

        /// <summary>
        /// Computes an (approximate) eigenvector corresponding to a provided eigenvalue.
        /// </summary>
        /// <param name="matrix">The matrix.</param>
        /// <param name="egienValue">The eigenvalue to use.</param>
        /// <returns>An approximate eigenvector (scaled/rounded) as a list.</returns>
        public static List<double> EigenVector(this Matrix matrix, double egienValue)
        {
            for (var i = 0; i < matrix.rowLength; i++)
            {
                for (var j = 0; j < matrix.columnLength; j++)
                {
                    if (i == j)
                    {
                        matrix.values[i, j] -= egienValue - 0.1;
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

        /// <summary>
        /// Computes an approximation of the dominant eigenvector using power iteration.
        /// </summary>
        /// <param name="matrix">The matrix.</param>
        /// <returns>An approximate dominant eigenvector (scaled/rounded) as a list.</returns>
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

        /// <summary>
        /// Builds solution functions for a linear ODE system using eigen-decomposition, given a single initial value repeated for each dimension.
        /// </summary>
        /// <param name="matrix">System matrix.</param>
        /// <param name="tZero">Initial value used for each component.</param>
        /// <returns>A list of solution functions (one per component).</returns>
        public static List<Func<double, double>> OdeSolver(this Matrix matrix, double tZero)
        {
            var list = Enumerable.Range(1, matrix.rowLength).Select(x => tZero).ToList();

            return matrix.OdeSolver(list);
        }

        /// <summary>
        /// Builds solution functions for a linear ODE system using eigenvalues/eigenvectors and initial conditions.
        /// </summary>
        /// <param name="matrix">System matrix.</param>
        /// <param name="tZeros">Initial values for each component.</param>
        /// <returns>A list of solution functions (one per component).</returns>
        public static List<Func<double, double>> OdeSolver(this Matrix matrix, List<double> tZeros)
        {
            var eigenValues = matrix.EigenValues();
            var eigenVectors = new Dictionary<double, List<double>>();

            eigenValues.Reverse();

            foreach (var value in eigenValues)
            {
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
                   (constants.Count > 2 ? Math.Exp(eigenValues[2] * t) * eigenVectors[eigenValues[2]][index] * constants[2] : 0.0);
        }

        /// <summary>
        /// Computes a sequence of eigenvalues using repeated dominant-eigenvector estimation and deflation.
        /// </summary>
        /// <param name="matrix">The matrix.</param>
        /// <returns>A list of eigenvalue estimates.</returns>
        public static List<double> EigenValues(this Matrix matrix)
        {
            var results = new List<double>();
            var dominantEigenvector = matrix.DominantEigenVector();

            results.Add((matrix * dominantEigenvector).Dot(dominantEigenvector) / dominantEigenvector.Dot(dominantEigenvector));

            for (var s = 0; s < matrix.rowLength - 1; s++)
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

        /// <summary>
        /// Computes the dominant eigenvalue using the Rayleigh quotient of an approximate dominant eigenvector.
        /// </summary>
        /// <param name="matrix">The matrix.</param>
        /// <returns>The dominant eigenvalue estimate.</returns>
        public static double DominantEigenValue(this Matrix matrix)
        {
            var result = matrix.DominantEigenVector();

            var eigenvalue = (matrix * result).Dot(result) / result.Dot(result);

            return eigenvalue;
        }
    }
}
