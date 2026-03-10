using CSharpNumerics.Numerics;
using CSharpNumerics.Numerics.FiniteDifference;
using CSharpNumerics.Numerics.FiniteDifference.TimeStepping;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Data;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Physics.Waves
{
    /// <summary>
    /// Solves the damped/driven 1D wave equation:
    /// ∂²u/∂t² + 2α ∂u/∂t = c² ∂²u/∂x² + S(x, t)
    /// <para>
    /// When α = 0 and S = 0 this reduces to the standard wave equation.
    /// Uses the same MOL approach as <see cref="WaveEquation1D"/> with an extra damping
    /// term and an optional source function.
    /// </para>
    /// </summary>
    public class DampedDrivenWaveEquation1D : IWaveField
    {
        #region Fields

        private readonly int _n;
        private readonly double _dx;
        private readonly double _c2;
        private readonly double _alpha;
        private readonly BoundaryCondition _bc;
        private readonly ITimeStepper _stepper;
        private Func<double, double, double> _source;

        private VectorN _state;
        private VectorN _initialState;
        private double _t;

        #endregion

        #region Constructor

        /// <summary>
        /// Creates a damped/driven 1D wave equation solver.
        /// </summary>
        /// <param name="n">Number of spatial grid points.</param>
        /// <param name="dx">Grid spacing (m).</param>
        /// <param name="c">Wave speed (m/s).</param>
        /// <param name="alpha">Damping coefficient α (1/s). 0 = undamped.</param>
        /// <param name="source">
        /// Source term S(x, t). If null, there is no driving force.
        /// </param>
        /// <param name="boundaryType">Boundary condition type.</param>
        /// <param name="stepper">Time stepper (defaults to RK4 since the system is dissipative).</param>
        public DampedDrivenWaveEquation1D(int n, double dx, double c,
            double alpha = 0,
            Func<double, double, double> source = null,
            BoundaryType boundaryType = BoundaryType.Fixed,
            ITimeStepper stepper = null)
        {
            if (n < 3) throw new ArgumentException("Need at least 3 grid points.", nameof(n));
            if (dx <= 0) throw new ArgumentException("Grid spacing must be positive.", nameof(dx));
            if (c <= 0) throw new ArgumentException("Wave speed must be positive.", nameof(c));
            if (alpha < 0) throw new ArgumentException("Damping must be non-negative.", nameof(alpha));

            _n = n;
            _dx = dx;
            _c2 = c * c;
            _alpha = alpha;
            _source = source;
            _bc = WaveEquation1D.MapBoundaryPublic(boundaryType);
            _stepper = stepper ?? new RK4Stepper();

            _state = new VectorN(new double[2 * _n]);
            _initialState = new VectorN(new double[2 * _n]);
            _t = 0.0;
        }

        #endregion

        #region Properties

        /// <summary>Number of spatial grid points.</summary>
        public int N => _n;

        /// <summary>Grid spacing (m).</summary>
        public double Dx => _dx;

        /// <summary>Damping coefficient α (1/s).</summary>
        public double Alpha => _alpha;

        /// <inheritdoc />
        public double Time => _t;

        /// <inheritdoc />
        public double TotalEnergy
        {
            get
            {
                var u = ExtractU(_state);
                var v = ExtractV(_state);
                var grad = GridOperators.Gradient1D(u, _dx, _bc);

                double ke = 0, pe = 0;
                for (int i = 0; i < _n; i++)
                {
                    ke += v[i] * v[i];
                    pe += grad[i] * grad[i];
                }

                return 0.5 * _dx * (ke + _c2 * pe);
            }
        }

        #endregion

        #region Initial Conditions

        /// <summary>
        /// Sets the initial displacement and velocity fields.
        /// </summary>
        public void SetInitialCondition(Func<double, double> u0 = null, Func<double, double> v0 = null)
        {
            var values = new double[2 * _n];
            for (int i = 0; i < _n; i++)
            {
                double x = i * _dx;
                values[i] = u0?.Invoke(x) ?? 0.0;
                values[_n + i] = v0?.Invoke(x) ?? 0.0;
            }

            _state = new VectorN(values);
            _initialState = new VectorN((double[])values.Clone());
            _t = 0.0;
        }

        /// <summary>
        /// Sets or replaces the source/driving function S(x, t).
        /// </summary>
        public void SetSource(Func<double, double, double> source)
        {
            _source = source;
        }

        #endregion

        #region Simulation

        /// <inheritdoc />
        public void Step(double dt)
        {
            if (dt <= 0) throw new ArgumentException("Time step must be positive.", nameof(dt));
            var result = _stepper.Solve(Rhs, _t, _t + dt, _state, dt);
            _state = result.Y;
            _t = result.T;
        }

        /// <inheritdoc />
        public void Reset()
        {
            _state = new VectorN((double[])_initialState.Values.Clone());
            _t = 0.0;
        }

        /// <summary>
        /// Runs the simulation to <paramref name="tEnd"/>.
        /// </summary>
        public TimeStepResult Simulate(double tEnd, double dt, bool recordTrajectory = false)
        {
            if (tEnd <= _t) throw new ArgumentException("End time must be greater than current time.");
            if (dt <= 0) throw new ArgumentException("Time step must be positive.", nameof(dt));

            var result = _stepper.Solve(Rhs, _t, tEnd, _state, dt, recordTrajectory);
            _state = result.Y;
            _t = result.T;
            return result;
        }

        #endregion

        #region Snapshots

        /// <inheritdoc />
        public List<Serie> Snapshot()
        {
            var result = new List<Serie>(_n);
            for (int i = 0; i < _n; i++)
                result.Add(new Serie { Index = i * _dx, Value = _state[i] });
            return result;
        }

        #endregion

        #region Internals

        private VectorN Rhs(double t, VectorN y)
        {
            var u = ExtractU(y);
            var v = ExtractV(y);
            var lap = GridOperators.Laplacian1D(u, _dx, _bc);

            var dydt = new double[2 * _n];
            for (int i = 0; i < _n; i++)
            {
                double x = i * _dx;
                double src = _source?.Invoke(x, t) ?? 0.0;
                dydt[i] = v[i];
                dydt[_n + i] = _c2 * lap[i] - 2.0 * _alpha * v[i] + src;
            }

            return new VectorN(dydt);
        }

        private VectorN ExtractU(VectorN y)
        {
            var vals = new double[_n];
            Array.Copy(y.Values, 0, vals, 0, _n);
            return new VectorN(vals);
        }

        private VectorN ExtractV(VectorN y)
        {
            var vals = new double[_n];
            Array.Copy(y.Values, _n, vals, 0, _n);
            return new VectorN(vals);
        }

        #endregion
    }
}
