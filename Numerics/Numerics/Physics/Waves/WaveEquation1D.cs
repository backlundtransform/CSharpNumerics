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
    /// Solves the 1D wave equation ∂²u/∂t² = c² ∂²u/∂x² using Method of Lines.
    /// <para>
    /// The spatial domain is discretised into N equally spaced points and the Laplacian
    /// is computed via <see cref="GridOperators.Laplacian1D"/>. The resulting ODE system
    /// is advanced with a configurable <see cref="ITimeStepper"/> (default: Velocity Verlet).
    /// </para>
    /// <para>
    /// State vector layout: y = [u₀ … u_{N-1} | v₀ … v_{N-1}]
    /// where u is displacement and v = ∂u/∂t.
    /// </para>
    /// </summary>
    public class WaveEquation1D : IWaveField
    {
        #region Fields

        private readonly int _n;
        private readonly double _dx;
        private readonly double _c;
        private readonly double _c2;
        private readonly BoundaryCondition _bc;
        private readonly ITimeStepper _stepper;

        private VectorN _state;       // [u | v]
        private VectorN _initialState;
        private double _t;

        #endregion

        #region Constructor

        /// <summary>
        /// Creates a 1D wave equation solver.
        /// </summary>
        /// <param name="n">Number of spatial grid points.</param>
        /// <param name="dx">Grid spacing (m).</param>
        /// <param name="c">Wave speed (m/s).</param>
        /// <param name="boundaryType">Boundary condition type.</param>
        /// <param name="stepper">
        /// Time-stepping method. If null, defaults to <see cref="VelocityVerletStepper"/>.
        /// </param>
        public WaveEquation1D(int n, double dx, double c,
            BoundaryType boundaryType = BoundaryType.Fixed,
            ITimeStepper stepper = null)
        {
            if (n < 3) throw new ArgumentException("Need at least 3 grid points.", nameof(n));
            if (dx <= 0) throw new ArgumentException("Grid spacing must be positive.", nameof(dx));
            if (c <= 0) throw new ArgumentException("Wave speed must be positive.", nameof(c));

            _n = n;
            _dx = dx;
            _c = c;
            _c2 = c * c;
            _bc = MapBoundary(boundaryType);
            _stepper = stepper ?? new VelocityVerletStepper();

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

        /// <summary>Wave speed (m/s).</summary>
        public double WaveSpeed => _c;

        /// <summary>Domain length L = (N-1) * dx.</summary>
        public double Length => (_n - 1) * _dx;

        /// <inheritdoc />
        public double Time => _t;

        /// <summary>
        /// CFL number for a given time step: C = c * dt / dx.
        /// Must be ≤ 1 for stability.
        /// </summary>
        public double CFL(double dt) => _c * dt / _dx;

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
        /// <param name="u0">Initial displacement u(x, 0). If null, zero.</param>
        /// <param name="v0">Initial velocity ∂u/∂t(x, 0). If null, zero.</param>
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
        /// Runs the simulation from the current time to <paramref name="tEnd"/>.
        /// </summary>
        /// <param name="tEnd">End time (s).</param>
        /// <param name="dt">Time step (s).</param>
        /// <param name="recordTrajectory">If true, records the full space-time evolution.</param>
        /// <returns>Integration result with optional trajectory.</returns>
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

        #region Snapshots & Output

        /// <summary>
        /// Returns the displacement at grid point <paramref name="index"/>.
        /// </summary>
        public double Displacement(int index) => _state[index];

        /// <summary>
        /// Returns the velocity at grid point <paramref name="index"/>.
        /// </summary>
        public double Velocity(int index) => _state[_n + index];

        /// <inheritdoc />
        public List<Serie> Snapshot()
        {
            var result = new List<Serie>(_n);
            for (int i = 0; i < _n; i++)
                result.Add(new Serie { Index = i * _dx, Value = _state[i] });
            return result;
        }

        /// <summary>
        /// Returns the full space-time field u(x, t) by running from t = 0 to <paramref name="tEnd"/>.
        /// Rows = spatial points, columns = time steps.
        /// The field is reset to initial conditions after the call.
        /// </summary>
        /// <param name="tEnd">End time (s).</param>
        /// <param name="dt">Time step (s).</param>
        public Matrix SpaceTimeField(double tEnd, double dt)
        {
            if (tEnd <= 0) throw new ArgumentException("End time must be positive.", nameof(tEnd));
            if (dt <= 0) throw new ArgumentException("Time step must be positive.", nameof(dt));

            Reset();
            int nSteps = (int)Math.Ceiling(tEnd / dt) + 1;
            var data = new double[_n, nSteps];

            // Record initial state
            for (int i = 0; i < _n; i++)
                data[i, 0] = _state[i];

            int col = 1;
            while (_t < tEnd && col < nSteps)
            {
                Step(dt);
                for (int i = 0; i < _n; i++)
                    data[i, col] = _state[i];
                col++;
            }

            Reset();
            return new Matrix(data);
        }

        /// <summary>
        /// Returns the energy density at each grid point: ½v² + ½c²(∂u/∂x)².
        /// </summary>
        public List<Serie> EnergyDensity()
        {
            var u = ExtractU(_state);
            var v = ExtractV(_state);
            var grad = GridOperators.Gradient1D(u, _dx, _bc);

            var result = new List<Serie>(_n);
            for (int i = 0; i < _n; i++)
            {
                double e = 0.5 * v[i] * v[i] + 0.5 * _c2 * grad[i] * grad[i];
                result.Add(new Serie { Index = i * _dx, Value = e });
            }

            return result;
        }

        /// <summary>
        /// Returns the frequency content at a given spatial point by recording
        /// the displacement over time and computing the FFT.
        /// The field is reset after the call.
        /// </summary>
        /// <param name="spatialIndex">Grid-point index to sample.</param>
        /// <param name="tEnd">Simulation duration (s).</param>
        /// <param name="dt">Time step (s).</param>
        public List<Serie> FrequencyContent(int spatialIndex, double tEnd, double dt)
        {
            if (spatialIndex < 0 || spatialIndex >= _n)
                throw new ArgumentOutOfRangeException(nameof(spatialIndex));
            if (tEnd <= 0) throw new ArgumentException("End time must be positive.", nameof(tEnd));
            if (dt <= 0) throw new ArgumentException("Time step must be positive.", nameof(dt));

            Reset();
            var samples = new List<double> { _state[spatialIndex] };

            while (_t < tEnd)
            {
                Step(dt);
                samples.Add(_state[spatialIndex]);
            }

            Reset();

            // Pad to next power of two
            int nFft = 1;
            while (nFft < samples.Count) nFft <<= 1;
            while (samples.Count < nFft) samples.Add(0.0);

            double samplingFrequency = 1.0 / dt;
            return samples.FastFourierTransform().ToFrequencyResolution(samplingFrequency);
        }

        #endregion

        #region Analytic Reference

        /// <summary>
        /// Returns the analytic mode shapes for a string with fixed ends: sin(nπx/L).
        /// Useful for verifying numerical solutions.
        /// </summary>
        /// <param name="modeNumber">Mode number (1, 2, 3, …).</param>
        public List<Serie> StandingWaveMode(int modeNumber)
        {
            if (modeNumber < 1)
                throw new ArgumentException("Mode number must be ≥ 1.", nameof(modeNumber));

            double l = Length;
            var result = new List<Serie>(_n);
            for (int i = 0; i < _n; i++)
            {
                double x = i * _dx;
                result.Add(new Serie { Index = x, Value = Math.Sin(modeNumber * Math.PI * x / l) });
            }

            return result;
        }

        /// <summary>
        /// Analytic frequency of the n-th standing wave mode for fixed/fixed boundaries:
        /// f_n = n * c / (2L).
        /// </summary>
        /// <param name="modeNumber">Mode number (1, 2, 3, …).</param>
        public double StandingWaveFrequency(int modeNumber)
        {
            if (modeNumber < 1)
                throw new ArgumentException("Mode number must be ≥ 1.", nameof(modeNumber));
            return modeNumber * _c / (2.0 * Length);
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
                dydt[i] = v[i];             // du/dt = v
                dydt[_n + i] = _c2 * lap[i]; // dv/dt = c² ∇²u
            }

            return new VectorN(dydt);
        }

        private VectorN ExtractU(VectorN y)
        {
            var u = new double[_n];
            Array.Copy(y.Values, 0, u, 0, _n);
            return new VectorN(u);
        }

        private VectorN ExtractV(VectorN y)
        {
            var v = new double[_n];
            Array.Copy(y.Values, _n, v, 0, _n);
            return new VectorN(v);
        }

        private static BoundaryCondition MapBoundary(BoundaryType bt) => MapBoundaryPublic(bt);

        /// <summary>
        /// Maps a <see cref="BoundaryType"/> to the corresponding <see cref="BoundaryCondition"/>
        /// used by the finite-difference operators.
        /// </summary>
        internal static BoundaryCondition MapBoundaryPublic(BoundaryType bt)
        {
            switch (bt)
            {
                case BoundaryType.Fixed: return BoundaryCondition.Dirichlet;
                case BoundaryType.Free: return BoundaryCondition.Neumann;
                case BoundaryType.Periodic: return BoundaryCondition.Periodic;
                case BoundaryType.Absorbing: return BoundaryCondition.Neumann; // simplified
                default: return BoundaryCondition.Dirichlet;
            }
        }

        #endregion
    }
}
