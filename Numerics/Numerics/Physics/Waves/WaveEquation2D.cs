using CSharpNumerics.Numerics.FiniteDifference;
using CSharpNumerics.Numerics.FiniteDifference.TimeStepping;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Data;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Physics.Waves
{
    /// <summary>
    /// Solves the 2D wave equation ∂²u/∂t² = c² ∇²u on a <see cref="Grid2D"/> using Method of Lines.
    /// <para>
    /// The 5-point Laplacian stencil from <see cref="GridOperators.Laplacian2D"/> is used for
    /// spatial discretisation. The resulting ODE system is advanced with an <see cref="ITimeStepper"/>.
    /// </para>
    /// <para>
    /// State vector layout: y = [u_flat | v_flat], each of length Nx × Ny.
    /// </para>
    /// </summary>
    public class WaveEquation2D : IWaveField
    {
        #region Fields

        private readonly Grid2D _grid;
        private readonly double _c;
        private readonly double _c2;
        private readonly BoundaryCondition _bc;
        private readonly ITimeStepper _stepper;
        private readonly int _fieldLen;

        private VectorN _state;
        private VectorN _initialState;
        private double _t;

        #endregion

        #region Constructor

        /// <summary>
        /// Creates a 2D wave equation solver on the given grid.
        /// </summary>
        /// <param name="grid">Spatial grid.</param>
        /// <param name="c">Wave speed (m/s).</param>
        /// <param name="boundaryType">Boundary condition type.</param>
        /// <param name="stepper">Time stepper (defaults to Velocity Verlet).</param>
        public WaveEquation2D(Grid2D grid, double c,
            BoundaryType boundaryType = BoundaryType.Fixed,
            ITimeStepper stepper = null)
        {
            if (grid == null) throw new ArgumentNullException(nameof(grid));
            if (c <= 0) throw new ArgumentException("Wave speed must be positive.", nameof(c));

            _grid = grid;
            _c = c;
            _c2 = c * c;
            _bc = WaveEquation1D.MapBoundaryPublic(boundaryType);
            _stepper = stepper ?? new VelocityVerletStepper();
            _fieldLen = grid.Length;

            _state = new VectorN(new double[2 * _fieldLen]);
            _initialState = new VectorN(new double[2 * _fieldLen]);
            _t = 0.0;
        }

        #endregion

        #region Properties

        /// <summary>The spatial grid.</summary>
        public Grid2D Grid => _grid;

        /// <summary>Wave speed (m/s).</summary>
        public double WaveSpeed => _c;

        /// <inheritdoc />
        public double Time => _t;

        /// <inheritdoc />
        public double TotalEnergy
        {
            get
            {
                var u = ExtractField(_state, 0);
                var v = ExtractField(_state, _fieldLen);
                var (dux, duy) = GridOperators.Gradient2D(u, _grid, _bc);

                double ke = 0, pe = 0;
                for (int i = 0; i < _fieldLen; i++)
                {
                    ke += v[i] * v[i];
                    pe += dux[i] * dux[i] + duy[i] * duy[i];
                }

                double cellArea = _grid.Dx * _grid.Dy;
                return 0.5 * cellArea * (ke + _c2 * pe);
            }
        }

        #endregion

        #region Initial Conditions

        /// <summary>
        /// Sets the initial displacement and velocity fields.
        /// </summary>
        /// <param name="u0">Initial displacement u(x, y, 0). If null, zero.</param>
        /// <param name="v0">Initial velocity ∂u/∂t(x, y, 0). If null, zero.</param>
        public void SetInitialCondition(
            Func<double, double, double> u0 = null,
            Func<double, double, double> v0 = null)
        {
            var values = new double[2 * _fieldLen];

            for (int iy = 0; iy < _grid.Ny; iy++)
            {
                for (int ix = 0; ix < _grid.Nx; ix++)
                {
                    int idx = _grid.Index(ix, iy);
                    double x = ix * _grid.Dx;
                    double y = iy * _grid.Dy;
                    values[idx] = u0?.Invoke(x, y) ?? 0.0;
                    values[_fieldLen + idx] = v0?.Invoke(x, y) ?? 0.0;
                }
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

        /// <inheritdoc />
        public List<Serie> Snapshot()
        {
            var result = new List<Serie>(_fieldLen);
            for (int i = 0; i < _fieldLen; i++)
                result.Add(new Serie { Index = i, Value = _state[i] });
            return result;
        }

        /// <summary>
        /// Returns the current displacement field as a 2D array [ix, iy].
        /// </summary>
        public double[,] SnapshotArray()
        {
            var u = ExtractField(_state, 0);
            return _grid.ToArray(u);
        }

        /// <summary>
        /// Returns the energy density at each grid cell.
        /// </summary>
        public double[,] EnergyDensityArray()
        {
            var u = ExtractField(_state, 0);
            var v = ExtractField(_state, _fieldLen);
            var (dux, duy) = GridOperators.Gradient2D(u, _grid, _bc);

            var e = new double[_fieldLen];
            for (int i = 0; i < _fieldLen; i++)
                e[i] = 0.5 * v[i] * v[i] + 0.5 * _c2 * (dux[i] * dux[i] + duy[i] * duy[i]);

            return _grid.ToArray(new VectorN(e));
        }

        #endregion

        #region Internals

        private VectorN Rhs(double t, VectorN y)
        {
            var u = ExtractField(y, 0);
            var v = ExtractField(y, _fieldLen);
            var lap = GridOperators.Laplacian2D(u, _grid, _bc);

            var dydt = new double[2 * _fieldLen];
            for (int i = 0; i < _fieldLen; i++)
            {
                dydt[i] = v[i];
                dydt[_fieldLen + i] = _c2 * lap[i];
            }

            return new VectorN(dydt);
        }

        private VectorN ExtractField(VectorN y, int offset)
        {
            var vals = new double[_fieldLen];
            Array.Copy(y.Values, offset, vals, 0, _fieldLen);
            return new VectorN(vals);
        }

        #endregion
    }
}
