using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Simulation;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Environmental.Enums;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Engines.GIS.RL
{
    /// <summary>
    /// Defines a mitigation action the agent can take each time step.
    /// </summary>
    public enum MitigationAction
    {
        /// <summary>Do nothing — let the plume evolve naturally.</summary>
        None = 0,
        /// <summary>Deploy absorbent barrier in zone 0 (north).</summary>
        BarrierNorth = 1,
        /// <summary>Deploy absorbent barrier in zone 1 (east).</summary>
        BarrierEast = 2,
        /// <summary>Deploy absorbent barrier in zone 2 (south).</summary>
        BarrierSouth = 3,
        /// <summary>Deploy absorbent barrier in zone 3 (west).</summary>
        BarrierWest = 4,
        /// <summary>Reduce emission rate by activating a filter.</summary>
        ActivateFilter = 5
    }

    /// <summary>
    /// RL environment that wraps the GIS plume simulator.
    /// The agent steps through a transient plume scenario and takes mitigation
    /// actions at each time step to minimise population exposure.
    /// <para>
    /// State vector (8 elements):
    /// [maxConcentration, meanConcentration, exposedCellFraction,
    ///  windSpeed, windDirX, windDirY, emissionRate, normalizedTime]
    /// </para>
    /// <para>
    /// Actions: 6 discrete — see <see cref="MitigationAction"/>.
    /// </para>
    /// <para>
    /// Reward: −(exposed cell fraction) − actionCost, with bonus for keeping
    /// exposure below threshold at episode end.
    /// </para>
    /// </summary>
    public class PlumeEnvironment : IGISEnvironment
    {
        // ── Configuration ────────────────────────────────────────────

        private readonly double _baseEmissionRate;
        private readonly double _baseWindSpeed;
        private readonly Vector _baseWindDirection;
        private readonly double _stackHeight;
        private readonly Vector _sourcePosition;
        private readonly StabilityClass _stability;
        private readonly GeoGrid _grid;
        private readonly TimeFrame _timeFrame;

        /// <summary>Concentration threshold (kg/m³) above which a cell is "exposed".</summary>
        public double Threshold { get; set; } = 1e-6;

        /// <summary>Cost subtracted from reward per non-None action.</summary>
        public double ActionCost { get; set; } = 0.05;

        /// <summary>Fraction by which a barrier reduces concentration in its quadrant.</summary>
        public double BarrierEfficiency { get; set; } = 0.4;

        /// <summary>Fraction by which the filter reduces emission rate.</summary>
        public double FilterEfficiency { get; set; } = 0.5;

        // ── Episode state ────────────────────────────────────────────

        private int _currentStep;
        private double _currentEmissionRate;
        private double _currentWindSpeed;
        private Vector _currentWindDirection;
        private bool[] _activeBarriers;       // N E S W
        private bool _filterActive;
        private GridSnapshot _lastSnapshot;
        private Random _rng;

        // ── IEnvironment ─────────────────────────────────────────────

        /// <inheritdoc/>
        public int ObservationSize => 8;

        /// <inheritdoc/>
        public int ActionSize => 6;

        /// <inheritdoc/>
        public bool IsDiscrete => true;

        /// <summary>Total number of time steps per episode.</summary>
        public int MaxSteps => _timeFrame.Count;

        /// <inheritdoc/>
        public GeoGrid Grid => _grid;

        /// <inheritdoc/>
        public TimeFrame TimeFrame => _timeFrame;

        /// <inheritdoc/>
        public GridSnapshot LastSnapshot => _lastSnapshot;

        /// <summary>
        /// Creates a plume RL environment.
        /// </summary>
        /// <param name="emissionRate">Baseline emission rate Q (kg/s).</param>
        /// <param name="windSpeed">Baseline wind speed (m/s).</param>
        /// <param name="windDirection">Baseline wind direction vector.</param>
        /// <param name="stackHeight">Effective stack height (metres).</param>
        /// <param name="sourcePosition">Source position in world coordinates.</param>
        /// <param name="grid">Spatial discretization.</param>
        /// <param name="timeFrame">Temporal discretization.</param>
        /// <param name="stability">Atmospheric stability class.</param>
        public PlumeEnvironment(
            double emissionRate,
            double windSpeed,
            Vector windDirection,
            double stackHeight,
            Vector sourcePosition,
            GeoGrid grid,
            TimeFrame timeFrame,
            StabilityClass stability = StabilityClass.D)
        {
            if (emissionRate <= 0) throw new ArgumentException("Emission rate must be positive.", nameof(emissionRate));
            if (windSpeed <= 0) throw new ArgumentException("Wind speed must be positive.", nameof(windSpeed));

            _baseEmissionRate = emissionRate;
            _baseWindSpeed = windSpeed;
            _baseWindDirection = windDirection;
            _stackHeight = stackHeight;
            _sourcePosition = sourcePosition;
            _grid = grid ?? throw new ArgumentNullException(nameof(grid));
            _timeFrame = timeFrame ?? throw new ArgumentNullException(nameof(timeFrame));
            _stability = stability;
            _rng = new Random();
        }

        // ═══════════════════════════════════════════════════════════════
        //  IEnvironment implementation
        // ═══════════════════════════════════════════════════════════════

        /// <inheritdoc/>
        public (VectorN state, Dictionary<string, object> info) Reset(int? seed = null)
        {
            _rng = seed.HasValue ? new Random(seed.Value) : new Random();
            _currentStep = 0;
            _currentEmissionRate = _baseEmissionRate;
            _currentWindSpeed = _baseWindSpeed;
            _currentWindDirection = _baseWindDirection;
            _activeBarriers = new bool[4];
            _filterActive = false;

            _lastSnapshot = SimulateStep(0);

            return (BuildObservation(), new Dictionary<string, object>
            {
                ["snapshot"] = _lastSnapshot
            });
        }

        /// <inheritdoc/>
        public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(int action)
        {
            if (action < 0 || action >= ActionSize)
                throw new ArgumentOutOfRangeException(nameof(action));

            // Apply action
            ApplyAction((MitigationAction)action);

            // Advance time
            _currentStep++;
            bool done = _currentStep >= _timeFrame.Count;

            // Simulate
            if (!done)
                _lastSnapshot = SimulateStep(_currentStep);

            // Compute reward
            double exposedFraction = ComputeExposedFraction(_lastSnapshot);
            double reward = -exposedFraction;

            if (action != (int)MitigationAction.None)
                reward -= ActionCost;

            // Terminal bonus: if exposure at end is zero, reward
            if (done && exposedFraction < 0.01)
                reward += 1.0;

            var info = new Dictionary<string, object>
            {
                ["snapshot"] = _lastSnapshot,
                ["exposedFraction"] = exposedFraction,
                ["emissionRate"] = _currentEmissionRate,
                ["filterActive"] = _filterActive,
                ["barriers"] = (bool[])_activeBarriers.Clone()
            };

            return (BuildObservation(), reward, done, info);
        }

        /// <inheritdoc/>
        public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(VectorN action)
            => Step((int)action[0]);

        // ═══════════════════════════════════════════════════════════════
        //  Internals
        // ═══════════════════════════════════════════════════════════════

        private void ApplyAction(MitigationAction action)
        {
            switch (action)
            {
                case MitigationAction.BarrierNorth:
                    _activeBarriers[0] = true;
                    break;
                case MitigationAction.BarrierEast:
                    _activeBarriers[1] = true;
                    break;
                case MitigationAction.BarrierSouth:
                    _activeBarriers[2] = true;
                    break;
                case MitigationAction.BarrierWest:
                    _activeBarriers[3] = true;
                    break;
                case MitigationAction.ActivateFilter:
                    if (!_filterActive)
                    {
                        _filterActive = true;
                        _currentEmissionRate = _baseEmissionRate * (1.0 - FilterEfficiency);
                    }
                    break;
            }
        }

        private GridSnapshot SimulateStep(int stepIndex)
        {
            double time = _timeFrame.TimeAt(stepIndex);

            var sim = new PlumeSimulator(
                _currentEmissionRate,
                _currentWindSpeed,
                _currentWindDirection,
                _stackHeight,
                _sourcePosition,
                _stability,
                PlumeMode.Transient);
            sim.ReleaseSeconds = _timeFrame.StepSeconds;

            var snapshot = sim.RunSingle(_grid, time, stepIndex);

            // Apply barrier effects: reduce concentration in the relevant quadrant
            if (_activeBarriers.Any(b => b))
                snapshot = ApplyBarriers(snapshot);

            return snapshot;
        }

        private GridSnapshot ApplyBarriers(GridSnapshot original)
        {
            double[] values = original.GetValues();
            double centreX = (_grid.XMin + _grid.XMax) / 2.0;
            double centreY = (_grid.YMin + _grid.YMax) / 2.0;

            for (int i = 0; i < values.Length; i++)
            {
                if (values[i] <= 0) continue;

                Vector cellPos = _grid.CellCentre(i);
                bool inNorth = cellPos.y > centreY;
                bool inSouth = cellPos.y <= centreY;
                bool inEast = cellPos.x > centreX;
                bool inWest = cellPos.x <= centreX;

                double reduction = 1.0;
                if (_activeBarriers[0] && inNorth) reduction *= (1.0 - BarrierEfficiency);
                if (_activeBarriers[1] && inEast) reduction *= (1.0 - BarrierEfficiency);
                if (_activeBarriers[2] && inSouth) reduction *= (1.0 - BarrierEfficiency);
                if (_activeBarriers[3] && inWest) reduction *= (1.0 - BarrierEfficiency);

                values[i] *= reduction;
            }

            return new GridSnapshot(_grid, values, original.Time, original.TimeIndex);
        }

        private double ComputeExposedFraction(GridSnapshot snapshot)
        {
            double[] values = snapshot.GetValues();
            int exposed = 0;
            for (int i = 0; i < values.Length; i++)
                if (values[i] > Threshold) exposed++;
            return values.Length > 0 ? (double)exposed / values.Length : 0;
        }

        private VectorN BuildObservation()
        {
            double[] values = _lastSnapshot.GetValues();

            double max = 0, sum = 0;
            int exposed = 0;
            for (int i = 0; i < values.Length; i++)
            {
                if (values[i] > max) max = values[i];
                sum += values[i];
                if (values[i] > Threshold) exposed++;
            }

            double mean = values.Length > 0 ? sum / values.Length : 0;
            double exposedFrac = values.Length > 0 ? (double)exposed / values.Length : 0;

            // Normalize wind direction to unit components
            Vector windUnit = _currentWindDirection.GetUnitVector();

            return new VectorN(new[]
            {
                max,
                mean,
                exposedFrac,
                _currentWindSpeed,
                windUnit.x,
                windUnit.y,
                _currentEmissionRate,
                (double)_currentStep / Math.Max(1, _timeFrame.Count - 1)
            });
        }
    }
}
