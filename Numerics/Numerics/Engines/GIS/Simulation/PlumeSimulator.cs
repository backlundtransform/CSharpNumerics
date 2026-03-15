using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics;
using CSharpNumerics.Physics.Enums;
using CSharpNumerics.Physics.Materials;
using NuclearDecay = CSharpNumerics.Physics.Materials.Nuclear.Decay.Decay;
using NuclearDose = CSharpNumerics.Physics.Materials.Nuclear.RadiationDose.RadiationDose;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Simulation
{
    /// <summary>
    /// Mode of physics evaluation used by <see cref="PlumeSimulator"/>.
    /// </summary>
    public enum PlumeMode
    {
        /// <summary>
        /// Steady-state Gaussian plume (same field at every time step).
        /// Quick to compute; suitable when wind/stability are constant.
        /// </summary>
        SteadyState,

        /// <summary>
        /// Transient Gaussian puff that advects downwind over time.
        /// Each time step produces a different snapshot.
        /// </summary>
        Transient
    }

    /// <summary>
    /// Evaluates a Gaussian plume (steady-state or transient puff) from
    /// <see cref="EnvironmentalExtensions"/> on a <see cref="GeoGrid"/>
    /// for every time step in a <see cref="TimeFrame"/>.
    /// <para>
    /// This is the single-scenario deterministic simulator used as the
    /// inner loop of the Monte Carlo engine.
    /// </para>
    /// </summary>
    public class PlumeSimulator
    {
        /// <summary>Source strength Q (kg/s).</summary>
        public double EmissionRate { get; }

        /// <summary>Mean wind speed u (m/s).</summary>
        public double WindSpeed { get; }

        /// <summary>Horizontal wind direction (z component ignored).</summary>
        public Vector WindDirection { get; }

        /// <summary>Effective stack height H (metres).</summary>
        public double StackHeight { get; }

        /// <summary>Source position in world coordinates.</summary>
        public Vector SourcePosition { get; }

        /// <summary>Pasquill–Gifford stability class.</summary>
        public StabilityClass Stability { get; }

        /// <summary>Physics mode: steady-state plume or transient puff.</summary>
        public PlumeMode Mode { get; }

        /// <summary>
        /// Duration of an instantaneous release in seconds (only used in
        /// <see cref="PlumeMode.Transient"/>). Defaults to the time-frame step.
        /// </summary>
        public double ReleaseSeconds { get; set; }

        /// <summary>
        /// Optional radioactive material. When set, each <see cref="GridSnapshot"/>
        /// will carry "activity" (Bq/m³) and "dose" (Sv) layers computed from
        /// the concentration field and the isotope properties.
        /// </summary>
        public MaterialDescriptor Material { get; set; }

        /// <summary>
        /// Creates a plume simulator with the given physical parameters.
        /// </summary>
        public PlumeSimulator(
            double emissionRate,
            double windSpeed,
            Vector windDirection,
            double stackHeight,
            Vector sourcePosition,
            StabilityClass stability = StabilityClass.D,
            PlumeMode mode = PlumeMode.SteadyState)
        {
            if (emissionRate <= 0) throw new ArgumentException("Emission rate must be positive.", nameof(emissionRate));
            if (windSpeed <= 0) throw new ArgumentException("Wind speed must be positive.", nameof(windSpeed));

            EmissionRate = emissionRate;
            WindSpeed = windSpeed;
            WindDirection = windDirection;
            StackHeight = stackHeight;
            SourcePosition = sourcePosition;
            Stability = stability;
            Mode = mode;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Run
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Evaluates the plume on every cell of <paramref name="grid"/>
        /// for every time step in <paramref name="timeFrame"/> and returns
        /// a list of <see cref="GridSnapshot"/> objects.
        /// </summary>
        public List<GridSnapshot> Run(GeoGrid grid, TimeFrame timeFrame)
        {
            var snapshots = new List<GridSnapshot>(timeFrame.Count);
            double[] times = timeFrame.ToArray();

            if (Mode == PlumeMode.SteadyState)
            {
                // Steady-state: evaluate once, replicate across time steps
                ScalarField field = EmissionRate.GaussianPlume(
                    WindSpeed, StackHeight, SourcePosition, WindDirection, Stability);

                double[] values = EvaluateField(field, grid);

                for (int t = 0; t < times.Length; t++)
                {
                    var snap = new GridSnapshot(grid, (double[])values.Clone(), times[t], t);
                    ApplyMaterialLayers(snap, timeFrame.StepSeconds);
                    snapshots.Add(snap);
                }
            }
            else // Transient
            {
                double release = ReleaseSeconds > 0 ? ReleaseSeconds : timeFrame.StepSeconds;

                for (int t = 0; t < times.Length; t++)
                {
                    double time = times[t];

                    ScalarField field = EmissionRate.GaussianPuff(
                        release, WindSpeed, StackHeight,
                        SourcePosition, WindDirection, time, Stability);

                    double[] values = EvaluateField(field, grid);
                    var snap = new GridSnapshot(grid, values, time, t);
                    ApplyMaterialLayers(snap, timeFrame.StepSeconds);
                    snapshots.Add(snap);
                }
            }

            return snapshots;
        }

        /// <summary>
        /// Evaluates the plume for a single time value and returns one snapshot.
        /// </summary>
        public GridSnapshot RunSingle(GeoGrid grid, double time, int timeIndex = 0)
        {
            ScalarField field;

            if (Mode == PlumeMode.SteadyState)
            {
                field = EmissionRate.GaussianPlume(
                    WindSpeed, StackHeight, SourcePosition, WindDirection, Stability);
            }
            else
            {
                double release = ReleaseSeconds > 0 ? ReleaseSeconds : 1.0;
                field = EmissionRate.GaussianPuff(
                    release, WindSpeed, StackHeight,
                    SourcePosition, WindDirection, time, Stability);
            }

            double[] values = EvaluateField(field, grid);
            var snapshot = new GridSnapshot(grid, values, time, timeIndex);
            ApplyMaterialLayers(snapshot, 1.0);
            return snapshot;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Helpers
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Evaluates a scalar field at every cell centre in the grid.
        /// </summary>
        private static double[] EvaluateField(ScalarField field, GeoGrid grid)
        {
            var values = new double[grid.CellCount];
            int i = 0;
            for (int iz = 0; iz < grid.Nz; iz++)
                for (int iy = 0; iy < grid.Ny; iy++)
                    for (int ix = 0; ix < grid.Nx; ix++)
                    {
                        values[i] = field.Evaluate(grid.CellCentre(ix, iy, iz));
                        i++;
                    }
            return values;
        }

        /// <summary>
        /// Computes activity and dose layers from the concentration field
        /// when a material is attached.
        /// </summary>
        private void ApplyMaterialLayers(GridSnapshot snapshot, double timeStepSeconds)
        {
            if (Material == null) return;

            var iso = Material.Isotope;
            var concentrations = snapshot.GetValues();
            int n = concentrations.Length;

            var activity = new double[n];
            var dose = new double[n];

            for (int i = 0; i < n; i++)
            {
                // Activity: concentration (kg/m³) × specific activity (Bq/kg) = Bq/m³
                // Apply decay correction for current time
                double a = NuclearDecay.ConcentrationToActivity(concentrations[i], iso);
                if (snapshot.Time > 0 && !iso.IsStable)
                    a = NuclearDecay.ActivityAtTime(a, snapshot.Time, iso);
                activity[i] = a;

                // Inhalation dose for the time step
                dose[i] = NuclearDose.InhalationDose(
                    a, timeStepSeconds, iso);
            }

            snapshot.SetLayer("activity", activity);
            snapshot.SetLayer("dose", dose);
        }
    }
}
