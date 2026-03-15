using CSharpNumerics.Physics.Enums;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Simulation
{
    /// <summary>
    /// Defines parameter variation ranges for Monte Carlo plume simulations.
    /// Each property pair (min/max or weights) describes the stochastic range
    /// from which the Monte Carlo runner will sample.
    /// <para>
    /// When a range is <c>null</c> the corresponding parameter is held at
    /// its deterministic baseline value from the builder.
    /// </para>
    /// </summary>
    public class ScenarioVariation
    {
        /// <summary>Minimum wind speed (m/s). Null = use baseline.</summary>
        public double? WindSpeedMin { get; private set; }

        /// <summary>Maximum wind speed (m/s). Null = use baseline.</summary>
        public double? WindSpeedMax { get; private set; }

        /// <summary>
        /// Standard deviation (degrees) of wind direction jitter
        /// applied as a random rotation around the z-axis.
        /// Null or 0 = no jitter.
        /// </summary>
        public double? WindDirectionJitterDeg { get; private set; }

        /// <summary>Minimum emission rate (kg/s). Null = use baseline.</summary>
        public double? EmissionRateMin { get; private set; }

        /// <summary>Maximum emission rate (kg/s). Null = use baseline.</summary>
        public double? EmissionRateMax { get; private set; }

        /// <summary>
        /// Categorical weights for each <see cref="StabilityClass"/>.
        /// The weights are normalised internally — they need not sum to 1.
        /// Null = always use the baseline stability class.
        /// </summary>
        public Dictionary<StabilityClass, double> StabilityWeights { get; private set; }

        // ═══════════════════════════════════════════════════════════════
        //  Fluent setters
        // ═══════════════════════════════════════════════════════════════

        /// <summary>Vary wind speed uniformly in [min, max] m/s.</summary>
        public ScenarioVariation WindSpeed(double min, double max)
        {
            if (min > max) throw new ArgumentException("min must be <= max.");
            WindSpeedMin = min;
            WindSpeedMax = max;
            return this;
        }

        /// <summary>Apply Gaussian jitter to the wind direction with the given standard deviation in degrees.</summary>
        public ScenarioVariation WindDirectionJitter(double stdDevDegrees)
        {
            WindDirectionJitterDeg = stdDevDegrees;
            return this;
        }

        /// <summary>Vary emission rate uniformly in [min, max] kg/s.</summary>
        public ScenarioVariation EmissionRate(double min, double max)
        {
            if (min > max) throw new ArgumentException("min must be <= max.");
            EmissionRateMin = min;
            EmissionRateMax = max;
            return this;
        }

        /// <summary>
        /// Set categorical weights for stability classes.
        /// Weights are normalised; they need not sum to 1.
        /// </summary>
        public ScenarioVariation SetStabilityWeights(
            double a = 0, double b = 0, double c = 0,
            double d = 0, double e = 0, double f = 0)
        {
            StabilityWeights = new Dictionary<StabilityClass, double>
            {
                [StabilityClass.A] = a,
                [StabilityClass.B] = b,
                [StabilityClass.C] = c,
                [StabilityClass.D] = d,
                [StabilityClass.E] = e,
                [StabilityClass.F] = f,
            };
            return this;
        }

        /// <summary>True when at least one variation range has been set.</summary>
        public bool HasVariation =>
            WindSpeedMin.HasValue ||
            (WindDirectionJitterDeg.HasValue && WindDirectionJitterDeg.Value > 0) ||
            EmissionRateMin.HasValue ||
            StabilityWeights != null;
    }
}
