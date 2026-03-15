using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Engines.GIS.Scenario
{
    /// <summary>
    /// Static entry point for the GIS fluent risk-scenario pipeline.
    /// <para>
    /// Usage:
    /// <code>
    /// var result = RiskScenario
    ///     .ForGaussianPlume(5.0)
    ///     .FromSource(new Vector(0, 0, 50))
    ///     .WithWind(10, new Vector(1, 0, 0))
    ///     .WithStability(StabilityClass.D)
    ///     .OverGrid(new GeoGrid(-500, 500, -500, 500, 0, 100, 10))
    ///     .OverTime(0, 3600, 60)
    ///     .RunMonteCarlo(1000)
    ///     .AnalyzeWith(new KMeans(), new SilhouetteEvaluator(), 2, 6)
    ///     .Build(threshold: 1e-6);
    /// </code>
    /// </para>
    /// </summary>
    public static class RiskScenario
    {
        /// <summary>
        /// Begin building a Gaussian plume risk scenario.
        /// </summary>
        /// <param name="emissionRate">Emission rate Q (kg/s).</param>
        public static RiskScenarioBuilder ForGaussianPlume(double emissionRate)
            => new RiskScenarioBuilder(emissionRate);
    }
}
