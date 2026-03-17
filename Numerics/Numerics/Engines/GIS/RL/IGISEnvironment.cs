using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;

namespace CSharpNumerics.Engines.GIS.RL
{
    /// <summary>
    /// Extends <see cref="IEnvironment"/> with GIS-specific metadata.
    /// Implement this interface to plug any spatial simulation into the
    /// <see cref="ScenarioRLAnalyzer"/> fluent pipeline.
    /// <para>
    /// Existing implementations: <see cref="PlumeEnvironment"/>.
    /// </para>
    /// </summary>
    public interface IGISEnvironment : IEnvironment
    {
        /// <summary>The spatial grid used by this environment.</summary>
        GeoGrid Grid { get; }

        /// <summary>The time discretization used by this environment.</summary>
        TimeFrame TimeFrame { get; }

        /// <summary>Concentration/value threshold for exposure detection.</summary>
        double Threshold { get; set; }

        /// <summary>Per-action cost in the reward function.</summary>
        double ActionCost { get; set; }

        /// <summary>Total time steps per episode.</summary>
        int MaxSteps { get; }

        /// <summary>
        /// The last computed snapshot (available after <see cref="IEnvironment.Reset"/>
        /// or <see cref="IEnvironment.Step(int)"/>).
        /// Returns null if no step has been taken yet.
        /// </summary>
        GridSnapshot LastSnapshot { get; }
    }
}
