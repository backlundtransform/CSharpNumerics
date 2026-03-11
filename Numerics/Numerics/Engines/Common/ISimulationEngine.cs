namespace CSharpNumerics.Engines.Common
{
    /// <summary>
    /// Common interface for all simulation engines (Game, Audio, GIS, …).
    /// </summary>
    public interface ISimulationEngine
    {
        /// <summary>Current simulation time in seconds.</summary>
        double Time { get; }

        /// <summary>True after <see cref="Init"/> has been called.</summary>
        bool IsInitialized { get; }

        /// <summary>Initialize / prepare the engine for simulation.</summary>
        void Init();

        /// <summary>Advance the simulation by <paramref name="dt"/> seconds.</summary>
        void Step(double dt);

        /// <summary>Reset to initial state.</summary>
        void Reset();
    }
}
