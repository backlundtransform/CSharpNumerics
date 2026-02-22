using CSharpNumerics.Statistics.Random;

namespace CSharpNumerics.Statistics.MonteCarlo
{
    /// <summary>
    /// Defines a model that can be evaluated by the Monte Carlo simulator.
    /// Each call to <see cref="Evaluate"/> represents one trial where the model
    /// draws from its internal distributions and returns a scalar outcome.
    /// </summary>
    public interface IMonteCarloModel
    {
        /// <summary>
        /// Runs a single stochastic trial and returns a scalar result.
        /// The implementation should use the provided <see cref="RandomGenerator"/>
        /// for all random decisions to ensure reproducibility when seeded.
        /// </summary>
        /// <param name="rng">The random generator to use for this trial.</param>
        /// <returns>A scalar outcome of the trial.</returns>
        double Evaluate(RandomGenerator rng);
    }
}
