using CSharpNumerics.Physics.Objects;
using System.Collections.Generic;

namespace CSharpNumerics.Physics.Applied.BroadPhase
{
    /// <summary>
    /// Interface for a broad-phase collision detection algorithm.
    /// Returns candidate pairs of body indices that may be colliding.
    /// The narrow phase then performs exact contact generation on these pairs.
    /// </summary>
    public interface IBroadPhase
    {
        /// <summary>
        /// Finds all pairs of bodies whose bounding spheres overlap.
        /// Results are appended to <paramref name="results"/> (cleared first).
        /// </summary>
        /// <param name="bodies">Array of rigid bodies.</param>
        /// <param name="radii">Parallel array of bounding radii per body.</param>
        /// <param name="count">Number of active bodies.</param>
        /// <param name="results">Output list of candidate pairs (body indices).</param>
        void FindPairs(RigidBody[] bodies, double[] radii, int count, List<(int a, int b)> results);
    }
}
