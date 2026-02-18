using CSharpNumerics.Physics.Applied.Constraints;
using CSharpNumerics.Physics.Objects;
using System.Collections.Generic;

namespace CSharpNumerics.Physics.Applied
{
    /// <summary>
    /// Iterative sequential impulse solver for constraints and joints.
    /// Each iteration applies corrective impulses for all constraints;
    /// more iterations yield tighter constraint satisfaction at higher CPU cost.
    /// Typical: 4–10 iterations for games, 20+ for precise simulations.
    /// </summary>
    public static class ConstraintSolver
    {
        /// <summary>
        /// Solves all constraints for the current time step.
        /// Call this after integrating velocities (v += a·dt) and before integrating positions (x += v·dt).
        /// </summary>
        /// <param name="bodies">Array of all rigid bodies in the simulation.</param>
        /// <param name="constraints">List of constraints to solve.</param>
        /// <param name="dt">Time step in seconds.</param>
        /// <param name="iterations">Number of solver iterations (higher = tighter). Default 10.</param>
        public static void Solve(RigidBody[] bodies, IReadOnlyList<IConstraint> constraints, double dt, int iterations = 10)
        {
            // Pre-step: compute Jacobians, effective masses, bias
            for (int c = 0; c < constraints.Count; c++)
                constraints[c].PreStep(bodies, dt);

            // Iterative solve: converge toward constraint satisfaction
            for (int i = 0; i < iterations; i++)
                for (int c = 0; c < constraints.Count; c++)
                    constraints[c].Solve(bodies);
        }
    }
}
