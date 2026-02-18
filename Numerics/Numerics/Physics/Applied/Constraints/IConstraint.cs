using CSharpNumerics.Physics.Objects;

namespace CSharpNumerics.Physics.Applied.Constraints
{
    /// <summary>
    /// Interface for a velocity-level constraint solved by the sequential impulse method.
    /// Each constraint corrects body velocities to satisfy its kinematic relationship.
    /// </summary>
    public interface IConstraint
    {
        /// <summary>
        /// Precomputes Jacobians, effective mass, and Baumgarte bias for the current time step.
        /// Called once per frame before the iterative solve loop.
        /// </summary>
        /// <param name="bodies">Array of all rigid bodies in the simulation.</param>
        /// <param name="dt">Time step in seconds.</param>
        void PreStep(RigidBody[] bodies, double dt);

        /// <summary>
        /// Applies one corrective impulse to the constrained bodies.
        /// Called multiple times per frame (typically 4–20 iterations) for convergence.
        /// </summary>
        /// <param name="bodies">Array of all rigid bodies in the simulation.</param>
        void Solve(RigidBody[] bodies);
    }
}
