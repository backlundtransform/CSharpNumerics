using CSharpNumerics.Physics.Applied.BroadPhase;
using CSharpNumerics.Physics.Applied.Constraints;
using CSharpNumerics.Physics.Applied.Objects;
using CSharpNumerics.Physics.Objects;
using Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Physics.Applied
{
    /// <summary>
    /// Container and orchestrator for a physics simulation.
    /// Manages bodies, constraints, and the full simulation pipeline:
    /// integrate velocities → solve constraints → integrate positions → detect &amp; resolve collisions.
    /// 
    /// Use <see cref="Step"/> for a single fixed step, or <see cref="Update"/> with the
    /// fixed-timestep accumulator for framerate-independent simulation.
    /// </summary>
    public class PhysicsWorld
    {
        private RigidBody[] _bodies;
        private double[] _radii;
        private int _bodyCount;

        private readonly List<IConstraint> _constraints = new();
        private readonly List<(int a, int b)> _pairs = new();
        private readonly IBroadPhase _broadPhase;
        private double _accumulator;

        #region Settings

        /// <summary>Global gravity applied to all dynamic bodies each step. Default: (0, 0, -9.80665).</summary>
        public Vector Gravity { get; set; } = new Vector(0, 0, -9.80665);

        /// <summary>Default coefficient of restitution for collisions (0 = sticky, 1 = elastic). Default 0.5.</summary>
        public double DefaultRestitution { get; set; } = 0.5;

        /// <summary>Default friction coefficient for collisions. Default 0.3.</summary>
        public double DefaultFriction { get; set; } = 0.3;

        /// <summary>Number of constraint solver iterations per step. Default 10.</summary>
        public int SolverIterations { get; set; } = 10;

        /// <summary>Fixed time step for the accumulator-based <see cref="Update"/> method. Default 1/60 s.</summary>
        public double FixedTimeStep { get; set; } = 1.0 / 60.0;

        #endregion

        #region Events

        /// <summary>
        /// Called for each collision detected during a step.
        /// Parameters: (bodyIndexA, bodyIndexB, contactPoint).
        /// </summary>
        public Action<int, int, ContactPoint>? OnCollision;

        #endregion

        /// <summary>
        /// Creates a new physics world.
        /// </summary>
        /// <param name="broadPhase">Broad phase algorithm. Defaults to <see cref="SweepAndPruneBroadPhase"/>.</param>
        /// <param name="initialCapacity">Initial body array capacity.</param>
        public PhysicsWorld(IBroadPhase? broadPhase = null, int initialCapacity = 64)
        {
            _broadPhase = broadPhase ?? new SweepAndPruneBroadPhase();
            _bodies = new RigidBody[initialCapacity];
            _radii = new double[initialCapacity];
        }

        #region Body Management

        /// <summary>
        /// Adds a rigid body to the world and returns its index.
        /// </summary>
        /// <param name="body">The rigid body to add.</param>
        /// <param name="boundingRadius">Bounding sphere radius for collision detection.</param>
        /// <returns>The body's index, used for constraints and queries.</returns>
        public int AddBody(RigidBody body, double boundingRadius = 1.0)
        {
            if (_bodyCount >= _bodies.Length)
            {
                int newSize = _bodies.Length * 2;
                Array.Resize(ref _bodies, newSize);
                Array.Resize(ref _radii, newSize);
            }
            int index = _bodyCount;
            _bodies[index] = body;
            _radii[index] = boundingRadius;
            _bodyCount++;
            return index;
        }

        /// <summary>Returns a reference to the body at the given index for direct modification.</summary>
        public ref RigidBody Body(int index) => ref _bodies[index];

        /// <summary>Returns the bounding radius of the body at the given index.</summary>
        public double BoundingRadius(int index) => _radii[index];

        /// <summary>Number of bodies in the world.</summary>
        public int BodyCount => _bodyCount;

        #endregion

        #region Constraint Management

        /// <summary>Adds a constraint or joint to the world.</summary>
        public void AddConstraint(IConstraint constraint) => _constraints.Add(constraint);

        /// <summary>Number of constraints in the world.</summary>
        public int ConstraintCount => _constraints.Count;

        #endregion

        #region Simulation

        /// <summary>
        /// Advances the simulation by one fixed time step.
        /// Pipeline: integrate velocities → solve constraints → integrate positions → detect &amp; resolve collisions.
        /// </summary>
        /// <param name="dt">Time step in seconds.</param>
        public void Step(double dt)
        {
            // 1. Integrate velocities (apply gravity)
            for (int i = 0; i < _bodyCount; i++)
            {
                if (_bodies[i].IsStatic) continue;
                _bodies[i].Velocity = _bodies[i].Velocity + dt * Gravity;
            }

            // 2. Solve constraints (velocity correction)
            if (_constraints.Count > 0)
                ConstraintSolver.Solve(_bodies, _constraints, dt, SolverIterations);

            // 3. Integrate positions
            for (int i = 0; i < _bodyCount; i++)
            {
                if (_bodies[i].IsStatic) continue;
                _bodies[i].Position = _bodies[i].Position + dt * _bodies[i].Velocity;
            }

            // 4. Collision detection & response (on integrated positions)
            _broadPhase.FindPairs(_bodies, _radii, _bodyCount, _pairs);

            for (int p = 0; p < _pairs.Count; p++)
            {
                var (a, b) = _pairs[p];
                var sA = new BoundingSphere(_bodies[a].Position, _radii[a]);
                var sB = new BoundingSphere(_bodies[b].Position, _radii[b]);
                var contact = sA.SphereSphereContact(sB);

                if (contact is ContactPoint c)
                {
                    CollisionResponse.ResolveCollision(
                        ref _bodies[a], ref _bodies[b], c,
                        DefaultRestitution, DefaultFriction);

                    CollisionResponse.CorrectPositions(
                        ref _bodies[a], ref _bodies[b], c);

                    OnCollision?.Invoke(a, b, c);
                }
            }
        }

        /// <summary>
        /// Advances the simulation using a fixed-timestep accumulator.
        /// Consumes <paramref name="elapsed"/> time in fixed-size steps, ensuring deterministic
        /// and framerate-independent behavior.
        /// </summary>
        /// <param name="elapsed">Wall-clock time elapsed since last call (seconds).</param>
        /// <returns>Number of physics steps executed.</returns>
        public int Update(double elapsed)
        {
            _accumulator += elapsed;
            int steps = 0;
            while (_accumulator >= FixedTimeStep)
            {
                Step(FixedTimeStep);
                _accumulator -= FixedTimeStep;
                steps++;
            }
            return steps;
        }

        /// <summary>
        /// Returns the interpolation factor for rendering between the last and next physics step.
        /// Use this to smooth visual positions: renderPos = lerp(prevPos, currPos, Alpha).
        /// </summary>
        public double Alpha => _accumulator / FixedTimeStep;

        #endregion
    }
}
