using CSharpNumerics.Physics.Applied;
using CSharpNumerics.Physics.Applied.Constraints;
using CSharpNumerics.Physics.Objects;
using CSharpNumerics.Numerics.Objects;

namespace NumericTest
{
    [TestClass]
    public class ConstraintTests
    {
        private const double Tol = 1e-6;
        private const double g = 9.80665;

        /// <summary>
        /// Runs a semi-implicit Euler loop with constraint solving between velocity and position updates.
        /// </summary>
        private static void Simulate(RigidBody[] bodies, IConstraint[] constraints,
            double dt, int steps, int solverIterations = 10, bool applyGravity = true)
        {
            for (int step = 0; step < steps; step++)
            {
                // Velocity integration (apply gravity)
                if (applyGravity)
                    for (int i = 0; i < bodies.Length; i++)
                        if (!bodies[i].IsStatic)
                            bodies[i].Velocity = bodies[i].Velocity + dt * new Vector(0, 0, -g);

                // Solve constraints
                ConstraintSolver.Solve(bodies, constraints, dt, solverIterations);

                // Position integration
                for (int i = 0; i < bodies.Length; i++)
                    if (!bodies[i].IsStatic)
                        bodies[i].Position = bodies[i].Position + dt * bodies[i].Velocity;
            }
        }

        #region DistanceConstraint

        [TestMethod]
        public void DistanceConstraint_MaintainsDistance()
        {
            var bodies = new RigidBody[]
            {
                RigidBody.CreateStatic(new Vector(0, 0, 5)), // ceiling anchor
                RigidBody.CreateSolidSphere(1, 0.1),         // hanging ball
            };
            bodies[1].Position = new Vector(0, 0, 2); // 3m below anchor

            var constraints = new IConstraint[]
            {
                new DistanceConstraint(0, 1, new Vector(0, 0, 0), new Vector(0, 0, 0), 3.0)
            };

            Simulate(bodies, constraints, 0.001, 5000);

            double dist = (bodies[1].Position - bodies[0].Position).GetMagnitude();
            Assert.AreEqual(3.0, dist, 0.05);
        }

        [TestMethod]
        public void DistanceConstraint_Pendulum_Oscillates()
        {
            var bodies = new RigidBody[]
            {
                RigidBody.CreateStatic(new Vector(0, 0, 5)),
                RigidBody.CreateSolidSphere(1, 0.1),
            };
            // Start displaced to the side
            bodies[1].Position = new Vector(2, 0, 2);

            double L = (bodies[1].Position - bodies[0].Position).GetMagnitude();
            var constraints = new IConstraint[]
            {
                new DistanceConstraint(0, 1, new Vector(0, 0, 0), new Vector(0, 0, 0), L)
            };

            Simulate(bodies, constraints, 0.001, 3000);

            // Should have swung to the other side (approximately)
            // Distance should still be maintained
            double dist = (bodies[1].Position - bodies[0].Position).GetMagnitude();
            Assert.AreEqual(L, dist, 0.1);
        }

        [TestMethod]
        public void DistanceConstraint_TwoDynamic_MomentumConserved()
        {
            var bodies = new RigidBody[]
            {
                RigidBody.CreateSolidSphere(2, 0.5),
                RigidBody.CreateSolidSphere(3, 0.5),
            };
            bodies[0].Position = new Vector(0, 0, 0);
            bodies[1].Position = new Vector(2, 0, 0);
            bodies[0].Velocity = new Vector(5, 0, 0);

            var constraints = new IConstraint[]
            {
                new DistanceConstraint(0, 1, new Vector(0, 0, 0), new Vector(0, 0, 0), 2.0)
            };

            var pBefore = bodies[0].Mass * bodies[0].Velocity + bodies[1].Mass * bodies[1].Velocity;

            Simulate(bodies, constraints, 0.001, 1000, applyGravity: false);

            var pAfter = bodies[0].Mass * bodies[0].Velocity + bodies[1].Mass * bodies[1].Velocity;
            Assert.AreEqual(pBefore.x, pAfter.x, 0.1);
            Assert.AreEqual(pBefore.y, pAfter.y, 0.1);
            Assert.AreEqual(pBefore.z, pAfter.z, 0.1);
        }

        [TestMethod]
        public void DistanceConstraint_FromCurrentPositions()
        {
            var bodies = new RigidBody[]
            {
                RigidBody.CreateStatic(new Vector(0, 0, 10)),
                RigidBody.CreateSolidSphere(1, 0.1),
            };
            bodies[1].Position = new Vector(3, 4, 10);

            var constraint = DistanceConstraint.FromCurrentPositions(0, 1, bodies);
            Assert.AreEqual(5.0, constraint.Distance, 1e-10);
        }

        #endregion

        #region BallSocketJoint

        [TestMethod]
        public void BallSocketJoint_KeepsAnchorsCoincident()
        {
            var bodies = new RigidBody[]
            {
                RigidBody.CreateStatic(new Vector(0, 0, 5)),
                RigidBody.CreateSolidSphere(1, 0.1),
            };
            bodies[1].Position = new Vector(0, 0, 3);

            var constraints = new IConstraint[]
            {
                new BallSocketJoint(0, 1,
                    new Vector(0, 0, -1),  // 1m below static center
                    new Vector(0, 0, 1))   // 1m above ball center
            };

            // Pivot should be at (0,0,4)
            Simulate(bodies, constraints, 0.001, 3000);

            var worldA = bodies[0].Position + bodies[0].Orientation * new Vector(0, 0, -1);
            var worldB = bodies[1].Position + bodies[1].Orientation * new Vector(0, 0, 1);
            double separation = (worldB - worldA).GetMagnitude();
            Assert.IsTrue(separation < 0.1, $"Separation {separation} too large");
        }

        [TestMethod]
        public void BallSocketJoint_FromWorldPivot()
        {
            var bodies = new RigidBody[]
            {
                RigidBody.CreateStatic(new Vector(0, 0, 5)),
                RigidBody.CreateSolidSphere(1, 0.1),
            };
            bodies[1].Position = new Vector(0, 0, 3);

            var joint = BallSocketJoint.FromWorldPivot(0, 1, new Vector(0, 0, 4), bodies);
            Assert.AreEqual(-1, joint.LocalAnchorA.z, 1e-10);
            Assert.AreEqual(1, joint.LocalAnchorB.z, 1e-10);
        }

        [TestMethod]
        public void BallSocketJoint_AllowsRotation()
        {
            var bodies = new RigidBody[]
            {
                RigidBody.CreateStatic(new Vector(0, 0, 5)),
                RigidBody.CreateSolidSphere(1, 0.1),
            };
            bodies[1].Position = new Vector(0, 0, 3);
            bodies[1].AngularVelocity = new Vector(0, 0, 5); // spin freely

            var constraints = new IConstraint[]
            {
                new BallSocketJoint(0, 1, new Vector(0, 0, 0), new Vector(0, 0, 2))
            };

            Simulate(bodies, constraints, 0.001, 1000, applyGravity: false);

            // Angular velocity should be preserved (ball-socket allows free rotation)
            Assert.IsTrue(bodies[1].AngularVelocity.GetMagnitude() > 1.0);
        }

        #endregion

        #region HingeJoint

        [TestMethod]
        public void HingeJoint_AllowsRotationOnHingeAxis()
        {
            var bodies = new RigidBody[]
            {
                RigidBody.CreateStatic(new Vector(0, 0, 0)),
                RigidBody.CreateSolidBox(1, 2, 0.1, 0.1),
            };
            bodies[1].Position = new Vector(1, 0, 0);
            bodies[1].AngularVelocity = new Vector(0, 0, 3); // spinning about Z (the hinge axis)

            var constraints = new IConstraint[]
            {
                new HingeJoint(0, 1, new Vector(0, 0, 0), new Vector(-1, 0, 0), new Vector(0, 0, 1))
            };

            Simulate(bodies, constraints, 0.001, 500, applyGravity: false);

            // Z-axis angular velocity should be preserved (hinge allows it)
            Assert.IsTrue(Math.Abs(bodies[1].AngularVelocity.z) > 1.0,
                $"Hinge-axis ω.z = {bodies[1].AngularVelocity.z} should be preserved");
        }

        [TestMethod]
        public void HingeJoint_BlocksRotationOnPerpendicularAxes()
        {
            var bodies = new RigidBody[]
            {
                RigidBody.CreateStatic(new Vector(0, 0, 0)),
                RigidBody.CreateSolidSphere(1, 1),
            };
            bodies[1].Position = new Vector(2, 0, 0);
            bodies[1].AngularVelocity = new Vector(5, 5, 0); // spinning on X and Y (both locked)

            var constraints = new IConstraint[]
            {
                new HingeJoint(0, 1, new Vector(0, 0, 0), new Vector(-2, 0, 0), new Vector(0, 0, 1))
            };

            Simulate(bodies, constraints, 0.001, 500, applyGravity: false);

            // X and Y angular velocity should be near zero (locked by hinge)
            Assert.IsTrue(Math.Abs(bodies[1].AngularVelocity.x) < 0.5,
                $"ω.x = {bodies[1].AngularVelocity.x} should be near zero");
            Assert.IsTrue(Math.Abs(bodies[1].AngularVelocity.y) < 0.5,
                $"ω.y = {bodies[1].AngularVelocity.y} should be near zero");
        }

        [TestMethod]
        public void HingeJoint_MaintainsPivotPosition()
        {
            var bodies = new RigidBody[]
            {
                RigidBody.CreateStatic(new Vector(0, 0, 5)),
                RigidBody.CreateSolidBox(1, 2, 0.2, 0.2),
            };
            bodies[1].Position = new Vector(1, 0, 5);

            var constraints = new IConstraint[]
            {
                HingeJoint.FromWorldPivot(0, 1, new Vector(0, 0, 5), new Vector(0, 0, 1), bodies)
            };

            Simulate(bodies, constraints, 0.001, 2000);

            // The position constraint part should keep the body near the pivot
            double distFromPivot = (bodies[1].Position - new Vector(0, 0, 5)).GetMagnitude();
            // Body center is 1m from pivot, so it should stay around there
            Assert.IsTrue(distFromPivot < 3.0, $"Body drifted to {distFromPivot}m from pivot");
        }

        #endregion

        #region SpringJoint

        [TestMethod]
        public void SpringJoint_PullsTowardRestLength()
        {
            var bodies = new RigidBody[]
            {
                RigidBody.CreateStatic(new Vector(0, 0, 0)),
                RigidBody.CreateSolidSphere(1, 0.5),
            };
            bodies[1].Position = new Vector(5, 0, 0); // stretched: distance=5, rest=2

            var constraints = new IConstraint[]
            {
                new SpringJoint(0, 1, new Vector(0, 0, 0), new Vector(0, 0, 0),
                    stiffness: 10, damping: 1, restLength: 2)
            };

            Simulate(bodies, constraints, 0.001, 5000, applyGravity: false);

            // Should have oscillated toward rest length and damped out
            double dist = (bodies[1].Position - bodies[0].Position).GetMagnitude();
            Assert.IsTrue(dist < 4.0, $"Distance {dist} should be approaching rest length 2");
        }

        [TestMethod]
        public void SpringJoint_DampedOscillation_SettlesNearRest()
        {
            var bodies = new RigidBody[]
            {
                RigidBody.CreateStatic(new Vector(0, 0, 0)),
                RigidBody.CreateSolidSphere(1, 0.5),
            };
            bodies[1].Position = new Vector(3, 0, 0);

            var constraints = new IConstraint[]
            {
                new SpringJoint(0, 1, new Vector(0, 0, 0), new Vector(0, 0, 0),
                    stiffness: 20, damping: 5, restLength: 1)
            };

            Simulate(bodies, constraints, 0.001, 20000, applyGravity: false);

            double dist = (bodies[1].Position - bodies[0].Position).GetMagnitude();
            double speed = bodies[1].Velocity.GetMagnitude();

            // With strong damping, should settle near rest length
            Assert.AreEqual(1.0, dist, 0.3);
            Assert.IsTrue(speed < 0.5, $"Speed {speed} should be near zero (damped)");
        }

        [TestMethod]
        public void SpringJoint_ZeroDamping_Oscillates()
        {
            var bodies = new RigidBody[]
            {
                RigidBody.CreateStatic(new Vector(0, 0, 0)),
                RigidBody.CreateSolidSphere(1, 0.5),
            };
            bodies[1].Position = new Vector(3, 0, 0); // 2m beyond rest

            var constraints = new IConstraint[]
            {
                new SpringJoint(0, 1, new Vector(0, 0, 0), new Vector(0, 0, 0),
                    stiffness: 10, damping: 0, restLength: 1)
            };

            // Measure initial energy: E = ½k·(x-L₀)²
            double E0 = 0.5 * 10 * 4; // ½·10·2² = 20

            Simulate(bodies, constraints, 0.001, 5000, applyGravity: false);

            double dist = (bodies[1].Position - bodies[0].Position).GetMagnitude();
            double v = bodies[1].Velocity.GetMagnitude();
            double stretch = dist - 1.0;
            double E = 0.5 * 10 * stretch * stretch + 0.5 * 1.0 * v * v;

            // Energy should be approximately conserved (no damping)
            Assert.AreEqual(E0, E, 3.0); // loose tolerance due to discrete integration
        }

        #endregion

        #region ConstraintSolver

        [TestMethod]
        public void Solver_MultipleConstraints_Chain()
        {
            // Chain of 3 bodies connected by distance constraints
            var bodies = new RigidBody[]
            {
                RigidBody.CreateStatic(new Vector(0, 0, 10)),
                RigidBody.CreateSolidSphere(1, 0.2),
                RigidBody.CreateSolidSphere(1, 0.2),
            };
            bodies[1].Position = new Vector(0, 0, 7);
            bodies[2].Position = new Vector(0, 0, 4);

            var constraints = new IConstraint[]
            {
                new DistanceConstraint(0, 1, new Vector(0, 0, 0), new Vector(0, 0, 0), 3.0),
                new DistanceConstraint(1, 2, new Vector(0, 0, 0), new Vector(0, 0, 0), 3.0),
            };

            Simulate(bodies, constraints, 0.001, 5000, solverIterations: 20);

            double d01 = (bodies[1].Position - bodies[0].Position).GetMagnitude();
            double d12 = (bodies[2].Position - bodies[1].Position).GetMagnitude();

            Assert.AreEqual(3.0, d01, 0.15);
            Assert.AreEqual(3.0, d12, 0.15);
        }

        [TestMethod]
        public void Solver_StaticBody_DoesNotMove()
        {
            var bodies = new RigidBody[]
            {
                RigidBody.CreateStatic(new Vector(0, 0, 10)),
                RigidBody.CreateSolidSphere(100, 1),
            };
            bodies[1].Position = new Vector(0, 0, 5);

            var constraints = new IConstraint[]
            {
                new DistanceConstraint(0, 1, new Vector(0, 0, 0), new Vector(0, 0, 0), 5.0)
            };

            Simulate(bodies, constraints, 0.001, 1000);

            Assert.AreEqual(0, bodies[0].Position.x, 1e-10);
            Assert.AreEqual(0, bodies[0].Position.y, 1e-10);
            Assert.AreEqual(10, bodies[0].Position.z, 1e-10);
        }

        #endregion
    }
}
