using CSharpNumerics.Physics.Applied;
using CSharpNumerics.Physics.Applied.BroadPhase;
using CSharpNumerics.Physics.Applied.Constraints;
using CSharpNumerics.Physics.Applied.Objects;
using CSharpNumerics.Physics.Objects;
using CSharpNumerics.Numerics.Objects;

namespace NumericTest
{
    [TestClass]
    public class PhysicsWorldTests
    {
        #region BroadPhase

        [TestMethod]
        public void BruteForce_FindsOverlappingPairs()
        {
            var bodies = new RigidBody[3];
            bodies[0] = RigidBody.CreateSolidSphere(1, 1); bodies[0].Position = new Vector(0, 0, 0);
            bodies[1] = RigidBody.CreateSolidSphere(1, 1); bodies[1].Position = new Vector(1.5, 0, 0);
            bodies[2] = RigidBody.CreateSolidSphere(1, 1); bodies[2].Position = new Vector(10, 0, 0);
            var radii = new double[] { 1, 1, 1 };

            var bp = new BruteForceBroadPhase();
            var pairs = new List<(int, int)>();
            bp.FindPairs(bodies, radii, 3, pairs);

            Assert.AreEqual(1, pairs.Count); // only 0-1 overlap
            Assert.IsTrue(pairs.Contains((0, 1)));
        }

        [TestMethod]
        public void BruteForce_SkipsStaticStaticPairs()
        {
            var bodies = new RigidBody[2];
            bodies[0] = RigidBody.CreateStatic(new Vector(0, 0, 0));
            bodies[1] = RigidBody.CreateStatic(new Vector(0.5, 0, 0));
            var radii = new double[] { 1, 1 };

            var bp = new BruteForceBroadPhase();
            var pairs = new List<(int, int)>();
            bp.FindPairs(bodies, radii, 2, pairs);

            Assert.AreEqual(0, pairs.Count);
        }

        [TestMethod]
        public void SweepAndPrune_SameResultsAsBruteForce()
        {
            var rng = new Random(42);
            var bodies = new RigidBody[20];
            var radii = new double[20];
            for (int i = 0; i < 20; i++)
            {
                bodies[i] = RigidBody.CreateSolidSphere(1, 0.5);
                bodies[i].Position = new Vector(rng.NextDouble() * 10, rng.NextDouble() * 10, rng.NextDouble() * 10);
                radii[i] = 1.0;
            }

            var brute = new BruteForceBroadPhase();
            var sap = new SweepAndPruneBroadPhase();
            var bruteResults = new List<(int, int)>();
            var sapResults = new List<(int, int)>();

            brute.FindPairs(bodies, radii, 20, bruteResults);
            sap.FindPairs(bodies, radii, 20, sapResults);

            // Normalize pairs for comparison (sort each pair, then sort list)
            var normalize = (List<(int, int)> list) =>
            {
                var sorted = list.Select(p => p.Item1 < p.Item2 ? p : (p.Item2, p.Item1)).ToList();
                sorted.Sort((a, b) => a.Item1 != b.Item1 ? a.Item1.CompareTo(b.Item1) : a.Item2.CompareTo(b.Item2));
                return sorted;
            };

            var nb = normalize(bruteResults);
            var ns = normalize(sapResults);

            Assert.AreEqual(nb.Count, ns.Count, $"Brute: {nb.Count}, SAP: {ns.Count}");
            for (int i = 0; i < nb.Count; i++)
                Assert.AreEqual(nb[i], ns[i]);
        }

        [TestMethod]
        public void SweepAndPrune_NoOverlap_EmptyResults()
        {
            var bodies = new RigidBody[3];
            bodies[0] = RigidBody.CreateSolidSphere(1, 1); bodies[0].Position = new Vector(0, 0, 0);
            bodies[1] = RigidBody.CreateSolidSphere(1, 1); bodies[1].Position = new Vector(10, 0, 0);
            bodies[2] = RigidBody.CreateSolidSphere(1, 1); bodies[2].Position = new Vector(20, 0, 0);
            var radii = new double[] { 1, 1, 1 };

            var sap = new SweepAndPruneBroadPhase();
            var pairs = new List<(int, int)>();
            sap.FindPairs(bodies, radii, 3, pairs);

            Assert.AreEqual(0, pairs.Count);
        }

        #endregion

        #region PhysicsWorld — Setup

        [TestMethod]
        public void World_AddBody_ReturnsSequentialIndices()
        {
            var world = new PhysicsWorld();

            int i0 = world.AddBody(RigidBody.CreateSolidSphere(1, 1), 1);
            int i1 = world.AddBody(RigidBody.CreateSolidSphere(2, 1), 1);
            int i2 = world.AddBody(RigidBody.CreateStatic(new Vector(0, 0, 0)), 1);

            Assert.AreEqual(0, i0);
            Assert.AreEqual(1, i1);
            Assert.AreEqual(2, i2);
            Assert.AreEqual(3, world.BodyCount);
        }

        [TestMethod]
        public void World_BodyRef_AllowsDirectModification()
        {
            var world = new PhysicsWorld();
            var body = RigidBody.CreateSolidSphere(1, 1);
            body.Position = new Vector(0, 0, 0);
            int idx = world.AddBody(body, 1);

            ref var b = ref world.Body(idx);
            b.Velocity = new Vector(5, 0, 0);

            Assert.AreEqual(5, world.Body(idx).Velocity.x, 1e-10);
        }

        #endregion

        #region PhysicsWorld — Gravity

        [TestMethod]
        public void World_Step_GravityAcceleratesBodies()
        {
            var world = new PhysicsWorld { Gravity = new Vector(0, 0, -10) };

            var ball = RigidBody.CreateSolidSphere(1, 1);
            ball.Position = new Vector(0, 0, 100);
            world.AddBody(ball, 1);

            double dt = 0.01;
            for (int i = 0; i < 100; i++) // 1 second
                world.Step(dt);

            // After 1s under g=-10: vz ≈ -10, z ≈ 100 - 5 = 95
            Assert.AreEqual(-10, world.Body(0).Velocity.z, 0.5);
            Assert.AreEqual(95, world.Body(0).Position.z, 1.0);
        }

        [TestMethod]
        public void World_Step_StaticBodyDoesNotMove()
        {
            var world = new PhysicsWorld();

            var wall = RigidBody.CreateStatic(new Vector(5, 0, 0));
            world.AddBody(wall, 10);

            for (int i = 0; i < 100; i++)
                world.Step(0.01);

            Assert.AreEqual(5, world.Body(0).Position.x, 1e-10);
            Assert.AreEqual(0, world.Body(0).Velocity.GetMagnitude(), 1e-10);
        }

        #endregion

        #region PhysicsWorld — Collisions

        [TestMethod]
        public void World_Step_SpheresCollide_AndBounce()
        {
            var world = new PhysicsWorld
            {
                Gravity = new Vector(0, 0, 0),
                DefaultRestitution = 1.0,
                DefaultFriction = 0
            };

            var a = RigidBody.CreateSolidSphere(1, 1);
            a.Position = new Vector(0, 0, 0);
            a.Velocity = new Vector(5, 0, 0);
            world.AddBody(a, 1);

            var b = RigidBody.CreateSolidSphere(1, 1);
            b.Position = new Vector(3, 0, 0);
            b.Velocity = new Vector(0, 0, 0);
            world.AddBody(b, 1);

            // Run a few steps to let them collide
            for (int i = 0; i < 200; i++)
                world.Step(0.001);

            // After elastic collision of equal masses: velocities should swap
            // Body A should be slow/stopped, body B should be moving +x
            Assert.IsTrue(world.Body(1).Velocity.x > 2.0,
                $"Body B vx = {world.Body(1).Velocity.x} should be moving right");
        }

        [TestMethod]
        public void World_OnCollision_CallbackFires()
        {
            var world = new PhysicsWorld
            {
                Gravity = new Vector(0, 0, 0),
                DefaultRestitution = 1.0
            };

            var a = RigidBody.CreateSolidSphere(1, 1);
            a.Position = new Vector(0, 0, 0);
            a.Velocity = new Vector(10, 0, 0);
            world.AddBody(a, 1);

            var b = RigidBody.CreateSolidSphere(1, 1);
            b.Position = new Vector(1.5, 0, 0);
            world.AddBody(b, 1);

            int collisionCount = 0;
            world.OnCollision = (ia, ib, contact) => collisionCount++;

            world.Step(0.01);

            Assert.IsTrue(collisionCount > 0, "OnCollision should have fired");
        }

        [TestMethod]
        public void World_BallOnFloor_Bounces()
        {
            var world = new PhysicsWorld
            {
                Gravity = new Vector(0, 0, -10),
                DefaultRestitution = 0.8,
                DefaultFriction = 0
            };

            // Floor
            var floor = RigidBody.CreateStatic(new Vector(0, 0, 0));
            world.AddBody(floor, 100); // large radius for floor

            // Ball dropping from height
            var ball = RigidBody.CreateSolidSphere(1, 1);
            ball.Position = new Vector(0, 0, 50);
            world.AddBody(ball, 1);

            for (int i = 0; i < 2000; i++)
                world.Step(0.001);

            // Ball should have bounced and still be above floor
            Assert.IsTrue(world.Body(1).Position.z > -2, $"Ball z = {world.Body(1).Position.z}");
        }

        #endregion

        #region PhysicsWorld — Constraints

        [TestMethod]
        public void World_Step_ConstraintMaintainsDistance()
        {
            var world = new PhysicsWorld
            {
                Gravity = new Vector(0, 0, -10),
                SolverIterations = 15
            };

            var anchor = RigidBody.CreateStatic(new Vector(0, 0, 10));
            int iA = world.AddBody(anchor, 0.1);

            var bob = RigidBody.CreateSolidSphere(1, 0.5);
            bob.Position = new Vector(3, 0, 10);
            int iB = world.AddBody(bob, 0.5);

            world.AddConstraint(new DistanceConstraint(iA, iB,
                new Vector(0, 0, 0), new Vector(0, 0, 0), 3.0));

            for (int i = 0; i < 3000; i++)
                world.Step(0.001);

            double dist = (world.Body(iB).Position - world.Body(iA).Position).GetMagnitude();
            Assert.AreEqual(3.0, dist, 0.15);
        }

        #endregion

        #region Fixed Timestep Accumulator

        [TestMethod]
        public void World_Update_FixedTimestep_Deterministic()
        {
            var world1 = new PhysicsWorld
            {
                Gravity = new Vector(0, 0, -10),
                FixedTimeStep = 0.01
            };
            var ball1 = RigidBody.CreateSolidSphere(1, 1);
            ball1.Position = new Vector(0, 0, 100);
            world1.AddBody(ball1, 1);

            var world2 = new PhysicsWorld
            {
                Gravity = new Vector(0, 0, -10),
                FixedTimeStep = 0.01
            };
            var ball2 = RigidBody.CreateSolidSphere(1, 1);
            ball2.Position = new Vector(0, 0, 100);
            world2.AddBody(ball2, 1);

            // Simulate same total time with different frame rates
            // World 1: 60 fps-ish (16.67ms frames)
            for (int i = 0; i < 60; i++)
                world1.Update(1.0 / 60.0);

            // World 2: 30 fps-ish (33.33ms frames)
            for (int i = 0; i < 30; i++)
                world2.Update(1.0 / 30.0);

            // Both should have run ~100 physics steps (1s / 0.01 = 100)
            // Positions should be identical (deterministic)
            Assert.AreEqual(world1.Body(0).Position.z, world2.Body(0).Position.z, 0.1);
            Assert.AreEqual(world1.Body(0).Velocity.z, world2.Body(0).Velocity.z, 0.1);
        }

        [TestMethod]
        public void World_Update_ReturnsStepCount()
        {
            var world = new PhysicsWorld { FixedTimeStep = 0.01 };
            world.AddBody(RigidBody.CreateSolidSphere(1, 1), 1);

            int steps = world.Update(0.035); // should run 3 steps (0.03), leave 0.005 in accumulator
            Assert.AreEqual(3, steps);

            steps = world.Update(0.006); // accumulator = 0.011 → 1 step
            Assert.AreEqual(1, steps);
        }

        [TestMethod]
        public void World_Alpha_InterpolationFactor()
        {
            var world = new PhysicsWorld { FixedTimeStep = 0.01 };
            world.AddBody(RigidBody.CreateSolidSphere(1, 1), 1);

            world.Update(0.015); // 1 step, 0.005 left
            Assert.AreEqual(0.5, world.Alpha, 0.01); // 0.005 / 0.01 = 0.5
        }

        #endregion

        #region PhysicsWorld — Full Scenario

        [TestMethod]
        public void World_FullScenario_ConstrainedPendulumWithCollision()
        {
            var world = new PhysicsWorld
            {
                Gravity = new Vector(0, 0, -10),
                DefaultRestitution = 0.8,
                SolverIterations = 15
            };

            // Static anchor
            var anchor = RigidBody.CreateStatic(new Vector(0, 0, 10));
            int iAnchor = world.AddBody(anchor, 0.1);

            // Pendulum bob
            var bob = RigidBody.CreateSolidSphere(1, 0.5);
            bob.Position = new Vector(3, 0, 10); // displaced sideways
            int iBob = world.AddBody(bob, 0.5);

            // Target on the ground
            var target = RigidBody.CreateSolidSphere(2, 0.5);
            target.Position = new Vector(0, 0, 7);
            int iTarget = world.AddBody(target, 0.5);

            // Pendulum constraint
            world.AddConstraint(new DistanceConstraint(iAnchor, iBob,
                new Vector(0, 0, 0), new Vector(0, 0, 0), 3.0));

            bool collisionDetected = false;
            world.OnCollision = (a, b, c) =>
            {
                if ((a == iBob && b == iTarget) || (a == iTarget && b == iBob))
                    collisionDetected = true;
            };

            // Run simulation
            for (int i = 0; i < 5000; i++)
                world.Step(0.001);

            // Verify constraint held
            double dist = (world.Body(iBob).Position - world.Body(iAnchor).Position).GetMagnitude();
            Assert.AreEqual(3.0, dist, 0.2);

            // Anchor didn't move
            Assert.AreEqual(10, world.Body(iAnchor).Position.z, 1e-10);
        }

        #endregion
    }
}
