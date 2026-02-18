using CSharpNumerics.Physics.Applied;
using CSharpNumerics.Physics.Applied.Objects;
using CSharpNumerics.Physics.Objects;
using Numerics.Objects;

namespace NumericTest
{
    [TestClass]
    public class CollisionTests
    {
        private const double Tol = 1e-10;

        #region AABB

        [TestMethod]
        public void AABB_CenterAndHalfExtents()
        {
            var box = new AABB(new Vector(-1, -2, -3), new Vector(1, 2, 3));

            Assert.AreEqual(0, box.Center.x, Tol);
            Assert.AreEqual(0, box.Center.y, Tol);
            Assert.AreEqual(0, box.Center.z, Tol);

            Assert.AreEqual(1, box.HalfExtents.x, Tol);
            Assert.AreEqual(2, box.HalfExtents.y, Tol);
            Assert.AreEqual(3, box.HalfExtents.z, Tol);
        }

        [TestMethod]
        public void AABB_FromCenterExtents()
        {
            var box = AABB.FromCenterExtents(new Vector(5, 5, 5), new Vector(1, 2, 3));

            Assert.AreEqual(4, box.Min.x, Tol);
            Assert.AreEqual(3, box.Min.y, Tol);
            Assert.AreEqual(2, box.Min.z, Tol);
            Assert.AreEqual(6, box.Max.x, Tol);
            Assert.AreEqual(7, box.Max.y, Tol);
            Assert.AreEqual(8, box.Max.z, Tol);
        }

        [TestMethod]
        public void AABB_Contains_PointInside()
        {
            var box = new AABB(new Vector(0, 0, 0), new Vector(10, 10, 10));

            Assert.IsTrue(box.Contains(new Vector(5, 5, 5)));
            Assert.IsTrue(box.Contains(new Vector(0, 0, 0))); // on boundary
            Assert.IsFalse(box.Contains(new Vector(11, 5, 5)));
        }

        [TestMethod]
        public void AABB_ClosestPoint_InsideReturnsSame()
        {
            var box = new AABB(new Vector(0, 0, 0), new Vector(10, 10, 10));
            var p = new Vector(5, 5, 5);
            var closest = box.ClosestPoint(p);

            Assert.AreEqual(5, closest.x, Tol);
            Assert.AreEqual(5, closest.y, Tol);
            Assert.AreEqual(5, closest.z, Tol);
        }

        [TestMethod]
        public void AABB_ClosestPoint_OutsideClampsToFace()
        {
            var box = new AABB(new Vector(0, 0, 0), new Vector(10, 10, 10));
            var p = new Vector(15, 5, 5);
            var closest = box.ClosestPoint(p);

            Assert.AreEqual(10, closest.x, Tol);
            Assert.AreEqual(5, closest.y, Tol);
            Assert.AreEqual(5, closest.z, Tol);
        }

        [TestMethod]
        public void AABB_Volume()
        {
            var box = new AABB(new Vector(0, 0, 0), new Vector(2, 3, 4));
            Assert.AreEqual(24, box.Volume, Tol);
        }

        #endregion

        #region BoundingSphere

        [TestMethod]
        public void BoundingSphere_Contains_PointInside()
        {
            var sphere = new BoundingSphere(new Vector(0, 0, 0), 5);

            Assert.IsTrue(sphere.Contains(new Vector(3, 0, 0)));
            Assert.IsTrue(sphere.Contains(new Vector(5, 0, 0))); // on surface
            Assert.IsFalse(sphere.Contains(new Vector(6, 0, 0)));
        }

        [TestMethod]
        public void BoundingSphere_Volume()
        {
            var sphere = new BoundingSphere(new Vector(0, 0, 0), 3);
            Assert.AreEqual(4.0 / 3.0 * Math.PI * 27, sphere.Volume, 1e-8);
        }

        #endregion

        #region Broad Phase — Overlap Tests

        [TestMethod]
        public void AABB_AABB_Intersects_Overlapping()
        {
            var a = new AABB(new Vector(0, 0, 0), new Vector(2, 2, 2));
            var b = new AABB(new Vector(1, 1, 1), new Vector(3, 3, 3));

            Assert.IsTrue(a.Intersects(b));
            Assert.IsTrue(b.Intersects(a));
        }

        [TestMethod]
        public void AABB_AABB_Intersects_Touching()
        {
            var a = new AABB(new Vector(0, 0, 0), new Vector(2, 2, 2));
            var b = new AABB(new Vector(2, 0, 0), new Vector(4, 2, 2));

            Assert.IsTrue(a.Intersects(b)); // touching faces count as intersection
        }

        [TestMethod]
        public void AABB_AABB_NoIntersect_Separated()
        {
            var a = new AABB(new Vector(0, 0, 0), new Vector(1, 1, 1));
            var b = new AABB(new Vector(5, 5, 5), new Vector(6, 6, 6));

            Assert.IsFalse(a.Intersects(b));
        }

        [TestMethod]
        public void AABB_AABB_NoIntersect_SingleAxisSeparation()
        {
            var a = new AABB(new Vector(0, 0, 0), new Vector(2, 2, 2));
            var b = new AABB(new Vector(0, 0, 3), new Vector(2, 2, 5));

            Assert.IsFalse(a.Intersects(b));
        }

        [TestMethod]
        public void Sphere_Sphere_Intersects_Overlapping()
        {
            var a = new BoundingSphere(new Vector(0, 0, 0), 2);
            var b = new BoundingSphere(new Vector(3, 0, 0), 2);

            Assert.IsTrue(a.Intersects(b));
        }

        [TestMethod]
        public void Sphere_Sphere_Intersects_Touching()
        {
            var a = new BoundingSphere(new Vector(0, 0, 0), 1);
            var b = new BoundingSphere(new Vector(2, 0, 0), 1);

            Assert.IsTrue(a.Intersects(b)); // exactly touching
        }

        [TestMethod]
        public void Sphere_Sphere_NoIntersect()
        {
            var a = new BoundingSphere(new Vector(0, 0, 0), 1);
            var b = new BoundingSphere(new Vector(5, 0, 0), 1);

            Assert.IsFalse(a.Intersects(b));
        }

        [TestMethod]
        public void AABB_Sphere_Intersects()
        {
            var box = new AABB(new Vector(0, 0, 0), new Vector(2, 2, 2));
            var sphere = new BoundingSphere(new Vector(3, 1, 1), 1.5);

            Assert.IsTrue(box.Intersects(sphere));
        }

        [TestMethod]
        public void AABB_Sphere_NoIntersect()
        {
            var box = new AABB(new Vector(0, 0, 0), new Vector(2, 2, 2));
            var sphere = new BoundingSphere(new Vector(5, 5, 5), 1);

            Assert.IsFalse(box.Intersects(sphere));
        }

        [TestMethod]
        public void Sphere_AABB_Symmetry()
        {
            var box = new AABB(new Vector(0, 0, 0), new Vector(2, 2, 2));
            var sphere = new BoundingSphere(new Vector(3, 1, 1), 1.5);

            Assert.AreEqual(box.Intersects(sphere), sphere.Intersects(box));
        }

        #endregion

        #region Narrow Phase — Contact Generation

        [TestMethod]
        public void SphereSphere_Contact_Normal()
        {
            var a = new BoundingSphere(new Vector(0, 0, 0), 2);
            var b = new BoundingSphere(new Vector(3, 0, 0), 2);

            var contact = a.SphereSphereContact(b);

            Assert.IsNotNull(contact);
            // Normal from A toward B: +X
            Assert.AreEqual(1, contact.Value.Normal.x, 1e-10);
            Assert.AreEqual(0, contact.Value.Normal.y, 1e-10);
            Assert.AreEqual(0, contact.Value.Normal.z, 1e-10);
        }

        [TestMethod]
        public void SphereSphere_Contact_PenetrationDepth()
        {
            var a = new BoundingSphere(new Vector(0, 0, 0), 2);
            var b = new BoundingSphere(new Vector(3, 0, 0), 2);

            var contact = a.SphereSphereContact(b);

            Assert.IsNotNull(contact);
            Assert.AreEqual(1, contact.Value.PenetrationDepth, 1e-10); // 2+2-3=1
        }

        [TestMethod]
        public void SphereSphere_Contact_Position()
        {
            var a = new BoundingSphere(new Vector(0, 0, 0), 2);
            var b = new BoundingSphere(new Vector(3, 0, 0), 2);

            var contact = a.SphereSphereContact(b);

            Assert.IsNotNull(contact);
            // Position on A's surface toward B: (2, 0, 0)
            Assert.AreEqual(2, contact.Value.Position.x, 1e-10);
            Assert.AreEqual(0, contact.Value.Position.y, 1e-10);
        }

        [TestMethod]
        public void SphereSphere_NoContact_Separated()
        {
            var a = new BoundingSphere(new Vector(0, 0, 0), 1);
            var b = new BoundingSphere(new Vector(5, 0, 0), 1);

            Assert.IsNull(a.SphereSphereContact(b));
        }

        [TestMethod]
        public void SphereAABB_Contact_SphereTouchingFace()
        {
            var sphere = new BoundingSphere(new Vector(4, 1, 1), 1.5);
            var box = new AABB(new Vector(0, 0, 0), new Vector(3, 3, 3));

            var contact = sphere.SphereAABBContact(box);

            Assert.IsNotNull(contact);
            Assert.IsTrue(contact.Value.PenetrationDepth > 0);
            // Normal points from box toward sphere: +X
            Assert.IsTrue(contact.Value.Normal.x > 0);
        }

        [TestMethod]
        public void SphereAABB_NoContact_Separated()
        {
            var sphere = new BoundingSphere(new Vector(10, 0, 0), 1);
            var box = new AABB(new Vector(0, 0, 0), new Vector(2, 2, 2));

            Assert.IsNull(sphere.SphereAABBContact(box));
        }

        #endregion

        #region Collision Response

        [TestMethod]
        public void ResolveCollision_EqualMasses_Elastic_SwapVelocities()
        {
            var a = RigidBody.CreateSolidSphere(1, 1);
            a.Position = new Vector(0, 0, 0);
            a.Velocity = new Vector(5, 0, 0);

            var b = RigidBody.CreateSolidSphere(1, 1);
            b.Position = new Vector(2, 0, 0);
            b.Velocity = new Vector(0, 0, 0);

            var contact = new ContactPoint(
                new Vector(1, 0, 0), // midpoint
                new Vector(1, 0, 0), // normal from A to B
                0.0);

            CollisionResponse.ResolveCollision(ref a, ref b, contact, restitution: 1.0);

            // Equal mass elastic head-on: velocities swap
            Assert.AreEqual(0, a.Velocity.x, 1e-10);
            Assert.AreEqual(5, b.Velocity.x, 1e-10);
        }

        [TestMethod]
        public void ResolveCollision_PerfectlyInelastic_StickTogether()
        {
            var a = RigidBody.CreateSolidSphere(2, 1);
            a.Position = new Vector(0, 0, 0);
            a.Velocity = new Vector(10, 0, 0);

            var b = RigidBody.CreateSolidSphere(2, 1);
            b.Position = new Vector(2, 0, 0);
            b.Velocity = new Vector(0, 0, 0);

            var contact = new ContactPoint(
                new Vector(1, 0, 0),
                new Vector(1, 0, 0),
                0.0);

            CollisionResponse.ResolveCollision(ref a, ref b, contact, restitution: 0.0);

            // e=0: final relative velocity along normal = 0
            Assert.AreEqual(a.Velocity.x, b.Velocity.x, 1e-10);
        }

        [TestMethod]
        public void ResolveCollision_MomentumConserved()
        {
            var a = RigidBody.CreateSolidSphere(3, 1);
            a.Position = new Vector(0, 0, 0);
            a.Velocity = new Vector(8, 0, 0);

            var b = RigidBody.CreateSolidSphere(5, 1);
            b.Position = new Vector(2, 0, 0);
            b.Velocity = new Vector(-2, 0, 0);

            var pBefore = a.Mass * a.Velocity + b.Mass * b.Velocity;

            var contact = new ContactPoint(
                new Vector(1, 0, 0),
                new Vector(1, 0, 0),
                0.0);

            CollisionResponse.ResolveCollision(ref a, ref b, contact, restitution: 0.7);

            var pAfter = a.Mass * a.Velocity + b.Mass * b.Velocity;

            Assert.AreEqual(pBefore.x, pAfter.x, 1e-10);
            Assert.AreEqual(pBefore.y, pAfter.y, 1e-10);
            Assert.AreEqual(pBefore.z, pAfter.z, 1e-10);
        }

        [TestMethod]
        public void ResolveCollision_EnergyConserved_WhenElastic()
        {
            var a = RigidBody.CreateSolidSphere(3, 1);
            a.Position = new Vector(0, 0, 0);
            a.Velocity = new Vector(8, 2, 0);

            var b = RigidBody.CreateSolidSphere(5, 1);
            b.Position = new Vector(2, 0, 0);
            b.Velocity = new Vector(-2, 1, 0);

            double keBefore = 0.5 * a.Mass * a.Velocity.Dot(a.Velocity)
                            + 0.5 * b.Mass * b.Velocity.Dot(b.Velocity);

            var contact = new ContactPoint(
                new Vector(1, 0, 0),
                new Vector(1, 0, 0),
                0.0);

            CollisionResponse.ResolveCollision(ref a, ref b, contact, restitution: 1.0);

            double keAfter = 0.5 * a.Mass * a.Velocity.Dot(a.Velocity)
                           + 0.5 * b.Mass * b.Velocity.Dot(b.Velocity);

            Assert.AreEqual(keBefore, keAfter, 1e-8);
        }

        [TestMethod]
        public void ResolveCollision_StaticBody_OnlyDynamicMoves()
        {
            var wall = RigidBody.CreateStatic(new Vector(5, 0, 0));

            var ball = RigidBody.CreateSolidSphere(1, 1);
            ball.Position = new Vector(3, 0, 0);
            ball.Velocity = new Vector(10, 0, 0);

            var contact = new ContactPoint(
                new Vector(5, 0, 0),
                new Vector(-1, 0, 0), // normal from ball toward wall → wall pushes back
                0.0);

            // Normal from A(wall) toward B(ball) should be (-1,0,0)
            // Actually: let's set up as A=ball, B=wall, normal from ball toward wall = (1,0,0)
            contact = new ContactPoint(
                new Vector(4, 0, 0),
                new Vector(1, 0, 0),
                0.0);

            CollisionResponse.ResolveCollision(ref ball, ref wall, contact, restitution: 1.0);

            Assert.AreEqual(5, wall.Position.x, Tol); // wall didn't move
            Assert.IsTrue(ball.Velocity.x < 0); // ball bounced back
        }

        [TestMethod]
        public void ResolveCollision_SeparatingBodies_NoChange()
        {
            var a = RigidBody.CreateSolidSphere(1, 1);
            a.Velocity = new Vector(-5, 0, 0); // moving away

            var b = RigidBody.CreateSolidSphere(1, 1);
            b.Position = new Vector(2, 0, 0);
            b.Velocity = new Vector(5, 0, 0); // moving away

            var contact = new ContactPoint(
                new Vector(1, 0, 0),
                new Vector(1, 0, 0),
                0.0);

            CollisionResponse.ResolveCollision(ref a, ref b, contact, restitution: 1.0);

            // Already separating — no impulse applied
            Assert.AreEqual(-5, a.Velocity.x, Tol);
            Assert.AreEqual(5, b.Velocity.x, Tol);
        }

        #endregion

        #region Positional Correction

        [TestMethod]
        public void CorrectPositions_PushesApart()
        {
            var a = RigidBody.CreateSolidSphere(1, 1);
            a.Position = new Vector(0, 0, 0);

            var b = RigidBody.CreateSolidSphere(1, 1);
            b.Position = new Vector(1.5, 0, 0);

            var contact = new ContactPoint(
                new Vector(0.75, 0, 0),
                new Vector(1, 0, 0),
                0.5); // 0.5m penetration

            double distBefore = (b.Position - a.Position).GetMagnitude();

            CollisionResponse.CorrectPositions(ref a, ref b, contact);

            double distAfter = (b.Position - a.Position).GetMagnitude();

            Assert.IsTrue(distAfter > distBefore);
        }

        [TestMethod]
        public void CorrectPositions_StaticBody_OnlyDynamicMoves()
        {
            var wall = RigidBody.CreateStatic(new Vector(5, 0, 0));

            var ball = RigidBody.CreateSolidSphere(1, 1);
            ball.Position = new Vector(4.5, 0, 0);

            var contact = new ContactPoint(
                new Vector(4.75, 0, 0),
                new Vector(1, 0, 0),
                0.5);

            CollisionResponse.CorrectPositions(ref ball, ref wall, contact);

            Assert.AreEqual(5, wall.Position.x, Tol); // wall didn't move
        }

        [TestMethod]
        public void CorrectPositions_BelowSlop_NoCorrection()
        {
            var a = RigidBody.CreateSolidSphere(1, 1);
            a.Position = new Vector(0, 0, 0);

            var b = RigidBody.CreateSolidSphere(1, 1);
            b.Position = new Vector(2, 0, 0);

            var contact = new ContactPoint(
                new Vector(1, 0, 0),
                new Vector(1, 0, 0),
                0.005); // below default slop of 0.01

            var posA = a.Position;
            var posB = b.Position;

            CollisionResponse.CorrectPositions(ref a, ref b, contact);

            Assert.AreEqual(posA.x, a.Position.x, Tol);
            Assert.AreEqual(posB.x, b.Position.x, Tol);
        }

        [TestMethod]
        public void CorrectPositions_EqualMass_SymmetricCorrection()
        {
            var a = RigidBody.CreateSolidSphere(5, 1);
            a.Position = new Vector(0, 0, 0);

            var b = RigidBody.CreateSolidSphere(5, 1);
            b.Position = new Vector(1, 0, 0);

            var contact = new ContactPoint(
                new Vector(0.5, 0, 0),
                new Vector(1, 0, 0),
                1.0);

            CollisionResponse.CorrectPositions(ref a, ref b, contact);

            // Equal mass: both move the same amount in opposite directions
            double moveA = Math.Abs(a.Position.x);
            double moveB = Math.Abs(b.Position.x - 1);
            Assert.AreEqual(moveA, moveB, 1e-10);
        }

        #endregion

        #region Integration — Full Collision Scenario

        [TestMethod]
        public void FullScenario_TwoSpheres_BounceApart()
        {
            // Set up two spheres approaching each other
            var a = RigidBody.CreateSolidSphere(1, 1);
            a.Position = new Vector(0, 0, 0);
            a.Velocity = new Vector(5, 0, 0);

            var b = RigidBody.CreateSolidSphere(1, 1);
            b.Position = new Vector(3, 0, 0);
            b.Velocity = new Vector(-3, 0, 0);

            var sA = new BoundingSphere(a.Position, 1);
            var sB = new BoundingSphere(b.Position, 1);

            // Detect
            var contact = sA.SphereSphereContact(sB);
            Assert.IsNotNull(contact);

            // Momentum before
            double pxBefore = a.Mass * a.Velocity.x + b.Mass * b.Velocity.x;

            // Resolve
            CollisionResponse.ResolveCollision(ref a, ref b, contact.Value, restitution: 1.0);

            // Momentum after
            double pxAfter = a.Mass * a.Velocity.x + b.Mass * b.Velocity.x;
            Assert.AreEqual(pxBefore, pxAfter, 1e-10);

            // Bodies should be separating now
            double relVel = b.Velocity.x - a.Velocity.x;
            Assert.IsTrue(relVel > 0); // moving apart
        }

        #endregion
    }
}
