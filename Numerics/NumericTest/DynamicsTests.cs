using CSharpNumerics.Physics;
using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Physics.Objects;
using CSharpNumerics.Numerics.Objects;

namespace NumericTest
{
    [TestClass]
    public class DynamicsTests
    {
        #region Newton's Laws

        [TestMethod]
        public void Acceleration_FromForceAndMass()
        {
            var force = new Vector(10, 0, 0);
            var a = force.Acceleration(mass: 5);

            Assert.AreEqual(2, a.x, 1e-10);
            Assert.AreEqual(0, a.y, 1e-10);
            Assert.AreEqual(0, a.z, 1e-10);
        }

        [TestMethod]
        public void Force_FromMassAndAcceleration()
        {
            var a = new Vector(0, 0, -9.8);
            var f = 10.0.Force(a);

            Assert.AreEqual(0, f.x, 1e-10);
            Assert.AreEqual(0, f.y, 1e-10);
            Assert.AreEqual(-98, f.z, 1e-10);
        }

        [TestMethod]
        public void NetForce_SumsAllForces()
        {
            var f1 = new Vector(10, 0, 0);
            var f2 = new Vector(0, 5, 0);
            var f3 = new Vector(-3, 0, 7);

            var net = f1.NetForce(f2, f3);

            Assert.AreEqual(7, net.x, 1e-10);
            Assert.AreEqual(5, net.y, 1e-10);
            Assert.AreEqual(7, net.z, 1e-10);
        }

        [TestMethod]
        public void Weight_DefaultDownward()
        {
            double mass = 10;
            var w = mass.Weight();

            Assert.AreEqual(0, w.x, 1e-10);
            Assert.AreEqual(0, w.y, 1e-10);
            Assert.AreEqual(-mass * PhysicsConstants.GravitationalAcceleration, w.z, 1e-10);
        }

        [TestMethod]
        public void Weight_CustomDirection()
        {
            double mass = 5;
            var w = mass.Weight(new Vector(0, -1, 0));

            Assert.AreEqual(0, w.x, 1e-10);
            Assert.AreEqual(-mass * PhysicsConstants.GravitationalAcceleration, w.y, 1e-10);
            Assert.AreEqual(0, w.z, 1e-10);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Acceleration_ZeroMass_Throws()
        {
            new Vector(1, 0, 0).Acceleration(0);
        }

        #endregion

        #region Momentum & Impulse

        [TestMethod]
        public void Momentum_Vector()
        {
            var p = 5.0.Momentum(new Vector(3, 4, 0));

            Assert.AreEqual(15, p.x, 1e-10);
            Assert.AreEqual(20, p.y, 1e-10);
        }

        [TestMethod]
        public void Momentum_Scalar()
        {
            double p = 5.0.Momentum(3.0);
            Assert.AreEqual(15, p, 1e-10);
        }

        [TestMethod]
        public void Impulse_ForceTimesDuration()
        {
            var force = new Vector(100, 0, 0);
            var J = force.Impulse(0.5);

            Assert.AreEqual(50, J.x, 1e-10);
        }

        [TestMethod]
        public void ImpulseFromVelocityChange_MatchesForceImpulse()
        {
            double mass = 2;
            var vBefore = new Vector(5, 0, 0);
            var vAfter = new Vector(15, 0, 0);

            var J = mass.ImpulseFromVelocityChange(vBefore, vAfter);

            Assert.AreEqual(20, J.x, 1e-10); // 2 * (15-5) = 20
        }

        [TestMethod]
        public void ApplyImpulse_UpdatesVelocity()
        {
            var impulse = new Vector(10, 0, 0);
            var v0 = new Vector(0, 0, 0);

            var v = impulse.ApplyImpulse(mass: 5, v0);

            Assert.AreEqual(2, v.x, 1e-10); // 10/5 = 2
        }

        [TestMethod]
        public void Momentum_Conservation_InelasticCollision()
        {
            double m1 = 3, m2 = 2;
            var v1 = new Vector(10, 0, 0);
            var v2 = new Vector(-5, 0, 0);

            var pBefore = m1 * v1 + m2 * v2;
            var vFinal = m1.InelasticCollisionVelocity(v1, m2, v2);
            var pAfter = (m1 + m2) * vFinal;

            Assert.AreEqual(pBefore.x, pAfter.x, 1e-10);
            Assert.AreEqual(pBefore.y, pAfter.y, 1e-10);
            Assert.AreEqual(pBefore.z, pAfter.z, 1e-10);
        }

        #endregion

        #region Energy

        [TestMethod]
        public void KineticEnergy_Vector()
        {
            // KE = 0.5 * 4 * (3² + 4²) = 0.5 * 4 * 25 = 50
            double ke = 4.0.KineticEnergy(new Vector(3, 4, 0));
            Assert.AreEqual(50, ke, 1e-10);
        }

        [TestMethod]
        public void KineticEnergy_Scalar()
        {
            double ke = 4.0.KineticEnergy(5.0);
            Assert.AreEqual(50, ke, 1e-10);
        }

        [TestMethod]
        public void PotentialEnergy_Surface()
        {
            // PE = 10 * 9.80665 * 5 = 490.3325
            double pe = 10.0.PotentialEnergy(5);
            Assert.AreEqual(10 * PhysicsConstants.GravitationalAcceleration * 5, pe, 1e-10);
        }

        [TestMethod]
        public void GravitationalPotentialEnergy_TwoMasses()
        {
            // U = -G * m1 * m2 / r
            double U = PhysicsConstants.EarthMass.GravitationalPotentialEnergy(PhysicsConstants.MoonMass, 3.844e8);
            Assert.IsTrue(U < 0); // always negative
            double expected = -PhysicsConstants.GravitationalConstant
                * PhysicsConstants.EarthMass * PhysicsConstants.MoonMass / 3.844e8;
            Assert.AreEqual(expected, U, Math.Abs(expected) * 1e-10);
        }

        [TestMethod]
        public void MechanicalEnergy_SumOfKEAndPE()
        {
            double mass = 2;
            var velocity = new Vector(3, 4, 0);
            double height = 10;

            double E = mass.MechanicalEnergy(velocity, height);
            double ke = mass.KineticEnergy(velocity);
            double pe = mass.PotentialEnergy(height);

            Assert.AreEqual(ke + pe, E, 1e-10);
        }

        [TestMethod]
        public void SpeedFromKineticEnergy_RoundTrip()
        {
            double mass = 5;
            double speed = 12;
            double ke = mass.KineticEnergy(speed);
            double recovered = mass.SpeedFromKineticEnergy(ke);

            Assert.AreEqual(speed, recovered, 1e-10);
        }

        [TestMethod]
        public void FreeFall_EnergyConservation()
        {
            // Drop from 20m: PE at top = KE at bottom
            double mass = 3;
            double height = 20;

            double peTop = mass.PotentialEnergy(height);
            double vBottom = height.FreeFallVelocity();
            double keBottom = mass.KineticEnergy(vBottom);

            Assert.AreEqual(peTop, keBottom, 1e-6);
        }

        #endregion

        #region Work & Power

        [TestMethod]
        public void Work_ParallelForce()
        {
            var force = new Vector(10, 0, 0);
            var displacement = new Vector(5, 0, 0);

            double W = force.Work(displacement);
            Assert.AreEqual(50, W, 1e-10);
        }

        [TestMethod]
        public void Work_PerpendicularForce_IsZero()
        {
            var force = new Vector(10, 0, 0);
            var displacement = new Vector(0, 5, 0);

            double W = force.Work(displacement);
            Assert.AreEqual(0, W, 1e-10);
        }

        [TestMethod]
        public void Work_Scalar_WithAngle()
        {
            // W = F * d * cos(60°) = 10 * 5 * 0.5 = 25
            double W = 10.0.Work(5, Math.PI / 3);
            Assert.AreEqual(25, W, 1e-10);
        }

        [TestMethod]
        public void Power_Instantaneous()
        {
            var force = new Vector(10, 0, 0);
            var velocity = new Vector(3, 0, 0);

            double P = force.Power(velocity);
            Assert.AreEqual(30, P, 1e-10);
        }

        [TestMethod]
        public void AveragePower_WorkOverTime()
        {
            double P = 100.0.AveragePower(duration: 5);
            Assert.AreEqual(20, P, 1e-10);
        }

        [TestMethod]
        public void WorkEnergyTheorem_EqualsKEChange()
        {
            double mass = 4;
            var v1 = new Vector(2, 0, 0);
            var v2 = new Vector(5, 0, 0);

            double deltaKE = mass.WorkEnergyTheorem(v1, v2);
            double expected = mass.KineticEnergy(v2) - mass.KineticEnergy(v1);

            Assert.AreEqual(expected, deltaKE, 1e-10);
        }

        [TestMethod]
        public void WorkEnergyTheorem_MatchesWorkDone()
        {
            // Constant force F=10N on 2kg object over 5m from rest
            double mass = 2;
            var force = new Vector(10, 0, 0);
            var displacement = new Vector(5, 0, 0);

            double W = force.Work(displacement);

            // v² = 2*a*s = 2*(10/2)*5 = 50 → v = √50
            var vAfter = new Vector(Math.Sqrt(50), 0, 0);
            double deltaKE = mass.WorkEnergyTheorem(new Vector(0, 0, 0), vAfter);

            Assert.AreEqual(W, deltaKE, 1e-10);
        }

        #endregion

        #region Collisions

        [TestMethod]
        public void ElasticCollision_MomentumConserved()
        {
            double m1 = 3, m2 = 5;
            double v1 = 10, v2 = -2;

            var (v1f, v2f) = m1.ElasticCollision(v1, m2, v2);

            double pBefore = m1 * v1 + m2 * v2;
            double pAfter = m1 * v1f + m2 * v2f;
            Assert.AreEqual(pBefore, pAfter, 1e-10);
        }

        [TestMethod]
        public void ElasticCollision_EnergyConserved()
        {
            double m1 = 3, m2 = 5;
            double v1 = 10, v2 = -2;

            var (v1f, v2f) = m1.ElasticCollision(v1, m2, v2);

            double keBefore = 0.5 * m1 * v1 * v1 + 0.5 * m2 * v2 * v2;
            double keAfter = 0.5 * m1 * v1f * v1f + 0.5 * m2 * v2f * v2f;
            Assert.AreEqual(keBefore, keAfter, 1e-10);
        }

        [TestMethod]
        public void ElasticCollision_EqualMasses_SwapVelocities()
        {
            double m = 5;
            double v1 = 8, v2 = 0;

            var (v1f, v2f) = m.ElasticCollision(v1, m, v2);

            Assert.AreEqual(0, v1f, 1e-10);
            Assert.AreEqual(8, v2f, 1e-10);
        }

        [TestMethod]
        public void InelasticCollision_MomentumConserved()
        {
            double m1 = 3, m2 = 2;
            double v1 = 10, v2 = -5;

            double vf = m1.InelasticCollisionVelocity(v1, m2, v2);
            double pBefore = m1 * v1 + m2 * v2;
            double pAfter = (m1 + m2) * vf;

            Assert.AreEqual(pBefore, pAfter, 1e-10);
        }

        [TestMethod]
        public void InelasticCollision_EnergyLost()
        {
            double m1 = 3, m2 = 2;
            double v1 = 10, v2 = -5;

            double loss = m1.InelasticCollisionEnergyLoss(v1, m2, v2);

            Assert.IsTrue(loss > 0); // energy is always lost
        }

        [TestMethod]
        public void InelasticCollision_Vector_3D()
        {
            double m1 = 2, m2 = 3;
            var v1 = new Vector(6, 0, 0);
            var v2 = new Vector(0, 4, 0);

            var vf = m1.InelasticCollisionVelocity(v1, m2, v2);

            // vf = (2*(6,0,0) + 3*(0,4,0)) / 5 = (12/5, 12/5, 0) = (2.4, 2.4, 0)
            Assert.AreEqual(2.4, vf.x, 1e-10);
            Assert.AreEqual(2.4, vf.y, 1e-10);
            Assert.AreEqual(0, vf.z, 1e-10);
        }

        [TestMethod]
        public void CoefficientOfRestitution_Elastic()
        {
            // Elastic: e = 1
            double m = 5;
            double v1 = 8, v2 = 0;
            var (v1f, v2f) = m.ElasticCollision(v1, m, v2);

            double e = v1.CoefficientOfRestitution(v2, v1f, v2f);
            Assert.AreEqual(1, e, 1e-10);
        }

        [TestMethod]
        public void CoefficientOfRestitution_PerfectlyInelastic()
        {
            double m1 = 3, m2 = 2;
            double v1 = 10, v2 = 0;
            double vf = m1.InelasticCollisionVelocity(v1, m2, v2);

            double e = v1.CoefficientOfRestitution(v2, vf, vf); // both stick together
            Assert.AreEqual(0, e, 1e-10);
        }

        #endregion

        #region RigidBody

        [TestMethod]
        public void RigidBody_CreateSolidSphere_CorrectInertia()
        {
            double mass = 10, radius = 2;
            var body = RigidBody.CreateSolidSphere(mass, radius);

            double expected = 2.0 / 5.0 * mass * radius * radius; // 16
            Assert.AreEqual(expected, body.InertiaTensorBody.values[0, 0], 1e-10);
            Assert.AreEqual(expected, body.InertiaTensorBody.values[1, 1], 1e-10);
            Assert.AreEqual(expected, body.InertiaTensorBody.values[2, 2], 1e-10);
            Assert.AreEqual(0, body.InertiaTensorBody.values[0, 1], 1e-10);
        }

        [TestMethod]
        public void RigidBody_CreateSolidBox_CorrectInertia()
        {
            double mass = 12, w = 2, h = 3, d = 4;
            var body = RigidBody.CreateSolidBox(mass, w, h, d);

            double Ix = mass / 12.0 * (h * h + d * d); // 12/12*(9+16) = 25
            double Iy = mass / 12.0 * (w * w + d * d); // 12/12*(4+16) = 20
            double Iz = mass / 12.0 * (w * w + h * h); // 12/12*(4+9) = 13

            Assert.AreEqual(Ix, body.InertiaTensorBody.values[0, 0], 1e-10);
            Assert.AreEqual(Iy, body.InertiaTensorBody.values[1, 1], 1e-10);
            Assert.AreEqual(Iz, body.InertiaTensorBody.values[2, 2], 1e-10);
        }

        [TestMethod]
        public void RigidBody_CreateSolidCylinder_CorrectInertia()
        {
            double mass = 6, radius = 1, height = 4;
            var body = RigidBody.CreateSolidCylinder(mass, radius, height);

            double Iaxis = 0.5 * mass * radius * radius; // 3
            double Idiameter = mass / 12.0 * (3 * radius * radius + height * height); // 6/12*(3+16) = 9.5

            Assert.AreEqual(Idiameter, body.InertiaTensorBody.values[0, 0], 1e-10);
            Assert.AreEqual(Idiameter, body.InertiaTensorBody.values[1, 1], 1e-10);
            Assert.AreEqual(Iaxis, body.InertiaTensorBody.values[2, 2], 1e-10);
        }

        [TestMethod]
        public void RigidBody_CreateHollowSphere_GreaterThanSolid()
        {
            double mass = 5, radius = 3;
            var solid = RigidBody.CreateSolidSphere(mass, radius);
            var hollow = RigidBody.CreateHollowSphere(mass, radius);

            Assert.IsTrue(hollow.InertiaTensorBody.values[0, 0] >
                          solid.InertiaTensorBody.values[0, 0]);
        }

        [TestMethod]
        public void RigidBody_CreateStatic_ZeroInverseMass()
        {
            var body = RigidBody.CreateStatic(new Vector(0, 5, 0));

            Assert.IsTrue(body.IsStatic);
            Assert.AreEqual(0, body.InverseMass, 1e-10);
            Assert.AreEqual(5, body.Position.y, 1e-10);
        }

        [TestMethod]
        public void RigidBody_ApplyForce_AccumulatesAndAccelerates()
        {
            var body = RigidBody.CreateSolidSphere(10, 1);
            body.ApplyForce(new Vector(20, 0, 0));
            body.ApplyForce(new Vector(0, 10, 0));

            Assert.AreEqual(20, body.AccumulatedForce.x, 1e-10);
            Assert.AreEqual(10, body.AccumulatedForce.y, 1e-10);

            var a = body.LinearAcceleration;
            Assert.AreEqual(2, a.x, 1e-10);
            Assert.AreEqual(1, a.y, 1e-10);
        }

        [TestMethod]
        public void RigidBody_ApplyForceAtPoint_GeneratesTorque()
        {
            var body = RigidBody.CreateSolidSphere(10, 1);
            body.Position = new Vector(0, 0, 0);

            // Force along +X at point (0, 1, 0) → torque along +Z
            body.ApplyForceAtPoint(new Vector(10, 0, 0), new Vector(0, 1, 0));

            // τ = r × F = (0,1,0) × (10,0,0) = (0,0,-10)
            Assert.AreEqual(0, body.AccumulatedTorque.x, 1e-10);
            Assert.AreEqual(0, body.AccumulatedTorque.y, 1e-10);
            Assert.AreEqual(-10, body.AccumulatedTorque.z, 1e-10);
        }

        [TestMethod]
        public void RigidBody_ClearForces_ResetsAccumulators()
        {
            var body = RigidBody.CreateSolidSphere(5, 1);
            body.ApplyForce(new Vector(100, 200, 300));
            body.ApplyTorque(new Vector(1, 2, 3));
            body.ClearForces();

            Assert.AreEqual(0, body.AccumulatedForce.GetMagnitude(), 1e-10);
            Assert.AreEqual(0, body.AccumulatedTorque.GetMagnitude(), 1e-10);
        }

        [TestMethod]
        public void RigidBody_KineticEnergy_TranslationalAndRotational()
        {
            var body = RigidBody.CreateSolidSphere(10, 2);
            body.Velocity = new Vector(3, 0, 0);
            body.AngularVelocity = new Vector(0, 0, 5);

            double keTranslational = 0.5 * 10 * 9; // 45
            double I = 2.0 / 5.0 * 10 * 4; // 16
            double keRotational = 0.5 * I * 25; // 200

            Assert.AreEqual(keTranslational + keRotational, body.KineticEnergy, 1e-10);
        }

        [TestMethod]
        public void RigidBody_LinearMomentum_EqualsM_Times_V()
        {
            var body = RigidBody.CreateSolidBox(5, 1, 1, 1);
            body.Velocity = new Vector(4, 0, 0);

            var p = body.LinearMomentum;
            Assert.AreEqual(20, p.x, 1e-10);
        }

        [TestMethod]
        public void RigidBody_AngularAcceleration_CorrectForSphere()
        {
            var body = RigidBody.CreateSolidSphere(10, 1);
            // I = 2/5*10*1 = 4, I_inv = 1/4 = 0.25
            body.ApplyTorque(new Vector(0, 0, 8));

            var alpha = body.AngularAcceleration;
            Assert.AreEqual(0, alpha.x, 1e-10);
            Assert.AreEqual(0, alpha.y, 1e-10);
            Assert.AreEqual(2, alpha.z, 1e-10); // 8 * 0.25 = 2
        }

        [TestMethod]
        public void RigidBody_Static_NoAcceleration()
        {
            var body = RigidBody.CreateStatic(new Vector(0, 0, 0));
            body.ApplyForce(new Vector(1000, 0, 0));
            body.ApplyTorque(new Vector(0, 0, 500));

            var a = body.LinearAcceleration;
            Assert.AreEqual(0, a.GetMagnitude(), 1e-10);
        }

        #endregion

        #region Moment of Inertia (Scalar)

        [TestMethod]
        public void MomentOfInertia_SolidSphere()
        {
            double I = 10.0.MomentOfInertiaSolidSphere(2);
            Assert.AreEqual(2.0 / 5.0 * 10 * 4, I, 1e-10);
        }

        [TestMethod]
        public void MomentOfInertia_HollowGreaterThanSolid()
        {
            double solid = 5.0.MomentOfInertiaSolidSphere(3);
            double hollow = 5.0.MomentOfInertiaHollowSphere(3);
            Assert.IsTrue(hollow > solid);
        }

        [TestMethod]
        public void MomentOfInertia_ThinRodEnd_FourTimesCenter()
        {
            // I_end = mL²/3, I_center = mL²/12 → I_end = 4·I_center
            // Also: I_end = I_center + m*(L/2)² via parallel axis
            double mass = 6, length = 4;
            double Icenter = mass.MomentOfInertiaThinRod(length);
            double Iend = mass.MomentOfInertiaThinRodEnd(length);

            Assert.AreEqual(4 * Icenter, Iend, 1e-10);
        }

        [TestMethod]
        public void ParallelAxis_Scalar_ThinRod()
        {
            double mass = 6, length = 4;
            double Icm = mass.MomentOfInertiaThinRod(length);
            double Iend = Icm.ParallelAxis(mass, length / 2.0);

            Assert.AreEqual(mass.MomentOfInertiaThinRodEnd(length), Iend, 1e-10);
        }

        [TestMethod]
        public void MomentOfInertia_SolidCylinder()
        {
            double I = 8.0.MomentOfInertiaSolidCylinder(3);
            Assert.AreEqual(0.5 * 8 * 9, I, 1e-10);
        }

        [TestMethod]
        public void MomentOfInertia_SolidBox()
        {
            double I = 12.0.MomentOfInertiaSolidBox(3, 4);
            Assert.AreEqual(12.0 / 12.0 * (9 + 16), I, 1e-10);
        }

        #endregion

        #region Inertia Tensor (3×3)

        [TestMethod]
        public void InertiaTensor_SolidSphere_IsotropicDiagonal()
        {
            var I = 10.0.InertiaTensorSolidSphere(2);
            double expected = 2.0 / 5.0 * 10 * 4;

            Assert.AreEqual(expected, I.values[0, 0], 1e-10);
            Assert.AreEqual(expected, I.values[1, 1], 1e-10);
            Assert.AreEqual(expected, I.values[2, 2], 1e-10);
            Assert.AreEqual(0, I.values[0, 1], 1e-10);
        }

        [TestMethod]
        public void InertiaTensor_SolidBox_DifferentAxes()
        {
            var I = 12.0.InertiaTensorSolidBox(1, 2, 3);

            // Ix = m/12*(h²+d²) = 12/12*(4+9) = 13
            // Iy = m/12*(w²+d²) = 12/12*(1+9) = 10
            // Iz = m/12*(w²+h²) = 12/12*(1+4) = 5
            Assert.AreEqual(13, I.values[0, 0], 1e-10);
            Assert.AreEqual(10, I.values[1, 1], 1e-10);
            Assert.AreEqual(5, I.values[2, 2], 1e-10);
        }

        [TestMethod]
        public void InertiaTensor_SolidCylinder_AxisLessThanDiameter()
        {
            var I = 6.0.InertiaTensorSolidCylinder(2, 4);

            // Iz (axis) = ½mr² = 12
            // Ix = Iy (diameter) = m/12*(3r²+h²) = 6/12*(12+16) = 14
            Assert.AreEqual(14, I.values[0, 0], 1e-10);
            Assert.AreEqual(14, I.values[1, 1], 1e-10);
            Assert.AreEqual(12, I.values[2, 2], 1e-10);
        }

        [TestMethod]
        public void ParallelAxis_Tensor_SphereOffset()
        {
            double mass = 5, radius = 1;
            var Icm = mass.InertiaTensorSolidSphere(radius);
            var offset = new Vector(3, 0, 0);
            var Inew = Icm.ParallelAxis(mass, offset);

            // Ixx should stay the same (axis through offset direction)
            // Iyy, Izz should increase by m*d² = 5*9 = 45
            double Icm_val = 2.0 / 5.0 * mass * radius * radius;
            Assert.AreEqual(Icm_val, Inew.values[0, 0], 1e-10);
            Assert.AreEqual(Icm_val + 45, Inew.values[1, 1], 1e-10);
            Assert.AreEqual(Icm_val + 45, Inew.values[2, 2], 1e-10);

            // Off-diagonal: -m*dx*dy = 0 since offset is (3,0,0)
            Assert.AreEqual(0, Inew.values[0, 1], 1e-10);
        }

        #endregion

        #region Torque & Rotational Dynamics

        [TestMethod]
        public void Torque_Vector_CrossProduct()
        {
            var r = new Vector(0, 1, 0);
            var F = new Vector(10, 0, 0);

            var tau = r.Torque(F);

            // (0,1,0) × (10,0,0) = (0,0,-10)
            Assert.AreEqual(0, tau.x, 1e-10);
            Assert.AreEqual(0, tau.y, 1e-10);
            Assert.AreEqual(-10, tau.z, 1e-10);
        }

        [TestMethod]
        public void Torque_Scalar_PerpendicularForce()
        {
            double tau = 20.0.Torque(momentArm: 3);
            Assert.AreEqual(60, tau, 1e-10);
        }

        [TestMethod]
        public void Torque_Scalar_WithAngle()
        {
            // τ = F·r·sin(30°) = 20·3·0.5 = 30
            double tau = 20.0.Torque(3, Math.PI / 6);
            Assert.AreEqual(30, tau, 1e-10);
        }

        [TestMethod]
        public void AngularMomentum_Matrix_TimesOmega()
        {
            var I = 10.0.InertiaTensorSolidSphere(1); // diag(4, 4, 4)
            var omega = new Vector(0, 0, 5);

            var L = I.AngularMomentum(omega);

            Assert.AreEqual(0, L.x, 1e-10);
            Assert.AreEqual(0, L.y, 1e-10);
            Assert.AreEqual(4 * 5, L.z, 1e-10);
        }

        [TestMethod]
        public void AngularMomentum_Scalar()
        {
            double I = 10.0.MomentOfInertiaSolidSphere(1);
            double L = I.AngularMomentum(5);
            Assert.AreEqual(I * 5, L, 1e-10);
        }

        [TestMethod]
        public void AngularAcceleration_FromInverseInertia()
        {
            var I = 10.0.InertiaTensorSolidSphere(1);
            var Iinv = I.Inverse();
            var torque = new Vector(0, 0, 8);

            var alpha = Iinv.AngularAcceleration(torque);

            Assert.AreEqual(2, alpha.z, 1e-10); // 8/4 = 2
        }

        [TestMethod]
        public void RotationalKineticEnergy_Matrix()
        {
            var I = 10.0.InertiaTensorSolidSphere(2); // diag(16,16,16)
            var omega = new Vector(0, 0, 3);

            double KE = I.RotationalKineticEnergy(omega);
            Assert.AreEqual(0.5 * 16 * 9, KE, 1e-10); // 72
        }

        [TestMethod]
        public void RotationalKineticEnergy_Scalar()
        {
            double I = 16;
            double KE = I.RotationalKineticEnergy(3);
            Assert.AreEqual(0.5 * 16 * 9, KE, 1e-10);
        }

        [TestMethod]
        public void RigidBody_AngularMomentum_MatchesExtension()
        {
            var body = RigidBody.CreateSolidSphere(10, 2);
            body.AngularVelocity = new Vector(1, 2, 3);

            var LBody = body.AngularMomentum;
            var LExtension = body.InertiaTensorBody.AngularMomentum(body.AngularVelocity);

            Assert.AreEqual(LExtension.x, LBody.x, 1e-10);
            Assert.AreEqual(LExtension.y, LBody.y, 1e-10);
           Assert.AreEqual(LExtension.z, LBody.z, 1e-10);
        }

        #endregion

        #region RigidBody Integration

        [TestMethod]
        public void IntegrateSemiImplicitEuler_FreeFall()
        {
            var body = RigidBody.CreateSolidSphere(10, 1);
            body.Position = new Vector(0, 0, 100);

            double g = PhysicsConstants.GravitationalAcceleration;
            double dt = 0.001;
            double totalTime = 2.0;
            int steps = (int)(totalTime / dt);

            for (int i = 0; i < steps; i++)
            {
                body.ApplyForce(new Vector(0, 0, -g * body.Mass));
                body.IntegrateSemiImplicitEuler(dt);
            }

            double expectedZ = 100 - 0.5 * g * totalTime * totalTime;
            Assert.AreEqual(expectedZ, body.Position.z, 0.1);
            Assert.AreEqual(-g * totalTime, body.Velocity.z, 0.1);
        }

        [TestMethod]
        public void IntegrateSemiImplicitEuler_StaticBody_DoesNotMove()
        {
            var body = RigidBody.CreateStatic(new Vector(5, 0, 0));
            body.ApplyForce(new Vector(1000, 0, 0));
            body.IntegrateSemiImplicitEuler(1.0);

            Assert.AreEqual(5, body.Position.x, 1e-10);
            Assert.AreEqual(0, body.Velocity.x, 1e-10);
        }

        [TestMethod]
        public void IntegrateSemiImplicitEuler_ConstantVelocity()
        {
            var body = RigidBody.CreateSolidSphere(5, 1);
            body.Velocity = new Vector(10, 0, 0);

            body.IntegrateSemiImplicitEuler(2.0);

            Assert.AreEqual(20, body.Position.x, 1e-10);
            Assert.AreEqual(10, body.Velocity.x, 1e-10);
        }

        [TestMethod]
        public void IntegrateSemiImplicitEuler_AngularVelocity()
        {
            var body = RigidBody.CreateSolidSphere(10, 1);
            // I = 2/5*10*1 = 4, torque = (0,0,8) → α = 2 rad/s²
            body.ApplyTorque(new Vector(0, 0, 8));
            body.IntegrateSemiImplicitEuler(1.0);

            Assert.AreEqual(2, body.AngularVelocity.z, 1e-10);
        }

        [TestMethod]
        public void IntegrateVelocityVerlet_FreeFall_ExactForConstantForce()
        {
            var body = RigidBody.CreateSolidSphere(10, 1);
            body.Position = new Vector(0, 0, 100);
            double g = PhysicsConstants.GravitationalAcceleration;

            Func<RigidBody, (Vector, Vector)> gravity = b =>
                (new Vector(0, 0, -g * b.Mass), new Vector(0, 0, 0));

            double dt = 0.01;
            int steps = 200; // 2 seconds
            for (int i = 0; i < steps; i++)
                body.IntegrateVelocityVerlet(gravity, dt);

            double t = steps * dt;
            double expectedZ = 100 - 0.5 * g * t * t;
            Assert.AreEqual(expectedZ, body.Position.z, 1e-10); // exact for constant force
        }

        [TestMethod]
        public void IntegrateVelocityVerlet_Spring_ConservesEnergy()
        {
            // 1D spring: F = -kx, k=1, m=1 → ω=1
            var body = RigidBody.CreateSolidSphere(1, 1);
            body.Position = new Vector(1, 0, 0); // stretched spring

            double k = 1.0;
            Func<RigidBody, (Vector, Vector)> spring = b =>
                (new Vector(-k * b.Position.x, 0, 0), new Vector(0, 0, 0));

            double E0 = 0.5 * k * 1.0; // initial PE, KE=0

            double dt = 0.01;
            for (int i = 0; i < 10000; i++) // 100 seconds
                body.IntegrateVelocityVerlet(spring, dt);

            double PE = 0.5 * k * body.Position.x * body.Position.x;
            double KE = 0.5 * body.Mass * body.Velocity.Dot(body.Velocity);
            Assert.AreEqual(E0, PE + KE, 1e-4);
        }

        [TestMethod]
        public void IntegrateVelocityVerlet_StaticBody_DoesNotMove()
        {
            var body = RigidBody.CreateStatic(new Vector(3, 0, 0));
            body.ApplyForce(new Vector(999, 0, 0));

            Func<RigidBody, (Vector, Vector)> force = _ =>
                (new Vector(999, 0, 0), new Vector(0, 0, 0));

            body.IntegrateVelocityVerlet(force, 1.0);

            Assert.AreEqual(3, body.Position.x, 1e-10);
        }

        [TestMethod]
        public void IntegrateEuler_FreeFall()
        {
            var body = RigidBody.CreateSolidSphere(10, 1);
            body.Position = new Vector(0, 0, 100);
            double g = PhysicsConstants.GravitationalAcceleration;

            double dt = 0.001;
            int steps = 2000;
            for (int i = 0; i < steps; i++)
            {
                body.ApplyForce(new Vector(0, 0, -g * body.Mass));
                body.IntegrateEuler(dt);
            }

            double t = steps * dt;
            double expectedZ = 100 - 0.5 * g * t * t;
            Assert.AreEqual(expectedZ, body.Position.z, 0.5);
        }

        [TestMethod]
        public void IntegrateSemiImplicitEuler_Spring_MoreStableThanEuler()
        {
            double k = 10.0;
            double dt = 0.05;
            int steps = 2000;

            // Semi-implicit
            var bodySI = RigidBody.CreateSolidSphere(1, 1);
            bodySI.Position = new Vector(1, 0, 0);
            for (int i = 0; i < steps; i++)
            {
                bodySI.ApplyForce(new Vector(-k * bodySI.Position.x, 0, 0));
                bodySI.IntegrateSemiImplicitEuler(dt);
            }

            // Explicit Euler
            var bodyE = RigidBody.CreateSolidSphere(1, 1);
            bodyE.Position = new Vector(1, 0, 0);
            for (int i = 0; i < steps; i++)
            {
                bodyE.ApplyForce(new Vector(-k * bodyE.Position.x, 0, 0));
                bodyE.IntegrateEuler(dt);
            }

            // Semi-implicit should have bounded energy, explicit Euler blows up
            double eSI = 0.5 * k * bodySI.Position.x * bodySI.Position.x
                + 0.5 * bodySI.Velocity.Dot(bodySI.Velocity);
            double eE = 0.5 * k * bodyE.Position.x * bodyE.Position.x
                + 0.5 * bodyE.Velocity.Dot(bodyE.Velocity);

            double E0 = 0.5 * k; // initial energy
            Assert.IsTrue(Math.Abs(eSI - E0) < Math.Abs(eE - E0));
        }

        #endregion

        #region Common Force Models

        [TestMethod]
        public void SpringForce_AtRest_ZeroForce()
        {
            double k = 50;
            var pos = new Vector(3, 0, 0);
            var anchor = new Vector(0, 0, 0);
            double restLength = 3.0; // exactly at rest length

            var F = k.SpringForce(restLength, pos, anchor);

            Assert.AreEqual(0, F.x, 1e-10);
            Assert.AreEqual(0, F.y, 1e-10);
            Assert.AreEqual(0, F.z, 1e-10);
        }

        [TestMethod]
        public void SpringForce_Stretched_PullsBack()
        {
            double k = 100;
            var pos = new Vector(5, 0, 0);
            var anchor = new Vector(0, 0, 0);
            double restLength = 3.0;

            var F = k.SpringForce(restLength, pos, anchor);

            // Stretched by 2m → F = -100*2 = -200 along x
            Assert.AreEqual(-200, F.x, 1e-10);
            Assert.AreEqual(0, F.y, 1e-10);
        }

        [TestMethod]
        public void SpringForce_Compressed_PushesAway()
        {
            double k = 100;
            var pos = new Vector(1, 0, 0);
            var anchor = new Vector(0, 0, 0);
            double restLength = 3.0;

            var F = k.SpringForce(restLength, pos, anchor);

            // Compressed by 2m → F = +200 along x (pushes away from anchor)
            Assert.AreEqual(200, F.x, 1e-10);
        }

        [TestMethod]
        public void SpringForce_3D_DirectionCorrect()
        {
            double k = 10;
            var pos = new Vector(3, 4, 0);  // distance = 5 from origin
            var anchor = new Vector(0, 0, 0);
            double restLength = 0;

            var F = k.SpringForce(restLength, pos, anchor);

            // F = -10 * 5 * (3/5, 4/5, 0) = (-30, -40, 0)
            Assert.AreEqual(-30, F.x, 1e-10);
            Assert.AreEqual(-40, F.y, 1e-10);
        }

        [TestMethod]
        public void SpringForce_WithVerlet_ConservesEnergy()
        {
            // 1D spring-mass: m=1, k=4, x0=2 → ω=2, E0 = ½k·x0² = 8
            var body = RigidBody.CreateSolidSphere(1, 1);
            body.Position = new Vector(2, 0, 0);
            double k = 4.0;

            Func<RigidBody, (Vector, Vector)> spring = b =>
                (k.SpringForce(0, b.Position, new Vector(0, 0, 0)), new Vector(0, 0, 0));

            double dt = 0.001;
            for (int i = 0; i < 50000; i++) // 50 seconds ≈ 16 oscillation periods
                body.IntegrateVelocityVerlet(spring, dt);

            double E = 0.5 * k * body.Position.Dot(body.Position)
                     + 0.5 * body.Velocity.Dot(body.Velocity);
            Assert.AreEqual(8.0, E, 1e-4);
        }

        [TestMethod]
        public void DampingForce_OpposesVelocity()
        {
            double c = 5.0;
            var v = new Vector(3, -4, 1);

            var F = c.DampingForce(v);

            Assert.AreEqual(-15, F.x, 1e-10);
            Assert.AreEqual(20, F.y, 1e-10);
            Assert.AreEqual(-5, F.z, 1e-10);
        }

        [TestMethod]
        public void DampingForce_ZeroVelocity_ZeroForce()
        {
            var F = 10.0.DampingForce(new Vector(0, 0, 0));

            Assert.AreEqual(0, F.x, 1e-10);
            Assert.AreEqual(0, F.y, 1e-10);
            Assert.AreEqual(0, F.z, 1e-10);
        }

        [TestMethod]
        public void DampingForce_DissipatesEnergy()
        {
            // Spring + damper → energy should decrease over time
            var body = RigidBody.CreateSolidSphere(1, 1);
            body.Position = new Vector(2, 0, 0);
            double k = 10.0, c = 0.5;

            double E0 = 0.5 * k * 4; // initial PE = 20

            Func<RigidBody, (Vector, Vector)> forces = b =>
            {
                var spring = k.SpringForce(0, b.Position, new Vector(0, 0, 0));
                var damp = c.DampingForce(b.Velocity);
                return (spring + damp, new Vector(0, 0, 0));
            };

            double dt = 0.001;
            for (int i = 0; i < 10000; i++)
                body.IntegrateVelocityVerlet(forces, dt);

            double E = 0.5 * k * body.Position.Dot(body.Position)
                     + 0.5 * body.Velocity.Dot(body.Velocity);
            Assert.IsTrue(E < E0); // energy must have decreased
        }

        [TestMethod]
        public void DragForce_Magnitude_MatchesFormula()
        {
            double Cd = 0.47; // sphere
            double rho = 1.225; // air
            double A = Math.PI * 0.05 * 0.05; // r=5cm
            var v = new Vector(10, 0, 0);

            var F = Cd.DragForce(rho, A, v);

            double expectedMag = 0.5 * Cd * rho * A * 100; // |v|²=100
            Assert.AreEqual(-expectedMag, F.x, 1e-10);
            Assert.AreEqual(0, F.y, 1e-10);
        }

        [TestMethod]
        public void DragForce_ZeroVelocity_ZeroForce()
        {
            var F = 0.47.DragForce(1.225, 1.0, new Vector(0, 0, 0));

            Assert.AreEqual(0, F.GetMagnitude(), 1e-10);
        }

        [TestMethod]
        public void DragForce_OpposesMotion()
        {
            var v = new Vector(3, 4, 0); // speed = 5
            var F = 1.0.DragForce(1.0, 1.0, v);

            // F should be anti-parallel to v
            double dot = F.Dot(v);
            Assert.IsTrue(dot < 0);
        }

        [TestMethod]
        public void DragForce_QuadraticScaling()
        {
            double Cd = 1.0, rho = 1.0, A = 1.0;

            var F1 = Cd.DragForce(rho, A, new Vector(1, 0, 0));
            var F2 = Cd.DragForce(rho, A, new Vector(2, 0, 0));

            // Drag scales as v²: doubling speed → 4× force
            Assert.AreEqual(4.0 * Math.Abs(F1.x), Math.Abs(F2.x), 1e-10);
        }

        [TestMethod]
        public void FrictionForce_Kinetic_OpposesVelocity()
        {
            double mu = 0.3;
            double N = 98.0; // ~10 kg on Earth
            var v = new Vector(5, 0, 0);

            var F = mu.FrictionForce(N, v);

            Assert.AreEqual(-0.3 * 98, F.x, 1e-10);
            Assert.AreEqual(0, F.y, 1e-10);
        }

        [TestMethod]
        public void FrictionForce_Kinetic_MagnitudeIsMuN()
        {
            double mu = 0.5;
            double N = 50;
            var v = new Vector(3, 4, 0); // speed = 5

            var F = mu.FrictionForce(N, v);

            Assert.AreEqual(mu * N, F.GetMagnitude(), 1e-10);
        }

        [TestMethod]
        public void FrictionForce_Static_CancelsAppliedForce_WhenBelowLimit()
        {
            double mu = 0.5;
            double N = 100;
            var v = new Vector(0, 0, 0); // stationary
            var applied = new Vector(20, 0, 0); // below μN = 50

            var F = mu.FrictionForce(N, v, applied);

            // Static friction fully cancels applied force
            Assert.AreEqual(-20, F.x, 1e-10);
        }

        [TestMethod]
        public void FrictionForce_Static_CappedAtMuN()
        {
            double mu = 0.5;
            double N = 100;
            var v = new Vector(0, 0, 0);
            var applied = new Vector(80, 0, 0); // above μN = 50

            var F = mu.FrictionForce(N, v, applied);

            // Capped at μN = 50
            Assert.AreEqual(-50, F.x, 1e-10);
        }

        [TestMethod]
        public void FrictionForce_ZeroNormal_ZeroFriction()
        {
            // No normal force → no friction (e.g. object in freefall)
            var F = 0.5.FrictionForce(0, new Vector(5, 0, 0));

            Assert.AreEqual(0, F.GetMagnitude(), 1e-10);
        }

        #endregion
    }
}
