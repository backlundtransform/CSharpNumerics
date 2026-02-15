using CSharpNumerics.Physics;
using CSharpNumerics.Physics.Constants;
using Numerics.Objects;

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
    }
}
