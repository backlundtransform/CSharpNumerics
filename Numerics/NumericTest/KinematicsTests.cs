using CSharpNumerics.Physics;
using Numerics.Objects;


namespace NumericTest
{
    [TestClass]
    public class KinematicsTests
    {
        [TestMethod]
        public void FreeFallVelocity_From10Meters_IsCorrect()
        {
            var v = 10.0.FreeFallVelocity();
            Assert.AreEqual(Math.Sqrt(2 * 9.80665 * 10), v, 6);
        }

        [TestMethod]
        public void PositionFromConstantAcceleration_Scalar()
        {
            var s = 2.0.PositionFromConstantAcceleration(
                time: 3,
                initialVelocity: 1,
                initialPosition: 0);

            // s = 0 + 1*3 + 0.5*2*9 = 12
            Assert.AreEqual(12.0, s, 6);
        }

        [TestMethod]
        public void PositionFromConstantAcceleration_Vector()
        {
            var a = new Vector(1, 0, 0);
            var v0 = new Vector(0, 0, 0);
            var s0 = new Vector(0, 0, 0);

            var s = a.PositionFromConstantAcceleration(2, v0, s0);

            // s = 0 + 0 + 0.5 * 1 * 4 = (2,0,0)
            Assert.AreEqual(2, s.x, 6);
        }

        [TestMethod]
        public void VelocityFromConstantAcceleration_Vector()
        {
            var a = new Vector(1, 0, 0);
            var v0 = new Vector(0, 0, 0);
            var v = a.VelocityFromConstantAcceleration(3, v0);
            // v = 0 + 1 * 3 = (3,0,0)
            Assert.AreEqual(3, v.x, 6);
        }
        [TestMethod]
        public void VelocityFromConstantAcceleration_Scalar()
        {
            var v = 2.0.VelocityFromConstantAcceleration(
                time: 4,
                initialVelocity: 1);
            // v = 1 + 2*4 = 9
            Assert.AreEqual(9.0, v, 6);
        }
        [TestMethod]
        public void PositionFromConstantVelocity_Scalar()
        {
            var s = 3.0.PositionFromConstantVelocity(
                time: 5,
                initialPosition: 2);
            // s = 2 + 3*5 = 17
            Assert.AreEqual(17.0, s, 6);
        }
        [TestMethod]
        public void PositionFromConstantVelocity_Vector()
        {
            var v = new Vector(2, 0, 0);
            var s0 = new Vector(1, 0, 0);
            var s = v.PositionFromConstantVelocity(4, s0);
            // s = (1,0,0) + (2,0,0)*4 = (9,0,0)
            Assert.AreEqual(9, s.x, 6);
        }

        [TestMethod]
        public void FreeFallVelocity_Vector_DefaultDirection()
        {
            var v = 10.0.FreeFallVelocity(); 
            var vVec = 10.0.FreeFallVelocity(new Vector(0, 0, -1));

        
            Assert.AreEqual(v, vVec.GetMagnitude(), 6);

          
            Assert.AreEqual(0, vVec.x, 6);
            Assert.AreEqual(0, vVec.y, 6);
            Assert.AreEqual(-1, vVec.z / vVec.GetMagnitude(), 6); // unit vector
        }

        [TestMethod]
        public void CentripetalAcceleration_Vector()
        {
            var velocity = new Vector(2, 0, 0);
            var radius = new Vector(0, 3, 0); 

            var a = velocity.CentripetalAcceleration(radius);

            var expected = -(2 * 2 / 3.0); 
            Assert.AreEqual(expected, a.y, 6);
            Assert.AreEqual(0, a.x, 6);
            Assert.AreEqual(0, a.z, 6);
        }

        #region Projectile Motion

        [TestMethod]
        public void ProjectileVelocityFromAngle_45Degrees()
        {
            double speed = 20.0;
            double angle = Math.PI / 4; // 45°
            var v0 = speed.ProjectileVelocityFromAngle(angle);

            double expected = 20 * Math.Cos(Math.PI / 4);
            Assert.AreEqual(expected, v0.x, 1e-10);
            Assert.AreEqual(0, v0.y, 1e-10);
            Assert.AreEqual(expected, v0.z, 1e-10);
        }

        [TestMethod]
        public void ProjectilePosition_AtLaunch_IsOrigin()
        {
            var v0 = new Vector(10, 0, 10);
            var pos = v0.ProjectilePosition(0);

            Assert.AreEqual(0, pos.x, 1e-10);
            Assert.AreEqual(0, pos.y, 1e-10);
            Assert.AreEqual(0, pos.z, 1e-10);
        }

        [TestMethod]
        public void ProjectilePosition_WithInitialHeight()
        {
            var v0 = new Vector(10, 0, 0); // horizontal launch
            var pos = v0.ProjectilePosition(0, initialHeight: 50);

            Assert.AreEqual(0, pos.x, 1e-10);
            Assert.AreEqual(50, pos.z, 1e-10);
        }

        [TestMethod]
        public void ProjectileVelocity_AtLaunch_EqualsInitial()
        {
            var v0 = new Vector(10, 5, 15);
            var v = v0.ProjectileVelocity(0);

            Assert.AreEqual(v0.x, v.x, 1e-10);
            Assert.AreEqual(v0.y, v.y, 1e-10);
            Assert.AreEqual(v0.z, v.z, 1e-10);
        }

        [TestMethod]
        public void ProjectileVelocity_HorizontalUnchanged()
        {
            var v0 = new Vector(10, 5, 0);
            var v = v0.ProjectileVelocity(2);

            Assert.AreEqual(10, v.x, 1e-10);
            Assert.AreEqual(5, v.y, 1e-10);
        }

        [TestMethod]
        public void ProjectileTimeOfFlight_Scalar_45Degrees()
        {
            // T = 2 * v0 * sin(45°) / g
            double speed = 20.0;
            double angle = Math.PI / 4;
            double g = 9.80665;

            double T = speed.ProjectileTimeOfFlight(angle);
            double expected = 2 * speed * Math.Sin(angle) / g;
            Assert.AreEqual(expected, T, 1e-10);
        }

        [TestMethod]
        public void ProjectileMaxHeight_Scalar_45Degrees()
        {
            // H = v0² sin²(θ) / (2g)
            double speed = 20.0;
            double angle = Math.PI / 4;
            double g = 9.80665;

            double H = speed.ProjectileMaxHeight(angle);
            double vz = speed * Math.Sin(angle);
            double expected = (vz * vz) / (2 * g);
            Assert.AreEqual(expected, H, 1e-10);
        }

        [TestMethod]
        public void ProjectileRange_Scalar_45Degrees_IsMaximum()
        {
            // Range is maximized at 45°
            double speed = 20.0;
            double range45 = speed.ProjectileRange(Math.PI / 4);
            double range30 = speed.ProjectileRange(Math.PI / 6);
            double range60 = speed.ProjectileRange(Math.PI / 3);

            Assert.IsTrue(range45 > range30);
            Assert.IsTrue(range45 > range60);
        }

        [TestMethod]
        public void ProjectileRange_Scalar_EqualsFormula()
        {
            // R = v₀² sin(2θ) / g
            double speed = 30.0;
            double angle = Math.PI / 6;
            double g = 9.80665;

            double R = speed.ProjectileRange(angle);
            double expected = (speed * speed * Math.Sin(2 * angle)) / g;
            Assert.AreEqual(expected, R, 1e-10);
        }

        [TestMethod]
        public void ProjectileTimeOfFlight_Vector_MatchesScalar()
        {
            double speed = 25.0;
            double angle = Math.PI / 3;

            var v0 = speed.ProjectileVelocityFromAngle(angle);
            double Tvec = v0.ProjectileTimeOfFlight();
            double Tscalar = speed.ProjectileTimeOfFlight(angle);

            Assert.AreEqual(Tscalar, Tvec, 1e-10);
        }

        [TestMethod]
        public void ProjectileMaxHeight_Vector_MatchesScalar()
        {
            double speed = 25.0;
            double angle = Math.PI / 3;

            var v0 = speed.ProjectileVelocityFromAngle(angle);
            double Hvec = v0.ProjectileMaxHeight();
            double Hscalar = speed.ProjectileMaxHeight(angle);

            Assert.AreEqual(Hscalar, Hvec, 1e-10);
        }

        [TestMethod]
        public void ProjectileRange_Vector_MatchesScalar()
        {
            double speed = 25.0;
            double angle = Math.PI / 3;

            var v0 = speed.ProjectileVelocityFromAngle(angle);
            double Rvec = v0.ProjectileRange();
            double Rscalar = speed.ProjectileRange(angle);

            Assert.AreEqual(Rscalar, Rvec, 1e-6);
        }

        [TestMethod]
        public void ProjectilePosition_AtTimeOfFlight_ReturnsToGround()
        {
            var v0 = 20.0.ProjectileVelocityFromAngle(Math.PI / 4);
            double T = v0.ProjectileTimeOfFlight();
            var pos = v0.ProjectilePosition(T);

            Assert.AreEqual(0, pos.z, 1e-6);
        }

        [TestMethod]
        public void ProjectilePosition_AtPeak_MatchesMaxHeight()
        {
            var v0 = 20.0.ProjectileVelocityFromAngle(Math.PI / 4);
            double g = 9.80665;
            double tPeak = v0.z / g;

            var pos = v0.ProjectilePosition(tPeak);
            double maxH = v0.ProjectileMaxHeight();

            Assert.AreEqual(maxH, pos.z, 1e-6);
        }

        [TestMethod]
        public void ProjectileVelocity_AtPeak_VerticalIsZero()
        {
            var v0 = 20.0.ProjectileVelocityFromAngle(Math.PI / 4);
            double g = 9.80665;
            double tPeak = v0.z / g;

            var v = v0.ProjectileVelocity(tPeak);

            Assert.AreEqual(0, v.z, 1e-6);
            Assert.AreEqual(v0.x, v.x, 1e-10); // horizontal unchanged
        }

        [TestMethod]
        public void ProjectileRange_WithInitialHeight_GreaterThanFlat()
        {
            var v0 = 15.0.ProjectileVelocityFromAngle(Math.PI / 6);
            double rangeFlat = v0.ProjectileRange(0);
            double rangeElevated = v0.ProjectileRange(10);

            Assert.IsTrue(rangeElevated > rangeFlat);
        }

        #endregion
    }
}
