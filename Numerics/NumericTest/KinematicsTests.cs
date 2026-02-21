using CSharpNumerics.Physics;
using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Numerics.Objects;


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

        [TestMethod]
        public void AngularSpeed_FromTangentialSpeed()
        {
            // v = 10 m/s, r = 5 m → ω = 2 rad/s
            double omega = 10.0.AngularSpeed(5);
            Assert.AreEqual(2, omega, 1e-10);
        }

        [TestMethod]
        public void AngularVelocity_InXYPlane()
        {
            // Object at (5,0,0) moving along +Y at 10 m/s → ω along +Z
            var vel = new Vector(0, 10, 0);
            var radius = new Vector(5, 0, 0);

            var omega = vel.AngularVelocity(radius);

            Assert.AreEqual(0, omega.x, 1e-10);
            Assert.AreEqual(0, omega.y, 1e-10);
            Assert.AreEqual(2, omega.z, 1e-10); // ω = v/r = 10/5 = 2
        }

        [TestMethod]
        public void AngularVelocity_MagnitudeEqualsScalar()
        {
            var vel = new Vector(0, 15, 0);
            var radius = new Vector(3, 0, 0);

            double scalarOmega = vel.GetMagnitude().AngularSpeed(radius.GetMagnitude());
            double vectorOmega = vel.AngularVelocity(radius).GetMagnitude();

            Assert.AreEqual(scalarOmega, vectorOmega, 1e-10);
        }

        [TestMethod]
        public void Period_FromSpeedAndRadius()
        {
            // v = 2π m/s, r = 1 m → T = 2π·1 / 2π = 1 s
            double T = (2 * Math.PI).Period(1);
            Assert.AreEqual(1, T, 1e-10);
        }

        [TestMethod]
        public void Period_GeneralCase()
        {
            // T = 2πr/v
            double v = 10.0;
            double r = 5.0;
            double T = v.Period(r);
            Assert.AreEqual(2 * Math.PI * r / v, T, 1e-10);
        }

        [TestMethod]
        public void Frequency_IsInverseOfPeriod()
        {
            double v = 12.0;
            double r = 4.0;
            double T = v.Period(r);
            double f = v.Frequency(r);
            Assert.AreEqual(1.0 / T, f, 1e-10);
        }

        [TestMethod]
        public void TangentialVelocity_FromAngularVelocity()
        {
            // ω along +Z, object at (3,0,0) → v = ω × r = (0,0,2) × (3,0,0) = (0,6,0)
            var omega = new Vector(0, 0, 2);
            var radius = new Vector(3, 0, 0);

            var v = omega.TangentialVelocity(radius);

            Assert.AreEqual(0, v.x, 1e-10);
            Assert.AreEqual(6, v.y, 1e-10);
            Assert.AreEqual(0, v.z, 1e-10);
        }

        [TestMethod]
        public void TangentialVelocity_MagnitudeEqualsOmegaTimesR()
        {
            var omega = new Vector(0, 0, 5);
            var radius = new Vector(4, 0, 0);

            var v = omega.TangentialVelocity(radius);

            Assert.AreEqual(20, v.GetMagnitude(), 1e-10); // |v| = ω·r = 5·4
        }

        [TestMethod]
        public void TangentialVelocity_PerpendicularToRadius()
        {
            var omega = new Vector(0, 0, 3);
            var radius = new Vector(2, 1, 0);

            var v = omega.TangentialVelocity(radius);

            double dot = v.Dot(radius);
            Assert.AreEqual(0, dot, 1e-10);
        }

        [TestMethod]
        public void AngularVelocity_TangentialVelocity_RoundTrip()
        {
            // Start with ω and r, get v, then recover ω from v and r
            var omegaIn = new Vector(0, 0, 4);
            var radius = new Vector(3, 0, 0);

            var v = omegaIn.TangentialVelocity(radius);
            var omegaOut = v.AngularVelocity(radius);

            Assert.AreEqual(omegaIn.x, omegaOut.x, 1e-10);
            Assert.AreEqual(omegaIn.y, omegaOut.y, 1e-10);
            Assert.AreEqual(omegaIn.z, omegaOut.z, 1e-10);
        }

        [TestMethod]
        public void AngularSpeed_Period_Consistent()
        {
            // ω = 2π/T ↔ T = 2π/ω
            double v = 8.0;
            double r = 2.0;
            double omega = v.AngularSpeed(r);
            double T = v.Period(r);
            Assert.AreEqual(2 * Math.PI / omega, T, 1e-10);
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

        #region Orbital Mechanics

        [TestMethod]
        public void GravitationalFieldStrength_EarthSurface()
        {
            // g = GM/R² ≈ 9.82 m/s² (slightly differs from standard g due to rounding)
            double g = PhysicsConstants.EarthMass.GravitationalFieldStrength(PhysicsConstants.EarthRadius);
            Assert.AreEqual(9.82, g, 0.01);
        }

        [TestMethod]
        public void GravitationalForce_EarthMoon()
        {
            // F = G·M·m/r²
            double r = 3.844e8; // Earth–Moon distance in meters
            double F = PhysicsConstants.EarthMass.GravitationalForce(PhysicsConstants.MoonMass, r);
            // ≈ 1.98e20 N
            Assert.AreEqual(1.98e20, F, 0.02e20);
        }

        [TestMethod]
        public void OrbitalSpeed_ISS()
        {
            // ISS orbits at ~408 km altitude → r = R_earth + 408km
            double r = PhysicsConstants.EarthRadius + 408000;
            double v = PhysicsConstants.EarthMass.OrbitalSpeed(r);
            // Expected ~7660 m/s
            Assert.AreEqual(7660, v, 20);
        }

        [TestMethod]
        public void OrbitalPeriod_ISS()
        {
            // ISS orbital period ≈ 92.6 min ≈ 5556 s
            double r = PhysicsConstants.EarthRadius + 408000;
            double T = PhysicsConstants.EarthMass.OrbitalPeriod(r);
            Assert.AreEqual(5556, T, 20);
        }

        [TestMethod]
        public void EscapeVelocity_Earth()
        {
            double v = PhysicsConstants.EarthMass.EscapeVelocity(PhysicsConstants.EarthRadius);
            // Expected ≈ 11186 m/s (matches PhysicsConstants.EscapeVelocityEarth)
            Assert.AreEqual(PhysicsConstants.EscapeVelocityEarth, v, 10);
        }

        [TestMethod]
        public void EscapeVelocity_IsSqrt2TimesOrbitalSpeed()
        {
            double r = PhysicsConstants.EarthRadius;
            double vOrb = PhysicsConstants.EarthMass.OrbitalSpeed(r);
            double vEsc = PhysicsConstants.EarthMass.EscapeVelocity(r);
            Assert.AreEqual(Math.Sqrt(2), vEsc / vOrb, 1e-10);
        }

        [TestMethod]
        public void OrbitalPosition_AtTimeZero_IsOnXAxis()
        {
            double r = 1e7;
            var pos = PhysicsConstants.EarthMass.OrbitalPosition(r, 0);

            Assert.AreEqual(r, pos.x, 1e-3);
            Assert.AreEqual(0, pos.y, 1e-3);
            Assert.AreEqual(0, pos.z, 1e-10);
        }

        [TestMethod]
        public void OrbitalPosition_MagnitudeIsConstant()
        {
            double r = 1e7;
            double M = PhysicsConstants.EarthMass;

            var pos1 = M.OrbitalPosition(r, 0);
            var pos2 = M.OrbitalPosition(r, 1000);
            var pos3 = M.OrbitalPosition(r, 2500);

            Assert.AreEqual(r, pos1.GetMagnitude(), 1e-3);
            Assert.AreEqual(r, pos2.GetMagnitude(), 1e-3);
            Assert.AreEqual(r, pos3.GetMagnitude(), 1e-3);
        }

        [TestMethod]
        public void OrbitalVelocity_MagnitudeMatchesOrbitalSpeed()
        {
            double r = 1e7;
            double M = PhysicsConstants.EarthMass;
            double expectedSpeed = M.OrbitalSpeed(r);

            var v = M.OrbitalVelocity(r, 1234);

            Assert.AreEqual(expectedSpeed, v.GetMagnitude(), 1e-3);
        }

        [TestMethod]
        public void OrbitalVelocity_PerpendicularToPosition()
        {
            double r = 1e7;
            double M = PhysicsConstants.EarthMass;

            var pos = M.OrbitalPosition(r, 500);
            var vel = M.OrbitalVelocity(r, 500);

            // Dot product should be zero (perpendicular)
            double dot = pos.Dot(vel);
            Assert.AreEqual(0, dot, 1e-3);
        }

        [TestMethod]
        public void OrbitalAcceleration_PointsTowardCenter()
        {
            double r = 1e7;
            double M = PhysicsConstants.EarthMass;

            var pos = M.OrbitalPosition(r, 800);
            var acc = M.OrbitalAcceleration(r, 800);

            // Acceleration should be antiparallel to position (points inward)
            // Normalize both and check: â = -r̂
            var rHat = pos.GetUnitVector();
            var aHat = acc.GetUnitVector();

            Assert.AreEqual(-rHat.x, aHat.x, 1e-10);
            Assert.AreEqual(-rHat.y, aHat.y, 1e-10);
        }

        [TestMethod]
        public void OrbitalAcceleration_MagnitudeEqualsGMOverR2()
        {
            double r = 1e7;
            double M = PhysicsConstants.EarthMass;

            var acc = M.OrbitalAcceleration(r, 600);
            double expectedMag = PhysicsConstants.GravitationalConstant * M / (r * r);

            Assert.AreEqual(expectedMag, acc.GetMagnitude(), 1e-6);
        }

        [TestMethod]
        public void OrbitalPosition_AfterFullPeriod_ReturnsToStart()
        {
            double r = 1e7;
            double M = PhysicsConstants.EarthMass;
            double T = M.OrbitalPeriod(r);

            var posStart = M.OrbitalPosition(r, 0);
            var posEnd = M.OrbitalPosition(r, T);

            Assert.AreEqual(posStart.x, posEnd.x, 1e-3);
            Assert.AreEqual(posStart.y, posEnd.y, 1e-3);
        }

        #endregion

        #region Relative Motion

        [TestMethod]
        public void RelativeVelocity_StationaryReference()
        {
            var vObj = new Vector(10, 5, 0);
            var vRef = new Vector(0, 0, 0);

            var vRel = vObj.RelativeVelocity(vRef);

            Assert.AreEqual(10, vRel.x, 1e-10);
            Assert.AreEqual(5, vRel.y, 1e-10);
            Assert.AreEqual(0, vRel.z, 1e-10);
        }

        [TestMethod]
        public void RelativeVelocity_SameVelocity_IsZero()
        {
            var v = new Vector(30, 20, 10);

            var vRel = v.RelativeVelocity(v);

            Assert.AreEqual(0, vRel.x, 1e-10);
            Assert.AreEqual(0, vRel.y, 1e-10);
            Assert.AreEqual(0, vRel.z, 1e-10);
        }

        [TestMethod]
        public void RelativeVelocity_OppositeDirections()
        {
            // Two cars heading toward each other at 30 m/s
            var vA = new Vector(30, 0, 0);
            var vB = new Vector(-30, 0, 0);

            var vRel = vA.RelativeVelocity(vB);

            Assert.AreEqual(60, vRel.x, 1e-10);
        }

        [TestMethod]
        public void RelativePosition_Separation()
        {
            var pObj = new Vector(100, 50, 0);
            var pRef = new Vector(20, 10, 0);

            var rRel = pObj.RelativePosition(pRef);

            Assert.AreEqual(80, rRel.x, 1e-10);
            Assert.AreEqual(40, rRel.y, 1e-10);
        }

        [TestMethod]
        public void RelativeAcceleration_Subtraction()
        {
            var aObj = new Vector(5, 0, -9.8);
            var aRef = new Vector(0, 0, -9.8);

            var aRel = aObj.RelativeAcceleration(aRef);

            Assert.AreEqual(5, aRel.x, 1e-10);
            Assert.AreEqual(0, aRel.y, 1e-10);
            Assert.AreEqual(0, aRel.z, 1e-10); // gravity cancels
        }

        [TestMethod]
        public void ClosingSpeed_HeadOn_IsPositive()
        {
            var vA = new Vector(30, 0, 0);
            var vB = new Vector(-20, 0, 0);
            var pA = new Vector(-100, 0, 0);
            var pB = new Vector(100, 0, 0);

            double cs = vA.ClosingSpeed(vB, pA, pB);

            // Objects approaching: closing speed should be positive
            Assert.AreEqual(50, cs, 1e-10);
        }

        [TestMethod]
        public void ClosingSpeed_MovingApart_IsNegative()
        {
            var vA = new Vector(-30, 0, 0);
            var vB = new Vector(20, 0, 0);
            var pA = new Vector(-100, 0, 0);
            var pB = new Vector(100, 0, 0);

            double cs = vA.ClosingSpeed(vB, pA, pB);

            Assert.IsTrue(cs < 0);
        }

        [TestMethod]
        public void RelativePositionAtTime_ConstantVelocity()
        {
            var vObj = new Vector(10, 0, 0);
            var vRef = new Vector(3, 0, 0);
            var r0Obj = new Vector(0, 0, 0);
            var r0Ref = new Vector(50, 0, 0);

            // At t=0: relative position = -50
            var rel0 = vObj.RelativePositionAtTime(vRef, 0, r0Obj, r0Ref);
            Assert.AreEqual(-50, rel0.x, 1e-10);

            // Relative velocity = 7 m/s, so gap closes in 50/7 s
            double tCatch = 50.0 / 7.0;
            var relCatch = vObj.RelativePositionAtTime(vRef, tCatch, r0Obj, r0Ref);
            Assert.AreEqual(0, relCatch.x, 1e-9);
        }

        [TestMethod]
        public void TimeOfClosestApproach_HeadOn_IsCorrect()
        {
            // Object A at x=-100 moving right at 10 m/s
            // Object B at x=+100 moving left at 10 m/s
            // They meet at t = 100/10 = 10 s (relative speed 20, gap 200 → t=10)
            var vA = new Vector(10, 0, 0);
            var vB = new Vector(-10, 0, 0);
            var pA = new Vector(-100, 0, 0);
            var pB = new Vector(100, 0, 0);

            double t = vA.TimeOfClosestApproach(vB, pA, pB);
            Assert.AreEqual(10, t, 1e-10);
        }

        [TestMethod]
        public void TimeOfClosestApproach_Perpendicular()
        {
            // Object at origin moving along +x
            // Reference at (0, 10, 0) moving along +x at same speed
            // Always parallel, closest approach is t=0
            var vObj = new Vector(5, 0, 0);
            var vRef = new Vector(5, 0, 0);
            var pObj = new Vector(0, 0, 0);
            var pRef = new Vector(0, 10, 0);

            double t = vObj.TimeOfClosestApproach(vRef, pObj, pRef);
            // No relative motion → t=0
            Assert.AreEqual(0, t, 1e-10);
        }

        [TestMethod]
        public void MinimumDistance_HeadOn_IsZero()
        {
            var vA = new Vector(10, 0, 0);
            var vB = new Vector(-10, 0, 0);
            var pA = new Vector(-50, 0, 0);
            var pB = new Vector(50, 0, 0);

            double d = vA.MinimumDistance(vB, pA, pB);
            Assert.AreEqual(0, d, 1e-9);
        }

        [TestMethod]
        public void MinimumDistance_MissingSideways()
        {
            // A moves along +x, B moves along +x with a lateral offset
            var vA = new Vector(10, 0, 0);
            var vB = new Vector(0, 0, 0);
            var pA = new Vector(-100, 0, 0);
            var pB = new Vector(0, 5, 0);

            double d = vA.MinimumDistance(vB, pA, pB);
            // Closest approach at t=10 (when A.x = 0), distance = 5 m (lateral offset)
            Assert.AreEqual(5, d, 1e-9);
        }

        [TestMethod]
        public void MinimumDistance_IsLessOrEqualToInitialDistance()
        {
            var vA = new Vector(3, 1, 0);
            var vB = new Vector(-2, 3, 0);
            var pA = new Vector(10, 20, 0);
            var pB = new Vector(50, 30, 0);

            double initialDist = (pA - pB).GetMagnitude();
            double minDist = vA.MinimumDistance(vB, pA, pB);

            Assert.IsTrue(minDist <= initialDist + 1e-10);
        }

        #endregion
    }
}
