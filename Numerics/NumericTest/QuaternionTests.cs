using Numerics.Objects;

namespace NumericTest
{
    [TestClass]
    public class QuaternionTests
    {
        private const double Tol = 1e-10;

        #region Construction & Identity

        [TestMethod]
        public void Identity_IsUnitQuaternion()
        {
            var q = Quaternion.Identity;
            Assert.AreEqual(1, q.w, Tol);
            Assert.AreEqual(0, q.x, Tol);
            Assert.AreEqual(0, q.y, Tol);
            Assert.AreEqual(0, q.z, Tol);
            Assert.IsTrue(q.IsUnit);
        }

        [TestMethod]
        public void Norm_MatchesExpected()
        {
            var q = new Quaternion(1, 2, 3, 4);
            Assert.AreEqual(Math.Sqrt(30), q.Norm, Tol);
        }

        [TestMethod]
        public void Normalize_ProducesUnitQuaternion()
        {
            var q = new Quaternion(1, 2, 3, 4).Normalize();
            Assert.AreEqual(1.0, q.Norm, Tol);
            Assert.IsTrue(q.IsUnit);
        }

        #endregion

        #region Arithmetic

        [TestMethod]
        public void Addition_ComponentWise()
        {
            var a = new Quaternion(1, 2, 3, 4);
            var b = new Quaternion(5, 6, 7, 8);
            var r = a + b;

            Assert.AreEqual(6, r.w, Tol);
            Assert.AreEqual(8, r.x, Tol);
            Assert.AreEqual(10, r.y, Tol);
            Assert.AreEqual(12, r.z, Tol);
        }

        [TestMethod]
        public void Subtraction_ComponentWise()
        {
            var a = new Quaternion(5, 6, 7, 8);
            var b = new Quaternion(1, 2, 3, 4);
            var r = a - b;

            Assert.AreEqual(4, r.w, Tol);
            Assert.AreEqual(4, r.x, Tol);
            Assert.AreEqual(4, r.y, Tol);
            Assert.AreEqual(4, r.z, Tol);
        }

        [TestMethod]
        public void ScalarMultiplication()
        {
            var q = new Quaternion(1, 2, 3, 4);
            var r = 2.0 * q;

            Assert.AreEqual(2, r.w, Tol);
            Assert.AreEqual(4, r.x, Tol);
            Assert.AreEqual(6, r.y, Tol);
            Assert.AreEqual(8, r.z, Tol);
        }

        [TestMethod]
        public void HamiltonProduct_ij_Equals_k()
        {
            // i·j = k
            var i = new Quaternion(0, 1, 0, 0);
            var j = new Quaternion(0, 0, 1, 0);
            var k = i * j;

            Assert.AreEqual(0, k.w, Tol);
            Assert.AreEqual(0, k.x, Tol);
            Assert.AreEqual(0, k.y, Tol);
            Assert.AreEqual(1, k.z, Tol);
        }

        [TestMethod]
        public void HamiltonProduct_ji_Equals_NegativeK()
        {
            // j·i = -k (non-commutative)
            var i = new Quaternion(0, 1, 0, 0);
            var j = new Quaternion(0, 0, 1, 0);
            var r = j * i;

            Assert.AreEqual(0, r.w, Tol);
            Assert.AreEqual(0, r.x, Tol);
            Assert.AreEqual(0, r.y, Tol);
            Assert.AreEqual(-1, r.z, Tol);
        }

        [TestMethod]
        public void HamiltonProduct_Associative()
        {
            var a = new Quaternion(1, 2, 3, 4).Normalize();
            var b = new Quaternion(5, 6, 7, 8).Normalize();
            var c = new Quaternion(9, 10, 11, 12).Normalize();

            var ab_c = (a * b) * c;
            var a_bc = a * (b * c);

            Assert.AreEqual(ab_c.w, a_bc.w, Tol);
            Assert.AreEqual(ab_c.x, a_bc.x, Tol);
            Assert.AreEqual(ab_c.y, a_bc.y, Tol);
            Assert.AreEqual(ab_c.z, a_bc.z, Tol);
        }

        #endregion

        #region Conjugate & Inverse

        [TestMethod]
        public void Conjugate_FlipsImaginaryParts()
        {
            var q = new Quaternion(1, 2, 3, 4);
            var c = q.Conjugate();

            Assert.AreEqual(1, c.w, Tol);
            Assert.AreEqual(-2, c.x, Tol);
            Assert.AreEqual(-3, c.y, Tol);
            Assert.AreEqual(-4, c.z, Tol);
        }

        [TestMethod]
        public void Inverse_TimesOriginal_IsIdentity()
        {
            var q = new Quaternion(1, 2, 3, 4);
            var r = q * q.Inverse();

            Assert.AreEqual(1, r.w, Tol);
            Assert.AreEqual(0, r.x, Tol);
            Assert.AreEqual(0, r.y, Tol);
            Assert.AreEqual(0, r.z, Tol);
        }

        [TestMethod]
        public void UnitQuaternion_Conjugate_EqualsInverse()
        {
            var q = new Quaternion(1, 2, 3, 4).Normalize();
            var conj = q.Conjugate();
            var inv = q.Inverse();

            Assert.AreEqual(conj.w, inv.w, Tol);
            Assert.AreEqual(conj.x, inv.x, Tol);
            Assert.AreEqual(conj.y, inv.y, Tol);
            Assert.AreEqual(conj.z, inv.z, Tol);
        }

        #endregion

        #region Rotation

        [TestMethod]
        public void Rotate_90DegreesAboutZ_RotatesXToY()
        {
            var q = Quaternion.FromAxisAngle(new Vector(0, 0, 1), Math.PI / 2);
            var v = new Vector(1, 0, 0);
            var rotated = q.Rotate(v);

            Assert.AreEqual(0, rotated.x, 1e-10);
            Assert.AreEqual(1, rotated.y, 1e-10);
            Assert.AreEqual(0, rotated.z, 1e-10);
        }

        [TestMethod]
        public void Rotate_180DegreesAboutY_FlipsX()
        {
            var q = Quaternion.FromAxisAngle(new Vector(0, 1, 0), Math.PI);
            var v = new Vector(1, 0, 0);
            var rotated = q.Rotate(v);

            Assert.AreEqual(-1, rotated.x, 1e-10);
            Assert.AreEqual(0, rotated.y, 1e-10);
            Assert.AreEqual(0, rotated.z, 1e-10);
        }

        [TestMethod]
        public void Rotate_Identity_NoChange()
        {
            var v = new Vector(3, 4, 5);
            var rotated = Quaternion.Identity.Rotate(v);

            Assert.AreEqual(3, rotated.x, Tol);
            Assert.AreEqual(4, rotated.y, Tol);
            Assert.AreEqual(5, rotated.z, Tol);
        }

        [TestMethod]
        public void Rotate_PreservesMagnitude()
        {
            var q = Quaternion.FromAxisAngle(new Vector(1, 1, 1), 1.23);
            var v = new Vector(3, 4, 5);
            var rotated = q.Rotate(v);

            Assert.AreEqual(v.GetMagnitude(), rotated.GetMagnitude(), 1e-10);
        }

        [TestMethod]
        public void ComposedRotations_MatchMultiplication()
        {
            // Two successive rotations via quaternion multiplication
            var q1 = Quaternion.FromAxisAngle(new Vector(0, 0, 1), Math.PI / 4);
            var q2 = Quaternion.FromAxisAngle(new Vector(0, 1, 0), Math.PI / 3);
            var combined = q2 * q1; // q2 applied after q1

            var v = new Vector(1, 0, 0);
            var step = q1.Rotate(v);
            var sequential = q2.Rotate(step);
            var direct = combined.Rotate(v);

            Assert.AreEqual(sequential.x, direct.x, 1e-10);
            Assert.AreEqual(sequential.y, direct.y, 1e-10);
            Assert.AreEqual(sequential.z, direct.z, 1e-10);
        }

        #endregion

        #region Matrix Conversion

        [TestMethod]
        public void ToMatrix_Identity_IsIdentityMatrix()
        {
            var m = Quaternion.Identity.ToMatrix();

            Assert.AreEqual(1, m.values[0, 0], Tol);
            Assert.AreEqual(0, m.values[0, 1], Tol);
            Assert.AreEqual(0, m.values[0, 2], Tol);
            Assert.AreEqual(0, m.values[1, 0], Tol);
            Assert.AreEqual(1, m.values[1, 1], Tol);
            Assert.AreEqual(0, m.values[1, 2], Tol);
            Assert.AreEqual(0, m.values[2, 0], Tol);
            Assert.AreEqual(0, m.values[2, 1], Tol);
            Assert.AreEqual(1, m.values[2, 2], Tol);
        }

        [TestMethod]
        public void ToMatrix_MatchesRotateVector()
        {
            var q = Quaternion.FromAxisAngle(new Vector(1, 1, 0), 0.7);
            var v = new Vector(3, -2, 5);

            var rotatedQ = q.Rotate(v);
            var m = q.ToMatrix();
            var rotatedM = m * v;

            Assert.AreEqual(rotatedQ.x, rotatedM.x, 1e-10);
            Assert.AreEqual(rotatedQ.y, rotatedM.y, 1e-10);
            Assert.AreEqual(rotatedQ.z, rotatedM.z, 1e-10);
        }

        [TestMethod]
        public void FromMatrix_RoundTrip()
        {
            var original = Quaternion.FromAxisAngle(new Vector(1, 2, 3), 1.0);
            var m = original.ToMatrix();
            var recovered = Quaternion.FromMatrix(m);

            // Quaternions q and -q represent the same rotation
            double dot = Math.Abs(original.Dot(recovered));
            Assert.AreEqual(1.0, dot, 1e-10);
        }

        #endregion

        #region Axis-Angle & Euler

        [TestMethod]
        public void FromAxisAngle_ToAxisAngle_RoundTrip()
        {
            var axis = new Vector(1, 2, 3).GetUnitVector();
            double angle = 1.23;
            var q = Quaternion.FromAxisAngle(axis, angle);
            var (rAxis, rAngle) = q.ToAxisAngle();

            Assert.AreEqual(angle, rAngle, 1e-10);
            Assert.AreEqual(axis.x, rAxis.x, 1e-10);
            Assert.AreEqual(axis.y, rAxis.y, 1e-10);
            Assert.AreEqual(axis.z, rAxis.z, 1e-10);
        }

        [TestMethod]
        public void FromEulerAngles_ToEulerAngles_RoundTrip()
        {
            double roll = 0.3, pitch = 0.5, yaw = 0.7;
            var q = Quaternion.FromEulerAngles(roll, pitch, yaw);
            var (rRoll, rPitch, rYaw) = q.ToEulerAngles();

            Assert.AreEqual(roll, rRoll, 1e-10);
            Assert.AreEqual(pitch, rPitch, 1e-10);
            Assert.AreEqual(yaw, rYaw, 1e-10);
        }

        #endregion

        #region ComplexNumber Bridge

        [TestMethod]
        public void FromComplexNumber_PreservesMultiplication()
        {
            var c1 = new ComplexNumber(3, 2);
            var c2 = new ComplexNumber(1, 4);
            var cProduct = c1 * c2; // (3+2i)(1+4i) = -5+14i

            var q1 = Quaternion.FromComplexNumber(c1);
            var q2 = Quaternion.FromComplexNumber(c2);
            var qProduct = q1 * q2;

            Assert.AreEqual(cProduct.realPart, qProduct.w, Tol);
            Assert.AreEqual(cProduct.imaginaryPart, qProduct.x, Tol);
            Assert.AreEqual(0, qProduct.y, Tol);
            Assert.AreEqual(0, qProduct.z, Tol);
        }

        #endregion

        #region Interpolation

        [TestMethod]
        public void Slerp_Endpoints()
        {
            var a = Quaternion.FromAxisAngle(new Vector(0, 0, 1), 0);
            var b = Quaternion.FromAxisAngle(new Vector(0, 0, 1), Math.PI / 2);

            var at0 = Quaternion.Slerp(a, b, 0);
            var at1 = Quaternion.Slerp(a, b, 1);

            Assert.AreEqual(1.0, Math.Abs(a.Dot(at0)), 1e-10);
            Assert.AreEqual(1.0, Math.Abs(b.Dot(at1)), 1e-10);
        }

        [TestMethod]
        public void Slerp_Midpoint_IsHalfRotation()
        {
            var a = Quaternion.Identity;
            var b = Quaternion.FromAxisAngle(new Vector(0, 0, 1), Math.PI / 2);
            var mid = Quaternion.Slerp(a, b, 0.5);

            // Midpoint should be 45° rotation about Z
            var expected = Quaternion.FromAxisAngle(new Vector(0, 0, 1), Math.PI / 4);
            Assert.AreEqual(1.0, Math.Abs(mid.Dot(expected)), 1e-10);
        }

        [TestMethod]
        public void Slerp_ProducesUnitQuaternions()
        {
            var a = Quaternion.FromAxisAngle(new Vector(1, 0, 0), 0.5);
            var b = Quaternion.FromAxisAngle(new Vector(0, 1, 0), 1.5);

            for (double t = 0; t <= 1.0; t += 0.1)
            {
                var q = Quaternion.Slerp(a, b, t);
                Assert.AreEqual(1.0, q.Norm, 1e-10);
            }
        }

        #endregion

        #region Integration

        [TestMethod]
        public void IntegrateOrientation_ConstantAngularVelocity()
        {
            // Spin at 1 rad/s about Z for π/2 seconds → 90° rotation
            var omega = new Vector(0, 0, 1);
            double dt = 0.0001;
            int steps = (int)(Math.PI / 2 / dt);
            var q = Quaternion.Identity;

            for (int i = 0; i < steps; i++)
                q = Quaternion.IntegrateOrientation(q, omega, dt);

            // Should have rotated (1,0,0) → approximately (0,1,0)
            var v = q.Rotate(new Vector(1, 0, 0));
            Assert.AreEqual(0, v.x, 1e-3);
            Assert.AreEqual(1, v.y, 1e-3);
            Assert.AreEqual(0, v.z, 1e-3);
        }

        [TestMethod]
        public void IntegrateOrientation_StaysUnit()
        {
            var omega = new Vector(1, 2, 3);
            var q = Quaternion.Identity;

            for (int i = 0; i < 10000; i++)
                q = Quaternion.IntegrateOrientation(q, omega, 0.001);

            Assert.AreEqual(1.0, q.Norm, 1e-10);
        }

        #endregion
    }
}
