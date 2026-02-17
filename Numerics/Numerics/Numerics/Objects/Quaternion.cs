using System;

namespace Numerics.Objects
{
    /// <summary>
    /// Represents a quaternion q = w + xi + yj + zk.
    /// Quaternions generalize complex numbers to 4 dimensions and are the
    /// standard representation for 3D rotations — compact, numerically stable,
    /// and free of gimbal lock.
    /// 
    /// Algebraic hierarchy:  ℝ (double) ⊂ ℂ (ComplexNumber) ⊂ ℍ (Quaternion)
    /// </summary>
    public struct Quaternion
    {
        /// <summary>Scalar (real) part.</summary>
        public double w;

        /// <summary>First imaginary component (i).</summary>
        public double x;

        /// <summary>Second imaginary component (j).</summary>
        public double y;

        /// <summary>Third imaginary component (k).</summary>
        public double z;

        /// <summary>
        /// Creates a quaternion from its four components.
        /// </summary>
        public Quaternion(double w, double x, double y, double z)
        {
            this.w = w;
            this.x = x;
            this.y = y;
            this.z = z;
        }

        #region Static Factories

        /// <summary>Identity quaternion (no rotation): (1, 0, 0, 0).</summary>
        public static Quaternion Identity => new Quaternion(1, 0, 0, 0);

        /// <summary>Zero quaternion: (0, 0, 0, 0).</summary>
        public static Quaternion Zero => new Quaternion(0, 0, 0, 0);

        /// <summary>
        /// Creates a unit quaternion from an axis and angle.
        /// q = cos(θ/2) + sin(θ/2)·(ux·i + uy·j + uz·k).
        /// </summary>
        /// <param name="axis">Rotation axis (will be normalized).</param>
        /// <param name="angleRadians">Rotation angle in radians.</param>
        public static Quaternion FromAxisAngle(Vector axis, double angleRadians)
        {
            var mag = axis.GetMagnitude();
            if (mag < 1e-15) return Identity;

            double halfAngle = angleRadians * 0.5;
            double s = Math.Sin(halfAngle) / mag;
            return new Quaternion(Math.Cos(halfAngle), axis.x * s, axis.y * s, axis.z * s);
        }

        /// <summary>
        /// Creates a unit quaternion from Euler angles (ZYX convention: yaw → pitch → roll).
        /// </summary>
        /// <param name="roll">Rotation about X in radians.</param>
        /// <param name="pitch">Rotation about Y in radians.</param>
        /// <param name="yaw">Rotation about Z in radians.</param>
        public static Quaternion FromEulerAngles(double roll, double pitch, double yaw)
        {
            double cr = Math.Cos(roll * 0.5), sr = Math.Sin(roll * 0.5);
            double cp = Math.Cos(pitch * 0.5), sp = Math.Sin(pitch * 0.5);
            double cy = Math.Cos(yaw * 0.5), sy = Math.Sin(yaw * 0.5);

            return new Quaternion(
                cr * cp * cy + sr * sp * sy,
                sr * cp * cy - cr * sp * sy,
                cr * sp * cy + sr * cp * sy,
                cr * cp * sy - sr * sp * cy);
        }

        /// <summary>
        /// Embeds a complex number into the quaternion algebra: a + bi → (a, b, 0, 0).
        /// This preserves complex multiplication as a subgroup of ℍ.
        /// </summary>
        /// <param name="c">The complex number to embed.</param>
        public static Quaternion FromComplexNumber(ComplexNumber c)
        {
            return new Quaternion(c.realPart, c.imaginaryPart, 0, 0);
        }

        /// <summary>
        /// Extracts a unit quaternion from a 3×3 rotation matrix using Shepperd's method.
        /// </summary>
        /// <param name="m">A 3×3 rotation matrix.</param>
        public static Quaternion FromMatrix(Matrix m)
        {
            double m00 = m.values[0, 0], m11 = m.values[1, 1], m22 = m.values[2, 2];
            double trace = m00 + m11 + m22;

            double w, x, y, z;

            if (trace > 0)
            {
                double s = 0.5 / Math.Sqrt(trace + 1.0);
                w = 0.25 / s;
                x = (m.values[2, 1] - m.values[1, 2]) * s;
                y = (m.values[0, 2] - m.values[2, 0]) * s;
                z = (m.values[1, 0] - m.values[0, 1]) * s;
            }
            else if (m00 > m11 && m00 > m22)
            {
                double s = 2.0 * Math.Sqrt(1.0 + m00 - m11 - m22);
                w = (m.values[2, 1] - m.values[1, 2]) / s;
                x = 0.25 * s;
                y = (m.values[0, 1] + m.values[1, 0]) / s;
                z = (m.values[0, 2] + m.values[2, 0]) / s;
            }
            else if (m11 > m22)
            {
                double s = 2.0 * Math.Sqrt(1.0 + m11 - m00 - m22);
                w = (m.values[0, 2] - m.values[2, 0]) / s;
                x = (m.values[0, 1] + m.values[1, 0]) / s;
                y = 0.25 * s;
                z = (m.values[1, 2] + m.values[2, 1]) / s;
            }
            else
            {
                double s = 2.0 * Math.Sqrt(1.0 + m22 - m00 - m11);
                w = (m.values[1, 0] - m.values[0, 1]) / s;
                x = (m.values[0, 2] + m.values[2, 0]) / s;
                y = (m.values[1, 2] + m.values[2, 1]) / s;
                z = 0.25 * s;
            }

            return new Quaternion(w, x, y, z).Normalize();
        }

        #endregion

        #region Properties

        /// <summary>The Euclidean norm: |q| = √(w² + x² + y² + z²).</summary>
        public double Norm => Math.Sqrt(w * w + x * x + y * y + z * z);

        /// <summary>Squared norm (avoids the square root).</summary>
        public double NormSquared => w * w + x * x + y * y + z * z;

        /// <summary>True if this is a unit quaternion (|q| ≈ 1).</summary>
        public bool IsUnit => Math.Abs(NormSquared - 1.0) < 1e-12;

        /// <summary>The vector (imaginary) part as a Vector.</summary>
        public Vector VectorPart => new Vector(x, y, z);

        #endregion

        #region Core Operations

        /// <summary>
        /// Returns the conjugate: q* = w - xi - yj - zk.
        /// For unit quaternions, the conjugate equals the inverse.
        /// </summary>
        public Quaternion Conjugate() => new Quaternion(w, -x, -y, -z);

        /// <summary>
        /// Returns the normalized (unit) quaternion: q / |q|.
        /// </summary>
        public Quaternion Normalize()
        {
            double n = Norm;
            if (n < 1e-15) return Identity;
            double inv = 1.0 / n;
            return new Quaternion(w * inv, x * inv, y * inv, z * inv);
        }

        /// <summary>
        /// Returns the multiplicative inverse: q⁻¹ = q* / |q|².
        /// </summary>
        public Quaternion Inverse()
        {
            double ns = NormSquared;
            if (ns < 1e-30) throw new InvalidOperationException("Cannot invert a zero quaternion.");
            double inv = 1.0 / ns;
            return new Quaternion(w * inv, -x * inv, -y * inv, -z * inv);
        }

        /// <summary>
        /// Dot product: q₁·q₂ = w₁w₂ + x₁x₂ + y₁y₂ + z₁z₂.
        /// Useful for measuring similarity between rotations.
        /// </summary>
        public double Dot(Quaternion other) => w * other.w + x * other.x + y * other.y + z * other.z;

        #endregion

        #region Rotation

        /// <summary>
        /// Rotates a vector by this unit quaternion: v' = q·v·q*.
        /// Uses the optimized formula to avoid full quaternion multiplication.
        /// </summary>
        /// <param name="v">The vector to rotate.</param>
        public Vector Rotate(Vector v)
        {
            // t = 2(q_vec × v)
            double tx = 2.0 * (y * v.z - z * v.y);
            double ty = 2.0 * (z * v.x - x * v.z);
            double tz = 2.0 * (x * v.y - y * v.x);

            // result = v + w·t + q_vec × t
            return new Vector(
                v.x + w * tx + (y * tz - z * ty),
                v.y + w * ty + (z * tx - x * tz),
                v.z + w * tz + (x * ty - y * tx));
        }

        /// <summary>
        /// Converts this unit quaternion to a 3×3 rotation matrix.
        /// </summary>
        public Matrix ToMatrix()
        {
            double xx = x * x, yy = y * y, zz = z * z;
            double xy = x * y, xz = x * z, yz = y * z;
            double wx = w * x, wy = w * y, wz = w * z;

            return new Matrix(new double[,]
            {
                { 1 - 2 * (yy + zz),     2 * (xy - wz),     2 * (xz + wy) },
                {     2 * (xy + wz), 1 - 2 * (xx + zz),     2 * (yz - wx) },
                {     2 * (xz - wy),     2 * (yz + wx), 1 - 2 * (xx + yy) }
            });
        }

        /// <summary>
        /// Extracts the rotation axis and angle from this unit quaternion.
        /// Returns (axis, angleRadians). If the angle is near zero, the axis defaults to (0, 0, 1).
        /// </summary>
        public (Vector axis, double angle) ToAxisAngle()
        {
            var q = NormSquared > 0 ? Normalize() : this;
            double angle = 2.0 * Math.Acos(Math.Clamp(q.w, -1.0, 1.0));
            double s = Math.Sqrt(1.0 - q.w * q.w);
            if (s < 1e-10)
                return (new Vector(0, 0, 1), angle);
            return (new Vector(q.x / s, q.y / s, q.z / s), angle);
        }

        /// <summary>
        /// Extracts Euler angles (ZYX convention) from this unit quaternion.
        /// Returns (roll, pitch, yaw) in radians.
        /// </summary>
        public (double roll, double pitch, double yaw) ToEulerAngles()
        {
            // Roll (X)
            double sinr = 2.0 * (w * x + y * z);
            double cosr = 1.0 - 2.0 * (x * x + y * y);
            double roll = Math.Atan2(sinr, cosr);

            // Pitch (Y) — clamp for safety
            double sinp = 2.0 * (w * y - z * x);
            double pitch = Math.Abs(sinp) >= 1.0
                ? Math.CopySign(Math.PI / 2, sinp)
                : Math.Asin(sinp);

            // Yaw (Z)
            double siny = 2.0 * (w * z + x * y);
            double cosy = 1.0 - 2.0 * (y * y + z * z);
            double yaw = Math.Atan2(siny, cosy);

            return (roll, pitch, yaw);
        }

        #endregion

        #region Interpolation

        /// <summary>
        /// Spherical linear interpolation between two unit quaternions.
        /// Follows the shortest arc (flips sign if dot product is negative).
        /// </summary>
        /// <param name="a">Start rotation (t = 0).</param>
        /// <param name="b">End rotation (t = 1).</param>
        /// <param name="t">Interpolation parameter [0, 1].</param>
        public static Quaternion Slerp(Quaternion a, Quaternion b, double t)
        {
            double dot = a.Dot(b);

            // Ensure shortest path
            if (dot < 0)
            {
                b = -b;
                dot = -dot;
            }

            // Fall back to lerp for nearly parallel quaternions
            if (dot > 0.9995)
                return Lerp(a, dot < 0 ? -b : b, t);

            double theta = Math.Acos(Math.Clamp(dot, -1.0, 1.0));
            double sinTheta = Math.Sin(theta);
            double wa = Math.Sin((1.0 - t) * theta) / sinTheta;
            double wb = Math.Sin(t * theta) / sinTheta;

            return new Quaternion(
                wa * a.w + wb * b.w,
                wa * a.x + wb * b.x,
                wa * a.y + wb * b.y,
                wa * a.z + wb * b.z);
        }

        /// <summary>
        /// Normalized linear interpolation — cheaper than Slerp, nearly identical for small angles.
        /// </summary>
        /// <param name="a">Start rotation.</param>
        /// <param name="b">End rotation.</param>
        /// <param name="t">Interpolation parameter [0, 1].</param>
        public static Quaternion Lerp(Quaternion a, Quaternion b, double t)
        {
            return new Quaternion(
                a.w + t * (b.w - a.w),
                a.x + t * (b.x - a.x),
                a.y + t * (b.y - a.y),
                a.z + t * (b.z - a.z)).Normalize();
        }

        #endregion

        #region Integration

        /// <summary>
        /// Integrates orientation by one time step using angular velocity.
        /// q(t+dt) = normalize(q + ½·dt·ω·q), where ω is treated as a pure quaternion (0, ωx, ωy, ωz).
        /// </summary>
        /// <param name="q">Current orientation quaternion.</param>
        /// <param name="angularVelocity">Angular velocity vector in rad/s.</param>
        /// <param name="dt">Time step in seconds.</param>
        public static Quaternion IntegrateOrientation(Quaternion q, Vector angularVelocity, double dt)
        {
            var omegaQ = new Quaternion(0, angularVelocity.x, angularVelocity.y, angularVelocity.z);
            var qdot = 0.5 * (omegaQ * q);
            return new Quaternion(
                q.w + qdot.w * dt,
                q.x + qdot.x * dt,
                q.y + qdot.y * dt,
                q.z + qdot.z * dt).Normalize();
        }

        #endregion

        #region Operators

        /// <summary>Hamilton product: q₁ · q₂ (non-commutative).</summary>
        public static Quaternion operator *(Quaternion a, Quaternion b)
        {
            return new Quaternion(
                a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
                a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
                a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
                a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w);
        }

        public static Quaternion operator +(Quaternion a, Quaternion b)
            => new Quaternion(a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z);

        public static Quaternion operator -(Quaternion a, Quaternion b)
            => new Quaternion(a.w - b.w, a.x - b.x, a.y - b.y, a.z - b.z);

        public static Quaternion operator -(Quaternion q)
            => new Quaternion(-q.w, -q.x, -q.y, -q.z);

        public static Quaternion operator *(double s, Quaternion q)
            => new Quaternion(s * q.w, s * q.x, s * q.y, s * q.z);

        public static Quaternion operator *(Quaternion q, double s)
            => new Quaternion(q.w * s, q.x * s, q.y * s, q.z * s);

        public static Quaternion operator /(Quaternion q, double s)
            => new Quaternion(q.w / s, q.x / s, q.y / s, q.z / s);

        #endregion

        public override string ToString() => $"({w}, {x}i, {y}j, {z}k)";
    }
}
