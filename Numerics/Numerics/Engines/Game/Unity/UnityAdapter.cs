using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Engines.Game.Unity;

/// <summary>
/// Surrogate types and converters for Unity integration.
///
/// This layer does NOT depend on Unity — it defines lightweight surrogate types
/// (UnityVector3, UnityQuaternion) that mirror Unity's types. At runtime in Unity,
/// a thin shim converts between these surrogates and actual Unity types.
///
/// This allows the core library to remain dependency-free while providing
/// convenient round-trip conversions for Unity consumers.
/// </summary>
public static class UnityAdapter
{
    // ════════════════════════════════════════════════════════════
    //  Surrogate Types
    // ════════════════════════════════════════════════════════════

    /// <summary>
    /// Surrogate for Unity's Vector3 (float precision, Y-up convention).
    /// Unity uses left-handed Y-up; CSharpNumerics uses right-handed Z-up.
    /// </summary>
    public struct UnityVector3
    {
        public float x, y, z;

        public UnityVector3(float x, float y, float z)
        {
            this.x = x; this.y = y; this.z = z;
        }

        public override string ToString() => $"({x:F3}, {y:F3}, {z:F3})";
    }

    /// <summary>
    /// Surrogate for Unity's Quaternion (float precision).
    /// </summary>
    public struct UnityQuaternion
    {
        public float x, y, z, w;

        public UnityQuaternion(float x, float y, float z, float w)
        {
            this.x = x; this.y = y; this.z = z; this.w = w;
        }

        /// <summary>Identity quaternion (no rotation).</summary>
        public static UnityQuaternion Identity => new UnityQuaternion(0, 0, 0, 1);

        public override string ToString() => $"({x:F3}, {y:F3}, {z:F3}, {w:F3})";
    }

    // ════════════════════════════════════════════════════════════
    //  Coordinate Convention Converters
    //  CSharpNumerics: right-handed, Z-up (X-right, Y-forward, Z-up)
    //  Unity:          left-handed, Y-up  (X-right, Y-up, Z-forward)
    //  Mapping: CSN(x,y,z) → Unity(x,z,y) with Z-flip for handedness
    // ════════════════════════════════════════════════════════════

    /// <summary>
    /// Convert CSharpNumerics Vector (Z-up) to Unity Vector3 (Y-up).
    /// CSN(x, y, z) → Unity(x, z, y).
    /// </summary>
    public static UnityVector3 ToUnityVector3(Vector v)
    {
        return new UnityVector3((float)v.x, (float)v.z, (float)v.y);
    }

    /// <summary>
    /// Convert Unity Vector3 (Y-up) to CSharpNumerics Vector (Z-up).
    /// Unity(x, y, z) → CSN(x, z, y).
    /// </summary>
    public static Vector FromUnityVector3(UnityVector3 v)
    {
        return new Vector(v.x, v.z, v.y);
    }

    /// <summary>
    /// Convert CSharpNumerics 3×3 rotation matrix (Z-up) to Unity quaternion (Y-up).
    /// </summary>
    public static UnityQuaternion ToUnityQuaternion(Matrix rotation)
    {
        if (rotation.rowLength != 3 || rotation.columnLength != 3)
            throw new ArgumentException("Rotation matrix must be 3×3.");

        // Convert rotation matrix from Z-up to Y-up coordinate system
        // by swapping Y and Z axes
        var swapped = new double[3, 3];
        // Map axes: X'=X, Y'=Z, Z'=Y
        // New matrix R' = S * R * S^T where S swaps Y,Z
        swapped[0, 0] = rotation.values[0, 0];
        swapped[0, 1] = rotation.values[0, 2];
        swapped[0, 2] = rotation.values[0, 1];
        swapped[1, 0] = rotation.values[2, 0];
        swapped[1, 1] = rotation.values[2, 2];
        swapped[1, 2] = rotation.values[2, 1];
        swapped[2, 0] = rotation.values[1, 0];
        swapped[2, 1] = rotation.values[1, 2];
        swapped[2, 2] = rotation.values[1, 1];

        return MatrixToQuaternion(swapped);
    }

    /// <summary>
    /// Convert Unity quaternion (Y-up) to CSharpNumerics 3×3 rotation matrix (Z-up).
    /// </summary>
    public static Matrix FromUnityQuaternion(UnityQuaternion q)
    {
        // Convert quaternion to rotation matrix in Y-up
        var yUpMatrix = QuaternionToMatrix(q.x, q.y, q.z, q.w);

        // Swap Y and Z back to Z-up
        var result = new double[3, 3];
        result[0, 0] = yUpMatrix[0, 0];
        result[0, 1] = yUpMatrix[0, 2];
        result[0, 2] = yUpMatrix[0, 1];
        result[1, 0] = yUpMatrix[2, 0];
        result[1, 1] = yUpMatrix[2, 2];
        result[1, 2] = yUpMatrix[2, 1];
        result[2, 0] = yUpMatrix[1, 0];
        result[2, 1] = yUpMatrix[1, 2];
        result[2, 2] = yUpMatrix[1, 1];

        return new Matrix(result);
    }

    /// <summary>
    /// Convert CSharpNumerics VectorN to Unity Vector3 (takes first 3 components).
    /// No coordinate swap — direct mapping for abstract data.
    /// </summary>
    public static UnityVector3 VectorNToUnityVector3(VectorN v)
    {
        return new UnityVector3(
            v.Length > 0 ? (float)v[0] : 0,
            v.Length > 1 ? (float)v[1] : 0,
            v.Length > 2 ? (float)v[2] : 0);
    }

    /// <summary>
    /// Convert Unity Vector3 to CSharpNumerics VectorN (3 components).
    /// No coordinate swap — direct mapping for abstract data.
    /// </summary>
    public static VectorN UnityVector3ToVectorN(UnityVector3 v)
    {
        return new VectorN(new double[] { v.x, v.y, v.z });
    }

    // ════════════════════════════════════════════════════════════
    //  Helpers
    // ════════════════════════════════════════════════════════════

    private static UnityQuaternion MatrixToQuaternion(double[,] m)
    {
        double trace = m[0, 0] + m[1, 1] + m[2, 2];
        double x, y, z, w;

        if (trace > 0)
        {
            double s = 0.5 / Math.Sqrt(trace + 1.0);
            w = 0.25 / s;
            x = (m[2, 1] - m[1, 2]) * s;
            y = (m[0, 2] - m[2, 0]) * s;
            z = (m[1, 0] - m[0, 1]) * s;
        }
        else if (m[0, 0] > m[1, 1] && m[0, 0] > m[2, 2])
        {
            double s = 2.0 * Math.Sqrt(1.0 + m[0, 0] - m[1, 1] - m[2, 2]);
            w = (m[2, 1] - m[1, 2]) / s;
            x = 0.25 * s;
            y = (m[0, 1] + m[1, 0]) / s;
            z = (m[0, 2] + m[2, 0]) / s;
        }
        else if (m[1, 1] > m[2, 2])
        {
            double s = 2.0 * Math.Sqrt(1.0 + m[1, 1] - m[0, 0] - m[2, 2]);
            w = (m[0, 2] - m[2, 0]) / s;
            x = (m[0, 1] + m[1, 0]) / s;
            y = 0.25 * s;
            z = (m[1, 2] + m[2, 1]) / s;
        }
        else
        {
            double s = 2.0 * Math.Sqrt(1.0 + m[2, 2] - m[0, 0] - m[1, 1]);
            w = (m[1, 0] - m[0, 1]) / s;
            x = (m[0, 2] + m[2, 0]) / s;
            y = (m[1, 2] + m[2, 1]) / s;
            z = 0.25 * s;
        }

        // Normalize
        double mag = Math.Sqrt(x * x + y * y + z * z + w * w);
        if (mag > 1e-10)
        {
            x /= mag; y /= mag; z /= mag; w /= mag;
        }

        return new UnityQuaternion((float)x, (float)y, (float)z, (float)w);
    }

    private static double[,] QuaternionToMatrix(double x, double y, double z, double w)
    {
        double xx = x * x, yy = y * y, zz = z * z;
        double xy = x * y, xz = x * z, yz = y * z;
        double wx = w * x, wy = w * y, wz = w * z;

        return new double[,]
        {
            { 1 - 2 * (yy + zz), 2 * (xy - wz),     2 * (xz + wy) },
            { 2 * (xy + wz),     1 - 2 * (xx + zz),  2 * (yz - wx) },
            { 2 * (xz - wy),     2 * (yz + wx),      1 - 2 * (xx + yy) }
        };
    }
}
