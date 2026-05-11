using System;

namespace CSharpNumerics.Numerics.Objects;

/// <summary>
/// Utility methods for coordinate frame transformations commonly used in
/// flight dynamics and rigid-body simulation.
/// 
/// Frames:
/// • World (inertial / Earth-fixed NED or ENU)
/// • Body  (fixed to rigid body, origin at CG)
/// • Wind  (aligned with velocity vector: x = drag direction, z = lift direction)
/// 
/// Conventions follow aerospace standard: right-hand rule, ZYX Euler order.
/// </summary>
public static class FrameTransforms
{
    /// <summary>
    /// Transforms a vector from body frame to world frame using a unit quaternion.
    /// v_world = q · v_body · q*
    /// </summary>
    /// <param name="attitude">Unit quaternion representing body orientation in world frame.</param>
    /// <param name="bodyVector">Vector expressed in body-fixed coordinates.</param>
    public static Vector BodyToWorld(Quaternion attitude, Vector bodyVector)
    {
        return attitude.Rotate(bodyVector);
    }

    /// <summary>
    /// Transforms a vector from world frame to body frame using a unit quaternion.
    /// v_body = q* · v_world · q
    /// </summary>
    /// <param name="attitude">Unit quaternion representing body orientation in world frame.</param>
    /// <param name="worldVector">Vector expressed in world-fixed coordinates.</param>
    public static Vector WorldToBody(Quaternion attitude, Vector worldVector)
    {
        return attitude.Conjugate().Rotate(worldVector);
    }

    /// <summary>
    /// Transforms a vector from body frame to wind frame.
    /// The wind frame is obtained from the body frame by rotating:
    ///   1. About body Y-axis by −α (angle of attack)
    ///   2. About body Z-axis by  β (sideslip angle)
    /// 
    /// Wind frame axes: x_w = opposite to drag, y_w = side force, z_w = opposite to lift.
    /// </summary>
    /// <param name="alpha">Angle of attack in radians.</param>
    /// <param name="beta">Sideslip angle in radians.</param>
    /// <param name="bodyVector">Vector in body-frame coordinates.</param>
    public static Vector BodyToWind(double alpha, double beta, Vector bodyVector)
    {
        double ca = Math.Cos(alpha), sa = Math.Sin(alpha);
        double cb = Math.Cos(beta), sb = Math.Sin(beta);

        // Rotation matrix: R_wind←body = R_beta(Z) · R_-alpha(Y)
        double x = cb * ca * bodyVector.x + cb * sa * bodyVector.z + sb * bodyVector.y;
        double y = -sb * ca * bodyVector.x - sb * sa * bodyVector.z + cb * bodyVector.y;
        double z = -sa * bodyVector.x + ca * bodyVector.z;

        return new Vector(x, y, z);
    }

    /// <summary>
    /// Transforms a vector from wind frame to body frame.
    /// Inverse of <see cref="BodyToWind"/>.
    /// </summary>
    /// <param name="alpha">Angle of attack in radians.</param>
    /// <param name="beta">Sideslip angle in radians.</param>
    /// <param name="windVector">Vector in wind-frame coordinates.</param>
    public static Vector WindToBody(double alpha, double beta, Vector windVector)
    {
        double ca = Math.Cos(alpha), sa = Math.Sin(alpha);
        double cb = Math.Cos(beta), sb = Math.Sin(beta);

        // Transpose of R_wind←body
        double x = cb * ca * windVector.x - sb * ca * windVector.y - sa * windVector.z;
        double y = sb * windVector.x + cb * windVector.y;
        double z = cb * sa * windVector.x - sb * sa * windVector.y + ca * windVector.z;

        return new Vector(x, y, z);
    }

    /// <summary>
    /// Computes the angle of attack (α) from the body-frame velocity vector.
    /// α = atan2(w, u) where u = forward speed, w = downward speed.
    /// </summary>
    /// <param name="bodyVelocity">Velocity vector in body frame (x = forward, y = right, z = down).</param>
    public static double AngleOfAttack(Vector bodyVelocity)
    {
        return Math.Atan2(bodyVelocity.z, bodyVelocity.x);
    }

    /// <summary>
    /// Computes the sideslip angle (β) from the body-frame velocity vector.
    /// β = asin(v / |V|) where v = lateral speed.
    /// </summary>
    /// <param name="bodyVelocity">Velocity vector in body frame (x = forward, y = right, z = down).</param>
    public static double SideslipAngle(Vector bodyVelocity)
    {
        double speed = bodyVelocity.GetMagnitude();
        if (speed < 1e-12) return 0.0;
        return Math.Asin(Math.Clamp(bodyVelocity.y / speed, -1.0, 1.0));
    }

    /// <summary>
    /// Builds a 3×3 Direction Cosine Matrix (DCM) from a unit quaternion.
    /// Equivalent to <see cref="Quaternion.ToMatrix"/> but provided here for discoverability.
    /// </summary>
    /// <param name="attitude">Unit quaternion.</param>
    public static Matrix DCM(Quaternion attitude)
    {
        return attitude.ToMatrix();
    }
}
