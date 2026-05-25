using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Mechanics.Propulsion;

/// <summary>
/// Variable-mass dynamics for rocket flight.
/// Extends rigid-body equations of motion with the dm/dt term:
///   F_net = m·a + ṁ·v_exhaust  (rocket equation in non-inertial form)
/// 
/// In inertial frame: F_thrust = ṁ · Ve (exhaust velocity),
/// so net force = thrust + external forces, a = F_net / m(t).
/// </summary>
public static class VariableMassDynamics
{
    private const double G0 = 9.80665;

    /// <summary>
    /// Computes linear acceleration for a variable-mass body.
    /// a = (F_thrust + F_external) / m
    /// where F_thrust = ṁ · v_exhaust = ṁ · Isp · g₀ (in thrust direction).
    /// </summary>
    /// <param name="totalMass">Current total mass in kg.</param>
    /// <param name="thrustForce">Thrust force vector in world frame (N).</param>
    /// <param name="externalForces">Sum of all external forces: gravity, drag, etc. (N).</param>
    /// <returns>Acceleration vector in world frame (m/s²).</returns>
    public static Vector ComputeAcceleration(double totalMass, Vector thrustForce, Vector externalForces)
    {
        if (totalMass <= 0) throw new ArgumentOutOfRangeException(nameof(totalMass));
        Vector netForce = thrustForce + externalForces;
        return (1.0 / totalMass) * netForce;
    }

    /// <summary>
    /// Computes angular acceleration for a body with changing inertia.
    /// Euler's rotational equation: τ = I·α + ω × (I·ω)
    /// α = I⁻¹ · (τ - ω × (I·ω))
    /// </summary>
    /// <param name="inertiaTensor">Current 3×3 inertia tensor in body frame.</param>
    /// <param name="torque">Net torque vector in body frame (N·m).</param>
    /// <param name="angularVelocity">Angular velocity in body frame (rad/s).</param>
    /// <returns>Angular acceleration in body frame (rad/s²).</returns>
    public static Vector ComputeAngularAcceleration(Matrix inertiaTensor, Vector torque, Vector angularVelocity)
    {
        Matrix invI = inertiaTensor.Inverse();
        // I·ω
        Vector Iw = MultiplyMatrixVector(inertiaTensor, angularVelocity);
        // ω × (I·ω)
        Vector gyroscopic = angularVelocity.Cross(Iw);
        // α = I⁻¹ · (τ - ω × I·ω)
        Vector rhs = torque - gyroscopic;
        return MultiplyMatrixVector(invI, rhs);
    }

    /// <summary>
    /// Computes the thrust moment about center of mass from gimbaled engine.
    /// τ = r × F where r is the vector from CG to engine and F is thrust force in body frame.
    /// </summary>
    /// <param name="cgToEngine">Vector from CG to engine mount point (body frame, m).</param>
    /// <param name="thrustBody">Thrust force vector in body frame (N).</param>
    public static Vector ComputeThrustTorque(Vector cgToEngine, Vector thrustBody)
    {
        return cgToEngine.Cross(thrustBody);
    }

    private static Vector MultiplyMatrixVector(Matrix m, Vector v)
    {
        double x = m.values[0, 0] * v.x + m.values[0, 1] * v.y + m.values[0, 2] * v.z;
        double y = m.values[1, 0] * v.x + m.values[1, 1] * v.y + m.values[1, 2] * v.z;
        double z = m.values[2, 0] * v.x + m.values[2, 1] * v.y + m.values[2, 2] * v.z;
        return new Vector(x, y, z);
    }
}
