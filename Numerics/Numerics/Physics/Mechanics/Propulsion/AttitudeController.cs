using System;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.Mechanics.Propulsion;

/// <summary>
/// Quaternion-feedback attitude controller with PID-like gains.
/// Computes torque commands to drive the vehicle attitude toward a commanded quaternion.
/// Uses the quaternion error method for singularity-free control.
/// </summary>
public class AttitudeController
{
    /// <summary>Proportional gain (rad/s² per rad of error). Default 2.0.</summary>
    public double Kp { get; set; } = 2.0;

    /// <summary>Derivative gain (rad/s² per rad/s of rate). Default 1.5.</summary>
    public double Kd { get; set; } = 1.5;

    /// <summary>Integral gain. Default 0 (typically not needed for rocket attitude).</summary>
    public double Ki { get; set; } = 0.0;

    /// <summary>Maximum commanded angular rate (rad/s). Default 5°/s.</summary>
    public double MaxRate { get; set; } = 5.0 * Math.PI / 180.0;

    /// <summary>Maximum commanded torque magnitude (N·m). Default 100000.</summary>
    public double MaxTorque { get; set; } = 100000.0;

    private Vector _integralError = new Vector(0, 0, 0);

    /// <summary>
    /// Computes the torque command to drive current attitude toward the commanded attitude.
    /// </summary>
    /// <param name="currentAttitude">Current vehicle attitude quaternion.</param>
    /// <param name="commandedAttitude">Desired attitude quaternion.</param>
    /// <param name="angularRate">Current angular rate in body frame (rad/s).</param>
    /// <param name="dt">Time step for integral term (seconds).</param>
    /// <returns>Torque vector in body frame (N·m).</returns>
    public Vector ComputeTorque(Quaternion currentAttitude, Quaternion commandedAttitude,
        Vector angularRate, double dt)
    {
        // Quaternion error: q_err = q_commanded * q_current^-1
        Quaternion qInv = currentAttitude.Conjugate();
        Quaternion qErr = commandedAttitude * qInv;

        // Ensure short rotation path (w > 0)
        if (qErr.w < 0)
            qErr = new Quaternion(-qErr.w, -qErr.x, -qErr.y, -qErr.z);

        // Error angle in body frame (vector part of quaternion ≈ 0.5 * axis * angle for small angles)
        Vector attitudeError = new Vector(2.0 * qErr.x, 2.0 * qErr.y, 2.0 * qErr.z);

        // Integral term
        if (Ki > 0 && dt > 0)
        {
            _integralError = new Vector(
                _integralError.x + attitudeError.x * dt,
                _integralError.y + attitudeError.y * dt,
                _integralError.z + attitudeError.z * dt);

            // Anti-windup clamp
            double intMag = _integralError.GetMagnitude();
            if (intMag > 1.0)
                _integralError = (1.0 / intMag) * _integralError;
        }

        // PID torque command
        double tx = Kp * attitudeError.x - Kd * angularRate.x + Ki * _integralError.x;
        double ty = Kp * attitudeError.y - Kd * angularRate.y + Ki * _integralError.y;
        double tz = Kp * attitudeError.z - Kd * angularRate.z + Ki * _integralError.z;

        Vector torque = new Vector(tx, ty, tz);

        // Clamp to max torque
        double mag = torque.GetMagnitude();
        if (mag > MaxTorque)
            torque = (MaxTorque / mag) * torque;

        return torque;
    }

    /// <summary>
    /// Computes the desired angular rate to reach the commanded attitude (rate-limited).
    /// </summary>
    public Vector ComputeDesiredRate(Quaternion currentAttitude, Quaternion commandedAttitude)
    {
        Quaternion qInv = currentAttitude.Conjugate();
        Quaternion qErr = commandedAttitude * qInv;

        if (qErr.w < 0)
            qErr = new Quaternion(-qErr.w, -qErr.x, -qErr.y, -qErr.z);

        Vector errorAxis = new Vector(2.0 * qErr.x, 2.0 * qErr.y, 2.0 * qErr.z);
        double errorMag = errorAxis.GetMagnitude();

        // Rate-limited proportional command
        Vector desiredRate = Kp * errorAxis;
        double rateMag = desiredRate.GetMagnitude();
        if (rateMag > MaxRate)
            desiredRate = (MaxRate / rateMag) * desiredRate;

        return desiredRate;
    }

    /// <summary>Resets integral error state.</summary>
    public void Reset()
    {
        _integralError = new Vector(0, 0, 0);
    }
}
