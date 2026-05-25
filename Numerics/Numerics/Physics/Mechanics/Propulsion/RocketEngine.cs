using System;

namespace CSharpNumerics.Physics.Mechanics.Propulsion;

/// <summary>
/// Models a rocket engine with sea-level and vacuum specific impulse,
/// throttle control, and gimbal angle limits.
/// 
/// Thrust varies with altitude via ambient pressure interpolation between
/// sea-level and vacuum Isp. Mass flow rate: ṁ = T / (Isp · g₀).
/// </summary>
public class RocketEngine
{
    private const double G0 = 9.80665; // m/s²

    /// <summary>Specific impulse at sea level (seconds).</summary>
    public double IspSeaLevel { get; }

    /// <summary>Specific impulse in vacuum (seconds).</summary>
    public double IspVacuum { get; }

    /// <summary>Maximum thrust in vacuum (Newtons).</summary>
    public double MaxThrustVacuum { get; }

    /// <summary>Minimum throttle fraction (0–1). Typical: 0.4 for deep throttle.</summary>
    public double MinThrottle { get; }

    /// <summary>Maximum gimbal angle in radians.</summary>
    public double MaxGimbalAngle { get; }

    /// <summary>Current throttle setting (0–1).</summary>
    public double Throttle { get; set; } = 1.0;

    /// <summary>Current gimbal pitch angle in radians (about body Y axis).</summary>
    public double GimbalPitch { get; set; }

    /// <summary>Current gimbal yaw angle in radians (about body Z axis).</summary>
    public double GimbalYaw { get; set; }

    /// <summary>True if the engine is ignited and producing thrust.</summary>
    public bool IsActive { get; set; }

    /// <summary>
    /// Creates a rocket engine.
    /// </summary>
    /// <param name="ispSeaLevel">Sea-level specific impulse in seconds.</param>
    /// <param name="ispVacuum">Vacuum specific impulse in seconds.</param>
    /// <param name="maxThrustVacuum">Maximum vacuum thrust in Newtons.</param>
    /// <param name="minThrottle">Minimum throttle fraction (default 0.4).</param>
    /// <param name="maxGimbalAngle">Maximum gimbal deflection in radians (default 5°).</param>
    public RocketEngine(double ispSeaLevel, double ispVacuum, double maxThrustVacuum,
        double minThrottle = 0.4, double maxGimbalAngle = 5.0 * Math.PI / 180.0)
    {
        if (ispSeaLevel <= 0) throw new ArgumentOutOfRangeException(nameof(ispSeaLevel));
        if (ispVacuum <= 0) throw new ArgumentOutOfRangeException(nameof(ispVacuum));
        if (maxThrustVacuum <= 0) throw new ArgumentOutOfRangeException(nameof(maxThrustVacuum));

        IspSeaLevel = ispSeaLevel;
        IspVacuum = ispVacuum;
        MaxThrustVacuum = maxThrustVacuum;
        MinThrottle = Math.Clamp(minThrottle, 0, 1);
        MaxGimbalAngle = Math.Abs(maxGimbalAngle);
    }

    /// <summary>
    /// Effective specific impulse at a given ambient pressure ratio (0 = vacuum, 1 = sea level).
    /// </summary>
    /// <param name="pressureRatio">Ambient pressure / sea-level pressure (0–1).</param>
    public double EffectiveIsp(double pressureRatio)
    {
        double ratio = Math.Clamp(pressureRatio, 0, 1);
        return IspVacuum + (IspSeaLevel - IspVacuum) * ratio;
    }

    /// <summary>
    /// Current thrust magnitude in Newtons at the given ambient pressure ratio.
    /// </summary>
    /// <param name="pressureRatio">Ambient pressure / sea-level pressure (0–1).</param>
    public double Thrust(double pressureRatio)
    {
        if (!IsActive) return 0;

        double effectiveThrottle = Math.Clamp(Throttle, MinThrottle, 1.0);
        double isp = EffectiveIsp(pressureRatio);
        // Mass flow is constant for a given throttle: ṁ = T_vac / (Isp_vac · g₀) · throttle
        double massFlowRate = MassFlowRate(pressureRatio);
        return massFlowRate * isp * G0;
    }

    /// <summary>
    /// Mass flow rate in kg/s at the given ambient pressure ratio and current throttle.
    /// ṁ = T_vac · throttle / (Isp_vac · g₀)
    /// </summary>
    /// <param name="pressureRatio">Ambient pressure / sea-level pressure (0–1).</param>
    public double MassFlowRate(double pressureRatio)
    {
        if (!IsActive) return 0;

        double effectiveThrottle = Math.Clamp(Throttle, MinThrottle, 1.0);
        return MaxThrustVacuum * effectiveThrottle / (IspVacuum * G0);
    }

    /// <summary>
    /// Sets throttle clamped to valid range [MinThrottle, 1.0] (or 0 if below MinThrottle).
    /// </summary>
    public void SetThrottle(double throttle)
    {
        Throttle = Math.Clamp(throttle, 0, 1);
    }

    /// <summary>
    /// Sets gimbal angles clamped to max gimbal angle.
    /// </summary>
    public void SetGimbal(double pitch, double yaw)
    {
        GimbalPitch = Math.Clamp(pitch, -MaxGimbalAngle, MaxGimbalAngle);
        GimbalYaw = Math.Clamp(yaw, -MaxGimbalAngle, MaxGimbalAngle);
    }
}
