using System;

namespace CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

/// <summary>
/// Models the aerodynamic effect of a control surface (elevator, aileron, rudder, flap).
/// A deflection δ produces incremental changes in lift and drag coefficients:
///   ΔCl = τ · Cl_δ · δ
///   ΔCd = Cd_δ · δ²
/// where τ is the flap effectiveness factor (function of chord ratio).
/// </summary>
public class ControlSurface
{
    /// <summary>Name of the control surface (e.g., "Elevator", "Left Aileron").</summary>
    public string Name { get; }

    /// <summary>
    /// Lift coefficient sensitivity to deflection in 1/rad.
    /// Typical values: elevator ≈ 3.0–5.0, aileron ≈ 3.0–4.5.
    /// </summary>
    public double ClDelta { get; }

    /// <summary>
    /// Drag coefficient sensitivity to deflection squared in 1/rad².
    /// Typical values: 0.01–0.05.
    /// </summary>
    public double CdDelta { get; }

    /// <summary>
    /// Flap effectiveness factor τ (0–1). Depends on the chord ratio of the control surface.
    /// τ = 1 means the full theoretical effect; typical real values ≈ 0.4–0.7.
    /// </summary>
    public double Tau { get; }

    /// <summary>Maximum positive deflection in radians.</summary>
    public double MaxDeflection { get; }

    /// <summary>Minimum (negative) deflection in radians.</summary>
    public double MinDeflection { get; }

    /// <summary>
    /// Which body axis this control surface primarily affects.
    /// Elevator → Pitch (Y), Aileron → Roll (X), Rudder → Yaw (Z).
    /// </summary>
    public ControlAxis Axis { get; }

    /// <summary>
    /// Creates a control surface model.
    /// </summary>
    /// <param name="name">Surface name.</param>
    /// <param name="axis">Primary control axis.</param>
    /// <param name="clDelta">Lift sensitivity dCl/dδ in 1/rad.</param>
    /// <param name="cdDelta">Drag sensitivity dCd/dδ² in 1/rad².</param>
    /// <param name="tau">Flap effectiveness factor (0–1).</param>
    /// <param name="maxDeflection">Max positive deflection in radians.</param>
    /// <param name="minDeflection">Min (negative) deflection in radians.</param>
    public ControlSurface(
        string name,
        ControlAxis axis,
        double clDelta,
        double cdDelta,
        double tau = 0.55,
        double maxDeflection = 0.4363,   // 25°
        double minDeflection = -0.4363)
    {
        Name = name;
        Axis = axis;
        ClDelta = clDelta;
        CdDelta = cdDelta;
        Tau = Math.Clamp(tau, 0, 1);
        MaxDeflection = maxDeflection;
        MinDeflection = minDeflection;
    }

    /// <summary>
    /// Computes the incremental lift coefficient from a control surface deflection.
    /// ΔCl = τ · Cl_δ · δ_clamped
    /// </summary>
    /// <param name="deflection">Commanded deflection in radians (will be clamped to limits).</param>
    public double DeltaCl(double deflection)
    {
        double d = ClampDeflection(deflection);
        return Tau * ClDelta * d;
    }

    /// <summary>
    /// Computes the incremental drag coefficient from a control surface deflection.
    /// ΔCd = Cd_δ · δ²
    /// </summary>
    /// <param name="deflection">Commanded deflection in radians (will be clamped to limits).</param>
    public double DeltaCd(double deflection)
    {
        double d = ClampDeflection(deflection);
        return CdDelta * d * d;
    }

    /// <summary>
    /// Clamps a deflection command to the mechanical limits of this surface.
    /// </summary>
    public double ClampDeflection(double deflection)
    {
        return Math.Clamp(deflection, MinDeflection, MaxDeflection);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Factory methods for common surfaces
    // ═══════════════════════════════════════════════════════════════

    /// <summary>Creates a typical elevator control surface.</summary>
    public static ControlSurface Elevator(
        double clDelta = 4.0, double cdDelta = 0.02, double tau = 0.55)
        => new ControlSurface("Elevator", ControlAxis.Pitch, clDelta, cdDelta, tau);

    /// <summary>Creates a typical aileron control surface.</summary>
    public static ControlSurface Aileron(
        double clDelta = 3.5, double cdDelta = 0.015, double tau = 0.50)
        => new ControlSurface("Aileron", ControlAxis.Roll, clDelta, cdDelta, tau);

    /// <summary>Creates a typical rudder control surface.</summary>
    public static ControlSurface Rudder(
        double clDelta = 3.0, double cdDelta = 0.02, double tau = 0.50)
        => new ControlSurface("Rudder", ControlAxis.Yaw, clDelta, cdDelta, tau);
}

/// <summary>
/// Identifies the primary body axis a control surface affects.
/// </summary>
public enum ControlAxis
{
    /// <summary>Roll axis (body X).</summary>
    Roll,
    /// <summary>Pitch axis (body Y).</summary>
    Pitch,
    /// <summary>Yaw axis (body Z).</summary>
    Yaw
}
