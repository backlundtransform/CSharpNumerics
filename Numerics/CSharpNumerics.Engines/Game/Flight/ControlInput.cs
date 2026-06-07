using System;

namespace CSharpNumerics.Engines.Game.Flight;

/// <summary>
/// Represents pilot/autopilot control inputs to the aircraft.
/// All deflections are normalized to [−1, +1] and mapped to physical
/// limits by <see cref="AircraftConfig"/>.
/// </summary>
public class ControlInput
{
    /// <summary>Throttle setting (0 = idle, 1 = full power).</summary>
    public double Throttle { get; set; }

    /// <summary>Pitch command (−1 = full nose-down, +1 = full nose-up). Maps to elevator deflection.</summary>
    public double Pitch { get; set; }

    /// <summary>Roll command (−1 = full left, +1 = full right). Maps to aileron deflection.</summary>
    public double Roll { get; set; }

    /// <summary>Yaw command (−1 = full left, +1 = full right). Maps to rudder deflection.</summary>
    public double Yaw { get; set; }

    /// <summary>Flap setting (0 = retracted, 1 = fully extended).</summary>
    public double Flaps { get; set; }

    /// <summary>Landing gear state (true = down).</summary>
    public bool GearDown { get; set; }

    /// <summary>Creates a neutral (zero-input) control state.</summary>
    public ControlInput()
    {
    }

    /// <summary>Creates a control input with specified values.</summary>
    public ControlInput(double throttle, double pitch, double roll, double yaw,
        double flaps = 0, bool gearDown = false)
    {
        Throttle = Math.Clamp(throttle, 0, 1);
        Pitch = Math.Clamp(pitch, -1, 1);
        Roll = Math.Clamp(roll, -1, 1);
        Yaw = Math.Clamp(yaw, -1, 1);
        Flaps = Math.Clamp(flaps, 0, 1);
        GearDown = gearDown;
    }
}
