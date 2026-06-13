using System;
using CSharpNumerics.Engines.Game.Rocket;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Engines.Game.Rocket.Examples;

/// <summary>
/// Example: Telemetry HUD prefab configuration.
/// 
/// Shows how to wire the TelemetryStream into a Unity UI display showing:
/// - Altitude (m / km)
/// - Velocity (m/s)
/// - Acceleration (m/s² / G)
/// - Fuel remaining (%)
/// - Apoapsis / Periapsis (km)
/// - Time to apoapsis
/// - Current mission phase
///
/// <code>
/// // In Unity MonoBehaviour:
/// var hud = new TelemetryHUDExample();
/// hud.Configure(rocketEngine, telemetryStream);
/// 
/// // Each frame:
/// hud.Update();
/// altitudeText.text = hud.AltitudeDisplay;
/// velocityText.text = hud.VelocityDisplay;
/// accelerationText.text = hud.AccelerationDisplay;
/// fuelBar.fillAmount = (float)hud.FuelFraction;
/// apoapsisText.text = hud.ApoapsisDisplay;
/// periapsisText.text = hud.PeriapsisDisplay;
/// </code>
/// </summary>
public class TelemetryHUDExample
{
    private TelemetryStream _stream;
    private RocketSimulationEngine _engine;

    /// <summary>Formatted altitude string.</summary>
    public string AltitudeDisplay { get; private set; } = "0 m";

    /// <summary>Formatted velocity string.</summary>
    public string VelocityDisplay { get; private set; } = "0 m/s";

    /// <summary>Formatted acceleration string.</summary>
    public string AccelerationDisplay { get; private set; } = "0.0 G";

    /// <summary>Fuel fraction (0–1) for UI fill bars.</summary>
    public double FuelFraction { get; private set; } = 1.0;

    /// <summary>Formatted apoapsis altitude string.</summary>
    public string ApoapsisDisplay { get; private set; } = "N/A";

    /// <summary>Formatted periapsis altitude string.</summary>
    public string PeriapsisDisplay { get; private set; } = "N/A";

    /// <summary>Formatted Mach number.</summary>
    public string MachDisplay { get; private set; } = "0.0";

    /// <summary>Formatted dynamic pressure.</summary>
    public string QDisplay { get; private set; } = "0 Pa";

    /// <summary>
    /// Configures the HUD with engine and telemetry stream references.
    /// </summary>
    public void Configure(RocketSimulationEngine engine, TelemetryStream stream)
    {
        _engine = engine;
        _stream = stream;
    }

    /// <summary>
    /// Updates all display strings from current telemetry.
    /// Call once per render frame.
    /// </summary>
    public void Update()
    {
        if (_stream == null) return;

        var t = _stream.Current;

        // Altitude
        if (t.Altitude > 100000)
            AltitudeDisplay = $"{t.Altitude / 1000.0:F1} km";
        else
            AltitudeDisplay = $"{t.Altitude:F0} m";

        // Velocity
        if (t.Speed > 1000)
            VelocityDisplay = $"{t.Speed / 1000.0:F2} km/s";
        else
            VelocityDisplay = $"{t.Speed:F1} m/s";

        // Acceleration in G
        double gForce = t.Acceleration / 9.80665;
        AccelerationDisplay = $"{gForce:F1} G";

        // Fuel
        FuelFraction = t.FuelPercent / 100.0;

        // Orbital elements
        if (t.ApoapsisAltitude > 0 && !double.IsInfinity(t.ApoapsisAltitude))
            ApoapsisDisplay = $"{t.ApoapsisAltitude / 1000.0:F1} km";
        else if (double.IsPositiveInfinity(t.ApoapsisAltitude))
            ApoapsisDisplay = "Escape";
        else
            ApoapsisDisplay = "Sub-orbital";

        if (t.PeriapsisAltitude > 0)
            PeriapsisDisplay = $"{t.PeriapsisAltitude / 1000.0:F1} km";
        else
            PeriapsisDisplay = $"{t.PeriapsisAltitude / 1000.0:F1} km (sub-surface)";

        // Mach (rough estimate: speed / 340 at sea level)
        double soundSpeed = 340.0; // Approximate
        if (t.Altitude < 80000)
        {
            // Rough temperature-based sound speed
            double temp = 288.15 - 0.0065 * Math.Min(t.Altitude, 11000);
            soundSpeed = Math.Sqrt(1.4 * 287.0 * temp);
        }
        double mach = t.Speed / soundSpeed;
        MachDisplay = $"{mach:F1}";

        // Dynamic pressure (approximate)
        double rho = 1.225 * Math.Exp(-t.Altitude / 8500.0);
        double q = 0.5 * rho * t.Speed * t.Speed;
        if (q > 1000)
            QDisplay = $"{q / 1000.0:F1} kPa";
        else
            QDisplay = $"{q:F0} Pa";
    }
}
