using System;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Mechanics;

namespace CSharpNumerics.Physics.OrbitalMechanics;

/// <summary>
/// Numerically propagates an orbit in the ECI frame using velocity-Verlet integration.
/// Supports J2 perturbation and optional atmospheric drag for low perigee orbits.
/// </summary>
public class OrbitalPropagator
{
    private Vector _position;
    private Vector _velocity;
    private double _time;

    /// <summary>Gravitational parameter μ (m³/s²). Default: Earth.</summary>
    public double Mu { get; set; } = EarthModel.GM;

    /// <summary>Whether to include J2 perturbation. Default true.</summary>
    public bool IncludeJ2 { get; set; } = true;

    /// <summary>
    /// Whether to include atmospheric drag for low orbits (perigee &lt; DragAltitudeThreshold).
    /// Default false.
    /// </summary>
    public bool IncludeAtmosphericDrag { get; set; }

    /// <summary>Altitude threshold (m) below which atmospheric drag is applied. Default 400 km.</summary>
    public double DragAltitudeThreshold { get; set; } = 400000.0;

    /// <summary>Ballistic coefficient for drag: Cd·A/m (m²/kg). Default 0.01.</summary>
    public double BallisticCoefficient { get; set; } = 0.01;

    /// <summary>Current position in ECI (meters).</summary>
    public Vector Position => _position;

    /// <summary>Current velocity in ECI (m/s).</summary>
    public Vector Velocity => _velocity;

    /// <summary>Current propagation time (seconds).</summary>
    public double Time => _time;

    /// <summary>Current orbital radius (meters).</summary>
    public double Radius => _position.GetMagnitude();

    /// <summary>Current speed (m/s).</summary>
    public double Speed => _velocity.GetMagnitude();

    /// <summary>Current specific orbital energy (J/kg).</summary>
    public double SpecificEnergy =>
        0.5 * Speed * Speed - Mu / Radius;

    /// <summary>Current orbital elements computed from state.</summary>
    public OrbitalElements Elements =>
        StateToElements.FromStateVector(_position, _velocity, Mu);

    /// <summary>
    /// Initializes the propagator with an ECI state vector.
    /// </summary>
    public void SetState(Vector position, Vector velocity)
    {
        _position = position;
        _velocity = velocity;
        _time = 0;
    }

    /// <summary>
    /// Initializes the propagator from orbital elements.
    /// </summary>
    public void SetState(OrbitalElements elements)
    {
        var (pos, vel) = ElementsToState.ToStateVector(elements);
        _position = pos;
        _velocity = vel;
        _time = 0;
        Mu = elements.Mu;
    }

    /// <summary>
    /// Advances the orbit by dt seconds using velocity-Verlet (symplectic) integration.
    /// </summary>
    public void Step(double dt)
    {
        // Velocity Verlet: symplectic, second-order, excellent energy conservation
        Vector acc = ComputeAcceleration(_position, _velocity);

        // Half-step velocity
        _velocity = new Vector(
            _velocity.x + 0.5 * acc.x * dt,
            _velocity.y + 0.5 * acc.y * dt,
            _velocity.z + 0.5 * acc.z * dt);

        // Full-step position
        _position = new Vector(
            _position.x + _velocity.x * dt,
            _position.y + _velocity.y * dt,
            _position.z + _velocity.z * dt);

        // New acceleration at updated position
        Vector accNew = ComputeAcceleration(_position, _velocity);

        // Half-step velocity
        _velocity = new Vector(
            _velocity.x + 0.5 * accNew.x * dt,
            _velocity.y + 0.5 * accNew.y * dt,
            _velocity.z + 0.5 * accNew.z * dt);

        _time += dt;
    }

    /// <summary>
    /// Propagates for a given total duration with the specified time step.
    /// </summary>
    public void Propagate(double duration, double dt)
    {
        int steps = (int)(duration / dt);
        for (int i = 0; i < steps; i++)
            Step(dt);

        double remaining = duration - steps * dt;
        if (remaining > 1e-10)
            Step(remaining);
    }

    private Vector ComputeAcceleration(Vector pos, Vector vel)
    {
        Vector acc;

        if (IncludeJ2)
        {
            acc = GravityModel.Acceleration(pos);
        }
        else
        {
            acc = GravityModel.PointMassAcceleration(pos);
        }

        // Atmospheric drag for low-altitude passes
        if (IncludeAtmosphericDrag)
        {
            double r = pos.GetMagnitude();
            double alt = r - EarthModel.SemiMajorAxis;
            if (alt < DragAltitudeThreshold && alt > 0)
            {
                double density = GetApproximateDensity(alt);
                double speed = vel.GetMagnitude();
                if (speed > 0.1)
                {
                    // a_drag = -0.5 * ρ * v² * (Cd·A/m) * v̂
                    double dragAccelMag = 0.5 * density * speed * speed * BallisticCoefficient;
                    Vector velUnit = (1.0 / speed) * vel;
                    acc = new Vector(
                        acc.x - dragAccelMag * velUnit.x,
                        acc.y - dragAccelMag * velUnit.y,
                        acc.z - dragAccelMag * velUnit.z);
                }
            }
        }

        return acc;
    }

    /// <summary>
    /// Exponential atmosphere approximation for orbital drag.
    /// ρ(h) = ρ₀ · exp(-h/H) where H ≈ 8500 m scale height.
    /// </summary>
    private static double GetApproximateDensity(double altitude)
    {
        if (altitude < 0) return 1.225;
        if (altitude > 600000) return 0;

        // Simple exponential model: adequate for drag estimation
        const double rho0 = 1.225; // sea level density
        const double H = 8500.0;   // scale height
        return rho0 * Math.Exp(-altitude / H);
    }
}
