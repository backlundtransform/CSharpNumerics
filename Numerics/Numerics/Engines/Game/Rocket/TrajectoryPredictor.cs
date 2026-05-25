using System;
using System.Collections.Generic;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Mechanics;
using CSharpNumerics.Physics.OrbitalMechanics;

namespace CSharpNumerics.Engines.Game.Rocket;

/// <summary>
/// Fast-forward trajectory propagation for rendering predicted orbit/trajectory lines.
/// Uses simplified 2-body or ballistic propagation without full sim fidelity
/// to produce real-time preview lines (e.g., KSP-style trajectory display).
/// </summary>
public class TrajectoryPredictor
{
    /// <summary>Number of prediction points to generate. Default 200.</summary>
    public int NumPoints { get; set; } = 200;

    /// <summary>Total prediction time horizon (seconds). Default one orbit period or 3600s.</summary>
    public double PredictionHorizon { get; set; } = 3600.0;

    /// <summary>Gravitational parameter (m³/s²). Default Earth.</summary>
    public double Mu { get; set; } = EarthModel.GM;

    /// <summary>Central body radius (m). Default Earth.</summary>
    public double BodyRadius { get; set; } = EarthModel.SemiMajorAxis;

    /// <summary>Whether to include atmospheric drag in prediction. Default false.</summary>
    public bool IncludeDrag { get; set; }

    /// <summary>Ballistic coefficient for drag prediction (kg/m²). Default 100.</summary>
    public double BallisticCoefficient { get; set; } = 100.0;

    /// <summary>Last computed prediction points.</summary>
    public IReadOnlyList<TrajectoryPoint> Points => _points;

    private readonly List<TrajectoryPoint> _points = new List<TrajectoryPoint>();

    /// <summary>
    /// Computes a predicted trajectory from current state using fast Kepler propagation.
    /// </summary>
    /// <param name="position">Current position (m) in inertial frame.</param>
    /// <param name="velocity">Current velocity (m/s) in inertial frame.</param>
    /// <returns>List of trajectory points.</returns>
    public IReadOnlyList<TrajectoryPoint> Predict(Vector position, Vector velocity)
    {
        _points.Clear();

        double r0 = position.GetMagnitude();
        double v0 = velocity.GetMagnitude();

        // Compute orbital elements for period estimation
        double energy = 0.5 * v0 * v0 - Mu / r0;
        double sma = -Mu / (2.0 * energy);

        double horizon = PredictionHorizon;
        if (energy < 0 && sma > 0)
        {
            // Bound orbit: limit prediction to one orbital period
            double period = 2.0 * Math.PI * Math.Sqrt(sma * sma * sma / Mu);
            if (period < horizon) horizon = period;
        }

        double dt = horizon / NumPoints;

        // Velocity-Verlet propagation (fast, no RK4 needed for preview)
        double px = position.x, py = position.y, pz = position.z;
        double vx = velocity.x, vy = velocity.y, vz = velocity.z;

        for (int i = 0; i <= NumPoints; i++)
        {
            double r = Math.Sqrt(px * px + py * py + pz * pz);
            double alt = r - BodyRadius;

            _points.Add(new TrajectoryPoint
            {
                Position = new Vector(px, py, pz),
                Velocity = new Vector(vx, vy, vz),
                Altitude = alt,
                Time = i * dt
            });

            if (alt < 0) break; // Impact

            // Gravity acceleration
            double r3 = r * r * r;
            double ax = -Mu * px / r3;
            double ay = -Mu * py / r3;
            double az = -Mu * pz / r3;

            // Optional drag (exponential atmosphere)
            if (IncludeDrag && alt < 200000 && alt > 0)
            {
                double rho = 1.225 * Math.Exp(-alt / 8500.0);
                double speed = Math.Sqrt(vx * vx + vy * vy + vz * vz);
                if (speed > 0)
                {
                    double dragAccel = 0.5 * rho * speed * speed / BallisticCoefficient;
                    ax -= dragAccel * vx / speed;
                    ay -= dragAccel * vy / speed;
                    az -= dragAccel * vz / speed;
                }
            }

            // Velocity-Verlet step
            double vxHalf = vx + 0.5 * ax * dt;
            double vyHalf = vy + 0.5 * ay * dt;
            double vzHalf = vz + 0.5 * az * dt;

            px += vxHalf * dt;
            py += vyHalf * dt;
            pz += vzHalf * dt;

            // Re-compute acceleration at new position
            double rNew = Math.Sqrt(px * px + py * py + pz * pz);
            double rNew3 = rNew * rNew * rNew;
            double axNew = -Mu * px / rNew3;
            double ayNew = -Mu * py / rNew3;
            double azNew = -Mu * pz / rNew3;

            if (IncludeDrag)
            {
                double altNew = rNew - BodyRadius;
                if (altNew < 200000 && altNew > 0)
                {
                    double rhoNew = 1.225 * Math.Exp(-altNew / 8500.0);
                    double speedNew = Math.Sqrt(vxHalf * vxHalf + vyHalf * vyHalf + vzHalf * vzHalf);
                    if (speedNew > 0)
                    {
                        double dragNew = 0.5 * rhoNew * speedNew * speedNew / BallisticCoefficient;
                        axNew -= dragNew * vxHalf / speedNew;
                        ayNew -= dragNew * vyHalf / speedNew;
                        azNew -= dragNew * vzHalf / speedNew;
                    }
                }
            }

            vx = vxHalf + 0.5 * axNew * dt;
            vy = vyHalf + 0.5 * ayNew * dt;
            vz = vzHalf + 0.5 * azNew * dt;
        }

        return _points;
    }

    /// <summary>
    /// Computes the predicted apoapsis and periapsis altitudes from current state.
    /// </summary>
    public (double apoapsisAltitude, double periapsisAltitude) PredictApsides(Vector position, Vector velocity)
    {
        double r = position.GetMagnitude();
        double v = velocity.GetMagnitude();
        double energy = 0.5 * v * v - Mu / r;

        if (energy >= 0) return (double.PositiveInfinity, r - BodyRadius); // Hyperbolic

        double sma = -Mu / (2.0 * energy);

        // Angular momentum
        Vector h = new Vector(
            position.y * velocity.z - position.z * velocity.y,
            position.z * velocity.x - position.x * velocity.z,
            position.x * velocity.y - position.y * velocity.x);
        double hMag = h.GetMagnitude();
        double ecc = Math.Sqrt(1.0 + 2.0 * energy * hMag * hMag / (Mu * Mu));

        double rApo = sma * (1.0 + ecc);
        double rPeri = sma * (1.0 - ecc);

        return (rApo - BodyRadius, rPeri - BodyRadius);
    }
}

/// <summary>
/// A single point on a predicted trajectory.
/// </summary>
public struct TrajectoryPoint
{
    /// <summary>Position in inertial frame (m).</summary>
    public Vector Position;
    /// <summary>Velocity in inertial frame (m/s).</summary>
    public Vector Velocity;
    /// <summary>Altitude above body surface (m).</summary>
    public double Altitude;
    /// <summary>Time from prediction start (s).</summary>
    public double Time;
}
