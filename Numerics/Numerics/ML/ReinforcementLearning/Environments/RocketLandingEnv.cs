using System;
using System.Collections.Generic;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.ML.ReinforcementLearning.Environments;

/// <summary>
/// RL environment for propulsive rocket landing (Falcon 9 style).
///
/// <b>Observation (8D):</b>
///   [altitude, vx, vy, vz, pitch_angle, pitch_rate, fuel_remaining, speed]
///
/// <b>Action (2D continuous):</b>
///   [throttle (0–1), gimbal_pitch (-1 to 1)]
///
/// <b>Reward:</b>
///   - Soft landing bonus: large reward for touchdown with speed &lt; 2 m/s
///   - Fuel efficiency: small positive reward per step for conserving fuel
///   - Penalty for crash (speed &gt; threshold at ground)
///   - Shaping reward: negative altitude rate (encouraging descent control)
///
/// <b>Done:</b>
///   - Altitude ≤ 0 (landed or crashed), or fuel exhausted, or max steps exceeded
/// </summary>
public class RocketLandingEnv : IEnvironment
{
    public int ObservationSize => 8;
    public int ActionSize => 2;
    public bool IsDiscrete => false;

    private readonly double _dt;
    private readonly int _maxSteps;
    private readonly double _dryMass;
    private readonly double _maxThrust;
    private readonly double _exhaustVelocity;
    private readonly double _gravity;
    private readonly double _maxGimbalAngle;

    private double _x, _y, _z;         // position (z = altitude, up)
    private double _vx, _vy, _vz;      // velocity
    private double _pitch;              // pitch angle from vertical (rad)
    private double _pitchRate;          // pitch rate (rad/s)
    private double _fuelMass;           // remaining fuel (kg)
    private double _initialFuel;
    private int _step;
    private Random _rng;
    private bool _done;

    /// <summary>
    /// Creates a rocket landing environment.
    /// </summary>
    /// <param name="dt">Simulation timestep (s). Default 0.1.</param>
    /// <param name="maxSteps">Max steps per episode. Default 500.</param>
    /// <param name="dryMass">Dry mass (kg). Default 25000.</param>
    /// <param name="fuelMass">Initial fuel mass (kg). Default 5000.</param>
    /// <param name="maxThrust">Max engine thrust (N). Default 800000.</param>
    /// <param name="exhaustVelocity">Effective exhaust velocity (m/s). Default 3000.</param>
    /// <param name="gravity">Gravitational acceleration (m/s²). Default 9.81.</param>
    public RocketLandingEnv(
        double dt = 0.1,
        int maxSteps = 500,
        double dryMass = 25000.0,
        double fuelMass = 5000.0,
        double maxThrust = 800000.0,
        double exhaustVelocity = 3000.0,
        double gravity = 9.81)
    {
        _dt = dt;
        _maxSteps = maxSteps;
        _dryMass = dryMass;
        _initialFuel = fuelMass;
        _maxThrust = maxThrust;
        _exhaustVelocity = exhaustVelocity;
        _gravity = gravity;
        _maxGimbalAngle = 5.0 * Math.PI / 180.0;
        _rng = new Random();
    }

    public (VectorN state, Dictionary<string, object> info) Reset(int? seed = null)
    {
        _rng = seed.HasValue ? new Random(seed.Value) : new Random();

        // Start at ~1000m altitude, descending
        _z = 800.0 + _rng.NextDouble() * 400.0; // 800-1200 m
        _x = (_rng.NextDouble() - 0.5) * 200.0; // ±100 m lateral
        _y = (_rng.NextDouble() - 0.5) * 200.0;
        _vx = (_rng.NextDouble() - 0.5) * 20.0; // ±10 m/s lateral
        _vy = (_rng.NextDouble() - 0.5) * 20.0;
        _vz = -(50.0 + _rng.NextDouble() * 50.0); // -50 to -100 m/s (descending)
        _pitch = (_rng.NextDouble() - 0.5) * 0.2; // small tilt
        _pitchRate = 0;
        _fuelMass = _initialFuel;
        _step = 0;
        _done = false;

        return (GetObservation(), new Dictionary<string, object>());
    }

    public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(int action)
    {
        // Map discrete to continuous: not primary use case
        double throttle = action == 0 ? 0.0 : 1.0;
        return Step(new VectorN(new double[] { throttle, 0.0 }));
    }

    public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(VectorN action)
    {
        if (_done)
            return (GetObservation(), 0, true, new Dictionary<string, object>());

        _step++;
        double throttle = Math.Max(0, Math.Min(1, action[0]));
        double gimbal = Math.Max(-1, Math.Min(1, action[1])) * _maxGimbalAngle;

        // Compute thrust
        double totalMass = _dryMass + _fuelMass;
        double thrust = 0;
        if (_fuelMass > 0)
        {
            thrust = throttle * _maxThrust;
            double massFlow = thrust / _exhaustVelocity;
            double fuelUsed = massFlow * _dt;
            if (fuelUsed > _fuelMass)
            {
                fuelUsed = _fuelMass;
                thrust = fuelUsed * _exhaustVelocity / _dt;
            }
            _fuelMass -= fuelUsed;
        }

        // Thrust direction: primarily vertical with gimbal offset
        double thrustAx = (thrust / totalMass) * Math.Sin(_pitch + gimbal);
        double thrustAz = (thrust / totalMass) * Math.Cos(_pitch + gimbal);

        // Attitude dynamics (simplified 2D pitch)
        double momentArm = 20.0;
        double inertia = totalMass * 100.0; // rough approximation
        double torque = thrust * momentArm * Math.Sin(gimbal);
        _pitchRate += (torque / inertia) * _dt;
        _pitch += _pitchRate * _dt;

        // Translation dynamics
        _vx += thrustAx * _dt;
        _vz += (thrustAz - _gravity) * _dt;
        _x += _vx * _dt;
        _y += _vy * _dt;
        _z += _vz * _dt;

        // Compute reward
        double speed = Math.Sqrt(_vx * _vx + _vy * _vy + _vz * _vz);
        double reward = 0;

        // Shaping: reward controlled descent
        reward += 0.01 * (1.0 - throttle); // fuel efficiency
        reward -= 0.001 * Math.Abs(_pitch); // stay upright

        var info = new Dictionary<string, object>();

        // Check termination
        if (_z <= 0)
        {
            _z = 0;
            _done = true;
            if (speed < 2.0 && Math.Abs(_pitch) < 0.1)
            {
                reward += 100.0; // soft landing
                info["landed"] = true;
            }
            else if (speed < 5.0)
            {
                reward += 20.0; // hard but survivable
                info["landed"] = true;
            }
            else
            {
                reward -= 50.0; // crash
                info["crashed"] = true;
            }
        }
        else if (_step >= _maxSteps)
        {
            _done = true;
            reward -= 20.0; // ran out of time
            info["timeout"] = true;
        }
        else if (_fuelMass <= 0 && _z > 50.0)
        {
            // Out of fuel at high altitude => inevitable crash
            reward -= 10.0;
        }

        return (GetObservation(), reward, _done, info);
    }

    private VectorN GetObservation()
    {
        double speed = Math.Sqrt(_vx * _vx + _vy * _vy + _vz * _vz);
        return new VectorN(new double[]
        {
            _z / 1000.0,          // normalized altitude
            _vx / 100.0,         // normalized horizontal vel
            _vy / 100.0,
            _vz / 100.0,         // normalized vertical vel
            _pitch,               // pitch angle
            _pitchRate,           // pitch rate
            _fuelMass / _initialFuel, // fuel fraction
            speed / 100.0         // normalized speed
        });
    }
}
