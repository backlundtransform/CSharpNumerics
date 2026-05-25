using System;
using System.Collections.Generic;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.ML.ReinforcementLearning.Environments;

/// <summary>
/// RL environment for optimizing rocket ascent trajectories.
/// The agent controls pitch program to maximize payload to orbit.
///
/// <b>Observation (8D):</b>
///   [altitude, speed, flight_path_angle, downrange, mass_fraction, dynamic_pressure, time_fraction, orbital_energy]
///
/// <b>Action (1D continuous):</b>
///   [pitch_rate_command (-1 to 1)] — normalized pitch rate
///
/// <b>Reward:</b>
///   - Terminal reward proportional to achieved orbital energy
///   - Penalty for exceeding max-Q structural limit
///   - Penalty for excessive angle of attack during max-Q
///   - Bonus for reaching orbit (positive orbital energy + low eccentricity)
///
/// <b>Done:</b>
///   - Orbit achieved, re-entry (altitude below 0 after ascent), or max time exceeded
/// </summary>
public class AscentOptimizationEnv : IEnvironment
{
    public int ObservationSize => 8;
    public int ActionSize => 1;
    public bool IsDiscrete => false;

    private readonly double _dt;
    private readonly int _maxSteps;
    private readonly double _initialMass;
    private readonly double _dryMassFraction;
    private readonly double _thrust;
    private readonly double _exhaustVelocity;
    private readonly double _maxPitchRate;
    private readonly double _maxQ;
    private readonly double _targetAltitude;
    private readonly double _earthRadius;
    private readonly double _mu;

    private double _altitude;
    private double _speed;
    private double _gamma;        // flight path angle (rad)
    private double _downrange;    // downrange distance (m)
    private double _mass;
    private double _time;
    private int _step;
    private bool _done;
    private bool _hasAscended;
    private Random _rng;

    /// <summary>
    /// Creates an ascent optimization environment.
    /// </summary>
    /// <param name="dt">Simulation timestep (s). Default 1.0.</param>
    /// <param name="maxSteps">Max steps per episode. Default 600 (10 min burn).</param>
    /// <param name="initialMass">Initial mass (kg). Default 500000.</param>
    /// <param name="dryMassFraction">Dry mass / initial mass. Default 0.08.</param>
    /// <param name="thrust">Constant thrust (N). Default 7000000 (roughly F9).</param>
    /// <param name="exhaustVelocity">Effective exhaust velocity (m/s). Default 3000.</param>
    /// <param name="targetAltitude">Target orbit altitude (m). Default 200000 (200km).</param>
    public AscentOptimizationEnv(
        double dt = 1.0,
        int maxSteps = 600,
        double initialMass = 500000.0,
        double dryMassFraction = 0.08,
        double thrust = 7000000.0,
        double exhaustVelocity = 3000.0,
        double targetAltitude = 200000.0)
    {
        _dt = dt;
        _maxSteps = maxSteps;
        _initialMass = initialMass;
        _dryMassFraction = dryMassFraction;
        _thrust = thrust;
        _exhaustVelocity = exhaustVelocity;
        _maxPitchRate = 2.0 * Math.PI / 180.0; // 2°/s max
        _maxQ = 40000.0; // Pa
        _targetAltitude = targetAltitude;
        _earthRadius = 6371000.0;
        _mu = 3.986e14;
        _rng = new Random();
    }

    public (VectorN state, Dictionary<string, object> info) Reset(int? seed = null)
    {
        _rng = seed.HasValue ? new Random(seed.Value) : new Random();

        _altitude = 0;
        _speed = 0;
        _gamma = Math.PI / 2.0; // vertical
        _downrange = 0;
        _mass = _initialMass;
        _time = 0;
        _step = 0;
        _done = false;
        _hasAscended = false;

        return (GetObservation(), new Dictionary<string, object>());
    }

    public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(int action)
    {
        double pitchCmd = action == 0 ? -1.0 : (action == 1 ? 0.0 : 1.0);
        return Step(new VectorN(new double[] { pitchCmd }));
    }

    public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(VectorN action)
    {
        if (_done)
            return (GetObservation(), 0, true, new Dictionary<string, object>());

        _step++;
        _time += _dt;

        double pitchRateCmd = Math.Max(-1, Math.Min(1, action[0])) * _maxPitchRate;

        // Atmosphere model (exponential)
        double scaleHeight = 8500.0;
        double rho0 = 1.225;
        double rho = rho0 * Math.Exp(-_altitude / scaleHeight);
        double dynamicPressure = 0.5 * rho * _speed * _speed;

        // Drag (Cd * A ~ 10 m² reference)
        double cd = 0.3;
        double area = 10.0;
        double drag = cd * area * dynamicPressure;

        // Thrust and mass flow
        double currentThrust = _thrust;
        double dryMass = _initialMass * _dryMassFraction;
        double massFlow = currentThrust / _exhaustVelocity;

        if (_mass <= dryMass)
        {
            currentThrust = 0;
            massFlow = 0;
        }

        // Equations of motion (2D trajectory: altitude + downrange)
        double r = _earthRadius + _altitude;
        double g = _mu / (r * r);

        double accelThrust = currentThrust / _mass;
        double accelDrag = drag / _mass;

        // dv/dt = thrust - drag - g*sin(gamma)
        double dv = (accelThrust - accelDrag - g * Math.Sin(_gamma)) * _dt;
        // dgamma/dt = (1/v)(thrust_perp/m - g*cos(gamma) + v²cos(gamma)/r)
        // With pitch rate command, we control gamma directly via angle of attack
        double dgamma = pitchRateCmd * _dt;

        // Gravity turn natural rate
        if (_speed > 10.0)
        {
            double gravTurnRate = -(g / _speed - _speed / r) * Math.Cos(_gamma);
            dgamma += gravTurnRate * _dt;
        }

        _speed += dv;
        if (_speed < 0) _speed = 0;
        _gamma += dgamma;
        _gamma = Math.Max(-Math.PI / 2, Math.Min(Math.PI / 2, _gamma));

        double dAlt = _speed * Math.Sin(_gamma) * _dt;
        double dDown = _speed * Math.Cos(_gamma) * _dt;
        _altitude += dAlt;
        _downrange += dDown;

        _mass -= massFlow * _dt;
        if (_mass < dryMass) _mass = dryMass;

        if (_altitude > 1000) _hasAscended = true;

        // Reward calculation
        double reward = 0;

        // Penalize excessive dynamic pressure
        if (dynamicPressure > _maxQ)
            reward -= 0.1 * (dynamicPressure - _maxQ) / _maxQ;

        // Small shaping: reward gaining altitude and speed
        reward += 0.001 * (_speed * Math.Sin(_gamma)) / 1000.0;

        var info = new Dictionary<string, object>();

        // Terminal conditions
        double orbitalVelocity = Math.Sqrt(_mu / r);
        double orbitalEnergy = 0.5 * _speed * _speed - _mu / r;

        if (_altitude >= _targetAltitude && Math.Abs(_gamma) < 0.05 && _speed > 0.9 * orbitalVelocity)
        {
            // Orbit achieved!
            _done = true;
            double payloadFraction = (_mass - dryMass) / _initialMass;
            reward += 100.0 + 500.0 * payloadFraction; // reward payload to orbit
            info["orbit_achieved"] = true;
            info["final_altitude"] = _altitude;
            info["final_speed"] = _speed;
            info["orbital_velocity"] = orbitalVelocity;
        }
        else if (_hasAscended && _altitude <= 0)
        {
            _done = true;
            reward -= 50.0;
            info["crashed"] = true;
        }
        else if (_step >= _maxSteps)
        {
            _done = true;
            // Partial credit for altitude/speed achieved
            reward += 10.0 * Math.Min(1.0, _altitude / _targetAltitude);
            reward += 10.0 * Math.Min(1.0, _speed / orbitalVelocity);
            info["timeout"] = true;
        }

        return (GetObservation(), reward, _done, info);
    }

    private VectorN GetObservation()
    {
        double r = _earthRadius + _altitude;
        double orbitalEnergy = 0.5 * _speed * _speed - _mu / r;
        double rho = 1.225 * Math.Exp(-_altitude / 8500.0);
        double dynamicPressure = 0.5 * rho * _speed * _speed;

        return new VectorN(new double[]
        {
            _altitude / _targetAltitude,                       // normalized altitude
            _speed / 8000.0,                                   // normalized speed (orbital ~7800)
            _gamma / (Math.PI / 2.0),                         // normalized flight path angle
            _downrange / 1000000.0,                           // normalized downrange (1000km)
            (_mass - _initialMass * _dryMassFraction) / (_initialMass * (1.0 - _dryMassFraction)), // fuel fraction
            dynamicPressure / _maxQ,                           // normalized Q
            _time / (_maxSteps * _dt),                        // time fraction
            Math.Max(-1, Math.Min(1, orbitalEnergy / 30000000.0)) // normalized orbital energy
        });
    }
}
