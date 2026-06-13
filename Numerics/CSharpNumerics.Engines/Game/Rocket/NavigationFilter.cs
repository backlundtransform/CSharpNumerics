using System;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Engines.Game.Rocket;

/// <summary>
/// Simple navigation state estimator for game use.
/// Provides either perfect state knowledge or adds configurable noise
/// to simulate real-world sensor uncertainty.
/// </summary>
public class NavigationFilter
{
    /// <summary>Whether to add measurement noise. Default false (perfect knowledge).</summary>
    public bool AddNoise { get; set; }

    /// <summary>Position noise standard deviation (meters). Default 10 m.</summary>
    public double PositionNoiseSigma { get; set; } = 10.0;

    /// <summary>Velocity noise standard deviation (m/s). Default 0.5 m/s.</summary>
    public double VelocityNoiseSigma { get; set; } = 0.5;

    /// <summary>Attitude noise standard deviation (radians). Default 0.1°.</summary>
    public double AttitudeNoiseSigma { get; set; } = 0.1 * Math.PI / 180.0;

    private Random _rng;

    public NavigationFilter(int? seed = null)
    {
        _rng = seed.HasValue ? new Random(seed.Value) : new Random();
    }

    /// <summary>
    /// Returns the estimated (possibly noisy) state given the true state.
    /// </summary>
    /// <param name="trueState">Actual rocket state.</param>
    /// <returns>Estimated rocket state (with noise if enabled).</returns>
    public RocketState Estimate(RocketState trueState)
    {
        if (!AddNoise) return trueState.Clone();

        Vector noisyPos = new Vector(
            trueState.Position.x + GaussianNoise(PositionNoiseSigma),
            trueState.Position.y + GaussianNoise(PositionNoiseSigma),
            trueState.Position.z + GaussianNoise(PositionNoiseSigma));

        Vector noisyVel = new Vector(
            trueState.Velocity.x + GaussianNoise(VelocityNoiseSigma),
            trueState.Velocity.y + GaussianNoise(VelocityNoiseSigma),
            trueState.Velocity.z + GaussianNoise(VelocityNoiseSigma));

        // Add small attitude perturbation
        double noiseAngle = GaussianNoise(AttitudeNoiseSigma);
        Vector noiseAxis = new Vector(
            GaussianNoise(1), GaussianNoise(1), GaussianNoise(1));
        double axisMag = noiseAxis.GetMagnitude();
        Quaternion noisyAtt = trueState.Attitude;
        if (axisMag > 1e-10 && Math.Abs(noiseAngle) > 1e-15)
        {
            noiseAxis = (1.0 / axisMag) * noiseAxis;
            Quaternion perturbation = Quaternion.FromAxisAngle(noiseAxis, noiseAngle);
            noisyAtt = perturbation * trueState.Attitude;
        }

        return new RocketState(noisyPos, noisyVel, noisyAtt, trueState.AngularRate, trueState.Mass);
    }

    private double GaussianNoise(double sigma)
    {
        // Box-Muller transform
        double u1 = 1.0 - _rng.NextDouble(); // avoid log(0)
        double u2 = _rng.NextDouble();
        return sigma * Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }
}
