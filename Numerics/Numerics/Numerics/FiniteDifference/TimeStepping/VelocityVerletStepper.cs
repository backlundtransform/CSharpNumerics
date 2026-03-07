using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Numerics.FiniteDifference.TimeStepping;

/// <summary>
/// Velocity Verlet time stepper for second-order ODE systems: x'' = a(t, x, v).
/// 
/// State vector y = [positions | velocities] (must have even length).
/// The rhs function must return [velocities | accelerations].
/// Internally uses the Verlet algorithm:
///   x_{n+1} = x_n + v_n·h + 0.5·a_n·h²
///   a_{n+1} = f(t_{n+1}, x_{n+1}, v_n)
///   v_{n+1} = v_n + 0.5·(a_n + a_{n+1})·h
/// 
/// Properties:
/// • Symplectic — preserves phase-space volume, excellent long-term energy conservation
/// • O(h²) accuracy
/// • Time-reversible
/// • Ideal for N-body, molecular dynamics, orbital mechanics, rigid body physics
/// </summary>
public class VelocityVerletStepper : ITimeStepper
{
    public string Name => "Velocity Verlet";

    /// <summary>
    /// Solve a second-order system.
    /// rhs(t, y) must return [velocities | accelerations] given y = [positions | velocities].
    /// </summary>
    public TimeStepResult Solve(
        Func<double, VectorN, VectorN> rhs,
        double t0, double tEnd, VectorN y0, double dt,
        bool recordTrajectory = false,
        Action<double, VectorN> onStep = null)
    {
        int n = y0.Length;
        if (n % 2 != 0)
            throw new ArgumentException("State vector must have even length [positions | velocities].");

        int half = n / 2;
        var y = new VectorN(y0.Values);
        double t = t0;
        int steps = 0;
        int evals = 0;

        List<(double t, VectorN y)> trajectory = null;
        if (recordTrajectory)
        {
            trajectory = new List<(double, VectorN)>();
            trajectory.Add((t, new VectorN(y.Values)));
        }

        // Initial acceleration: extract from rhs
        var dy = rhs(t, y);
        evals++;
        var a = new double[half];
        for (int i = 0; i < half; i++)
            a[i] = dy[half + i]; // acceleration part

        while (t < tEnd)
        {
            double h = Math.Min(dt, tEnd - t);

            // Position update: x += v·h + 0.5·a·h²
            for (int i = 0; i < half; i++)
                y[i] += h * y[half + i] + 0.5 * h * h * a[i];

            // New acceleration at updated position
            var dyNew = rhs(t + h, y);
            evals++;
            var aNew = new double[half];
            for (int i = 0; i < half; i++)
                aNew[i] = dyNew[half + i];

            // Velocity update: v += 0.5·(a + aNew)·h
            for (int i = 0; i < half; i++)
                y[half + i] += 0.5 * h * (a[i] + aNew[i]);

            a = aNew;
            t += h;
            steps++;

            onStep?.Invoke(t, y);
            if (recordTrajectory)
                trajectory.Add((t, new VectorN(y.Values)));
        }

        return new TimeStepResult
        {
            T = t,
            Y = y,
            Trajectory = trajectory,
            Steps = steps,
            FunctionEvaluations = evals,
            LastStepSize = dt
        };
    }
}
