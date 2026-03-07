using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Numerics.FiniteDifference.TimeStepping;

/// <summary>
/// Explicit (forward) Euler time stepper: y_{n+1} = y_n + h · f(t_n, y_n).
/// First-order accurate, O(h). Simple and fast but requires small step sizes for stability.
/// </summary>
public class EulerStepper : ITimeStepper
{
    public string Name => "Euler";

    public TimeStepResult Solve(
        Func<double, VectorN, VectorN> rhs,
        double t0, double tEnd, VectorN y0, double dt,
        bool recordTrajectory = false,
        Action<double, VectorN> onStep = null)
    {
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

        while (t < tEnd)
        {
            double h = Math.Min(dt, tEnd - t);

            var dy = rhs(t, y);
            evals++;

            y = y + h * dy;
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
