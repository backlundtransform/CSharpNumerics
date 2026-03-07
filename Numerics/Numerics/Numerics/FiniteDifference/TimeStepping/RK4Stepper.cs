using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Numerics.FiniteDifference.TimeStepping;

/// <summary>
/// Classical 4th-order Runge-Kutta time stepper (RK4).
/// Fourth-order accurate, O(h⁴). 4 function evaluations per step.
/// The workhorse for most non-stiff ODE/PDE problems with fixed step size.
/// </summary>
public class RK4Stepper : ITimeStepper
{
    public string Name => "RK4";

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

            var k1 = rhs(t, y);
            var k2 = rhs(t + 0.5 * h, y + 0.5 * h * k1);
            var k3 = rhs(t + 0.5 * h, y + 0.5 * h * k2);
            var k4 = rhs(t + h, y + h * k3);
            evals += 4;

            y = y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
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
