using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Numerics.FiniteDifference.TimeStepping;

/// <summary>
/// Dormand–Prince adaptive Runge-Kutta 4(5) time stepper.
/// Embeds a 4th-order and 5th-order solution in 6 stages (FSAL — first same as last)
/// to estimate local truncation error and automatically adjust the step size.
/// 
/// This is the method behind MATLAB's ode45 and SciPy's RK45.
/// 
/// Properties:
/// • Automatic step-size control — no manual dt tuning needed
/// • 4th-order local accuracy with 5th-order error estimator
/// • FSAL: only 6 new evaluations per accepted step (k7 reused as k1 of next step)
/// • Safety factor and min/max step ratios prevent runaway growth or collapse
/// </summary>
public class AdaptiveRK45Stepper : ITimeStepper
{
    /// <summary>Absolute error tolerance (default 1e-6).</summary>
    public double AbsoluteTolerance { get; set; } = 1e-6;

    /// <summary>Relative error tolerance (default 1e-6).</summary>
    public double RelativeTolerance { get; set; } = 1e-6;

    /// <summary>Safety factor for step size adjustment (default 0.9).</summary>
    public double SafetyFactor { get; set; } = 0.9;

    /// <summary>Maximum step growth factor per step (default 5.0).</summary>
    public double MaxGrowth { get; set; } = 5.0;

    /// <summary>Minimum step shrink factor per step (default 0.2).</summary>
    public double MinShrink { get; set; } = 0.2;

    /// <summary>Maximum number of steps before aborting (default 1_000_000).</summary>
    public int MaxSteps { get; set; } = 1_000_000;

    public string Name => "Dormand-Prince 4(5)";

    // ═══════════════════════════════════════════════════════════════
    //  Dormand–Prince coefficients
    // ═══════════════════════════════════════════════════════════════

    // Nodes (c)
    private static readonly double c2 = 1.0 / 5.0;
    private static readonly double c3 = 3.0 / 10.0;
    private static readonly double c4 = 4.0 / 5.0;
    private static readonly double c5 = 8.0 / 9.0;
    // c6 = 1, c7 = 1

    // A-matrix rows
    private static readonly double a21 = 1.0 / 5.0;
    private static readonly double a31 = 3.0 / 40.0, a32 = 9.0 / 40.0;
    private static readonly double a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0;
    private static readonly double a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0, a54 = -212.0 / 729.0;
    private static readonly double a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0;

    // 5th-order weights (b)
    private static readonly double b1 = 35.0 / 384.0, b3 = 500.0 / 1113.0, b4 = 125.0 / 192.0, b5 = -2187.0 / 6784.0, b6 = 11.0 / 84.0;
    // b2 = 0, b7 = 0

    // Error coefficients: e_i = b_i - b*_i (difference between 5th and 4th order)
    private static readonly double e1 = 71.0 / 57600.0, e3 = -71.0 / 16695.0, e4 = 71.0 / 1920.0, e5 = -17253.0 / 339200.0, e6 = 22.0 / 525.0, e7 = -1.0 / 40.0;

    // ═══════════════════════════════════════════════════════════════
    //  Solve
    // ═══════════════════════════════════════════════════════════════

    public TimeStepResult Solve(
        Func<double, VectorN, VectorN> rhs,
        double t0, double tEnd, VectorN y0, double dt,
        bool recordTrajectory = false,
        Action<double, VectorN> onStep = null)
    {
        var y = new VectorN(y0.Values);
        double t = t0;
        double h = dt;
        int steps = 0;
        int rejected = 0;
        int evals = 0;
        double lastError = 0;

        List<(double t, VectorN y)> trajectory = null;
        if (recordTrajectory)
        {
            trajectory = new List<(double, VectorN)>();
            trajectory.Add((t, new VectorN(y.Values)));
        }

        // First evaluation (FSAL: reuse k7 from previous step as k1)
        var k1 = rhs(t, y);
        evals++;

        while (t < tEnd)
        {
            if (steps >= MaxSteps)
                throw new InvalidOperationException(
                    $"AdaptiveRK45Stepper exceeded {MaxSteps} steps. " +
                    "Consider relaxing tolerances or increasing MaxSteps.");

            // Don't overshoot tEnd
            h = Math.Min(h, tEnd - t);

            // ── Compute stages ───────────────────────────────────
            var k2 = rhs(t + c2 * h, y + h * (a21 * k1));
            var k3 = rhs(t + c3 * h, y + h * (a31 * k1 + a32 * k2));
            var k4 = rhs(t + c4 * h, y + h * (a41 * k1 + a42 * k2 + a43 * k3));
            var k5 = rhs(t + c5 * h, y + h * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4));
            var k6 = rhs(t + h, y + h * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5));
            evals += 5;

            // ── 5th-order solution ───────────────────────────────
            var yNew = y + h * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6);

            // ── Error estimate (difference between 5th and 4th order) ──
            var errVec = h * (e1 * k1 + e3 * k3 + e4 * k4 + e5 * k5 + e6 * k6);

            // Need k7 for FSAL and error
            var k7 = rhs(t + h, yNew);
            evals++;
            errVec = errVec + h * e7 * k7;

            // ── Error norm (mixed absolute/relative) ─────────────
            double err = ErrorNorm(y, yNew, errVec);
            lastError = err;

            if (err <= 1.0)
            {
                // Accept step
                t += h;
                y = yNew;
                k1 = k7; // FSAL
                steps++;

                onStep?.Invoke(t, y);
                if (recordTrajectory)
                    trajectory.Add((t, new VectorN(y.Values)));

                // Grow step size
                double factor = err > 0
                    ? SafetyFactor * Math.Pow(err, -0.2)
                    : MaxGrowth;
                factor = Math.Min(factor, MaxGrowth);
                h *= factor;
            }
            else
            {
                // Reject step — shrink and retry
                rejected++;
                double factor = SafetyFactor * Math.Pow(err, -0.25);
                factor = Math.Max(factor, MinShrink);
                h *= factor;
            }
        }

        return new TimeStepResult
        {
            T = t,
            Y = y,
            Trajectory = trajectory,
            Steps = steps,
            RejectedSteps = rejected,
            FunctionEvaluations = evals,
            LastError = lastError,
            LastStepSize = h
        };
    }

    // ═══════════════════════════════════════════════════════════════
    //  Error norm: RMS of err_i / (atol + rtol * max(|y_i|, |yNew_i|))
    // ═══════════════════════════════════════════════════════════════

    private double ErrorNorm(VectorN y, VectorN yNew, VectorN err)
    {
        double sum = 0;
        int n = y.Length;

        for (int i = 0; i < n; i++)
        {
            double scale = AbsoluteTolerance + RelativeTolerance * Math.Max(Math.Abs(y[i]), Math.Abs(yNew[i]));
            double ratio = err[i] / scale;
            sum += ratio * ratio;
        }

        return Math.Sqrt(sum / n);
    }
}
