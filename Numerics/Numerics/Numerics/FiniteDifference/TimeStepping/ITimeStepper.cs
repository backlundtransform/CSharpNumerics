using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Numerics.FiniteDifference.TimeStepping;

/// <summary>
/// Contract for time-stepping integrators that advance an ODE system y' = f(t, y) from t0 to tEnd.
/// Complements the existing extension-method solvers in <see cref="DifferentialEquationExtensions"/>
/// with stateful, configurable integration including trajectory recording and per-step callbacks.
/// </summary>
public interface ITimeStepper
{
    /// <summary>Display name of the integration method (e.g. "RK4", "Dormand-Prince 4/5").</summary>
    string Name { get; }

    /// <summary>
    /// Solve y' = rhs(t, y) from t0 to tEnd.
    /// </summary>
    /// <param name="rhs">Right-hand side function f(t, y).</param>
    /// <param name="t0">Initial time.</param>
    /// <param name="tEnd">Final time.</param>
    /// <param name="y0">Initial state vector.</param>
    /// <param name="dt">Initial (or fixed) step size.</param>
    /// <param name="recordTrajectory">If true, the result contains the full trajectory.</param>
    /// <param name="onStep">Optional callback invoked after each accepted step.</param>
    /// <returns>A <see cref="TimeStepResult"/> with final state, diagnostics, and optional trajectory.</returns>
    TimeStepResult Solve(
        Func<double, VectorN, VectorN> rhs,
        double t0,
        double tEnd,
        VectorN y0,
        double dt,
        bool recordTrajectory = false,
        Action<double, VectorN> onStep = null);
}
