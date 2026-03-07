using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Numerics.FiniteDifference.TimeStepping;

/// <summary>
/// Result of a time-stepping integration.
/// Contains the final state plus optional trajectory, diagnostics, and step statistics.
/// </summary>
public class TimeStepResult
{
    /// <summary>Final time reached.</summary>
    public double T { get; set; }

    /// <summary>Final state vector.</summary>
    public VectorN Y { get; set; }

    /// <summary>
    /// Full trajectory as (time, state) pairs.
    /// Populated when the stepper is configured to record history.
    /// </summary>
    public List<(double t, VectorN y)> Trajectory { get; set; }

    /// <summary>Total number of steps taken (may exceed expected count for adaptive steppers).</summary>
    public int Steps { get; set; }

    /// <summary>Number of rejected steps (adaptive steppers only).</summary>
    public int RejectedSteps { get; set; }

    /// <summary>Number of right-hand side evaluations.</summary>
    public int FunctionEvaluations { get; set; }

    /// <summary>Estimated local error at the final step (adaptive steppers only).</summary>
    public double LastError { get; set; }

    /// <summary>Actual step size used at the final step.</summary>
    public double LastStepSize { get; set; }
}
