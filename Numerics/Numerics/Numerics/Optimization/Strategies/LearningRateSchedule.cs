using System;

namespace CSharpNumerics.Numerics.Optimization.Strategies;

/// <summary>
/// Adjusts the learning rate over the course of training according to a schedule.
/// </summary>
public class LearningRateSchedule
{
    private readonly ScheduleType _type;
    private readonly double _initialLr;
    private readonly double _decayRate;
    private readonly int _decaySteps;
    private readonly double _minLr;

    public enum ScheduleType
    {
        Constant,
        StepDecay,
        ExponentialDecay,
        InverseTimeDecay,
        CosineAnnealing
    }

    public LearningRateSchedule(double initialLr, ScheduleType type = ScheduleType.Constant,
        double decayRate = 0.1, int decaySteps = 100, double minLr = 1e-6)
    {
        _initialLr = initialLr;
        _type = type;
        _decayRate = decayRate;
        _decaySteps = Math.Max(decaySteps, 1);
        _minLr = minLr;
    }

    /// <summary>
    /// Get the learning rate for the given epoch/step.
    /// </summary>
    public double GetLearningRate(int step)
    {
        double lr = _type switch
        {
            ScheduleType.Constant => _initialLr,
            ScheduleType.StepDecay => _initialLr * Math.Pow(_decayRate, step / _decaySteps),
            ScheduleType.ExponentialDecay => _initialLr * Math.Pow(_decayRate, (double)step / _decaySteps),
            ScheduleType.InverseTimeDecay => _initialLr / (1.0 + _decayRate * step),
            ScheduleType.CosineAnnealing => _minLr + 0.5 * (_initialLr - _minLr)
                * (1 + Math.Cos(Math.PI * step / _decaySteps)),
            _ => _initialLr
        };

        return Math.Max(lr, _minLr);
    }
}
