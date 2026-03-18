using System;

namespace CSharpNumerics.Numerics.Optimization.Strategies;

/// <summary>
/// Monitors a metric (typically validation loss) and signals when training
/// should stop because the metric has stopped improving.
/// </summary>
public class EarlyStopping
{
    public int Patience { get; set; }
    public double MinDelta { get; set; }

    private double _bestValue;
    private int _counter;
    private bool _triggered;

    public EarlyStopping(int patience = 10, double minDelta = 1e-4)
    {
        Patience = patience;
        MinDelta = minDelta;
        Reset();
    }

    /// <summary>
    /// Report the current metric value.
    /// Returns true when patience is exhausted (training should stop).
    /// </summary>
    /// <param name="metricValue">Current value of the monitored metric (lower is better).</param>
    public bool Check(double metricValue)
    {
        if (_triggered) return true;

        if (metricValue < _bestValue - MinDelta)
        {
            _bestValue = metricValue;
            _counter = 0;
            return false;
        }

        _counter++;
        if (_counter >= Patience)
        {
            _triggered = true;
            return true;
        }

        return false;
    }

    /// <summary>True after <see cref="Check"/> returned true.</summary>
    public bool ShouldStop => _triggered;

    /// <summary>The best metric value observed so far.</summary>
    public double BestValue => _bestValue;

    public void Reset()
    {
        _bestValue = double.MaxValue;
        _counter = 0;
        _triggered = false;
    }
}
