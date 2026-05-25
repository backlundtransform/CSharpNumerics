using System;

namespace CSharpNumerics.Engines.Game.Rocket;

/// <summary>
/// Automatic throttle reduction near Max-Q to keep dynamic pressure below structural limits.
/// Computes a throttle multiplier based on current dynamic pressure relative to limits.
/// 
/// Typical use: reduce throttle to 60–80% around Max-Q, return to full throttle after.
/// This is standard practice on vehicles like Falcon 9 and Space Shuttle.
/// </summary>
public class ThrottleBucket
{
    /// <summary>Dynamic pressure threshold to begin throttle reduction (Pa).</summary>
    public double QOnset { get; }

    /// <summary>Target maximum dynamic pressure limit (Pa).</summary>
    public double QTarget { get; }

    /// <summary>Minimum throttle fraction during bucket (e.g. 0.6).</summary>
    public double MinThrottle { get; }

    /// <summary>True if currently in throttle bucket (reducing throttle).</summary>
    public bool IsActive { get; private set; }

    /// <summary>True once Max-Q has passed and throttle has returned to 1.0.</summary>
    public bool IsComplete { get; private set; }

    /// <summary>
    /// Creates a throttle bucket controller.
    /// </summary>
    /// <param name="qOnset">Q at which to begin reducing throttle (Pa). Default 25 kPa.</param>
    /// <param name="qTarget">Target Q limit not to exceed (Pa). Default 35 kPa.</param>
    /// <param name="minThrottle">Minimum throttle during bucket (default 0.6).</param>
    public ThrottleBucket(double qOnset = 25000, double qTarget = 35000, double minThrottle = 0.6)
    {
        if (qOnset <= 0) throw new ArgumentOutOfRangeException(nameof(qOnset));
        if (qTarget <= qOnset) throw new ArgumentOutOfRangeException(nameof(qTarget));
        if (minThrottle <= 0 || minThrottle >= 1) throw new ArgumentOutOfRangeException(nameof(minThrottle));

        QOnset = qOnset;
        QTarget = qTarget;
        MinThrottle = minThrottle;
    }

    /// <summary>
    /// Computes the throttle multiplier for the current dynamic pressure.
    /// Returns 1.0 (full throttle) below onset, reduces proportionally up to QTarget.
    /// After Max-Q has passed (q decreasing below onset), returns 1.0 permanently.
    /// </summary>
    /// <param name="currentQ">Current dynamic pressure in Pascals.</param>
    /// <param name="maxQPassed">True if Max-Q has already been detected.</param>
    /// <returns>Throttle multiplier in range [MinThrottle, 1.0].</returns>
    public double ComputeThrottle(double currentQ, bool maxQPassed)
    {
        if (IsComplete)
            return 1.0;

        if (maxQPassed)
        {
            IsActive = false;
            IsComplete = true;
            return 1.0;
        }

        if (currentQ < QOnset)
        {
            IsActive = false;
            return 1.0;
        }

        // In the bucket: interpolate between 1.0 at onset and MinThrottle at QTarget
        IsActive = true;
        double t = Math.Clamp((currentQ - QOnset) / (QTarget - QOnset), 0, 1);
        return 1.0 - t * (1.0 - MinThrottle);
    }

    /// <summary>
    /// Resets the throttle bucket to initial state.
    /// </summary>
    public void Reset()
    {
        IsActive = false;
        IsComplete = false;
    }
}
