using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace CSharpNumerics.Engines.Game.Performance;

/// <summary>
/// Lightweight profiling timer for measuring game engine performance.
///
/// Tracks named timing sections with running statistics (min, max, average, count).
/// Designed for per-frame profiling of physics, fluid, collision, and rendering subsystems.
///
/// Usage:
///   profiler.Begin("Physics");
///   world.Step(dt);
///   profiler.End("Physics");
///   
///   var stats = profiler.GetStats("Physics");
///   Console.WriteLine($"Physics: avg={stats.AverageMs:F2}ms, max={stats.MaxMs:F2}ms");
/// </summary>
public class ProfileTimer
{
    /// <summary>
    /// Timing statistics for a named section.
    /// </summary>
    public struct TimingStats
    {
        /// <summary>Section name.</summary>
        public string Name;

        /// <summary>Number of measurements.</summary>
        public int Count;

        /// <summary>Total time in milliseconds.</summary>
        public double TotalMs;

        /// <summary>Minimum time in milliseconds.</summary>
        public double MinMs;

        /// <summary>Maximum time in milliseconds.</summary>
        public double MaxMs;

        /// <summary>Average time in milliseconds.</summary>
        public double AverageMs => Count > 0 ? TotalMs / Count : 0;

        /// <summary>Last measured time in milliseconds.</summary>
        public double LastMs;

        public override string ToString() =>
            $"{Name}: avg={AverageMs:F3}ms, min={MinMs:F3}ms, max={MaxMs:F3}ms, count={Count}";
    }

    private struct ActiveTimer
    {
        public long StartTicks;
    }

    private readonly Dictionary<string, TimingStats> _stats = new();
    private readonly Dictionary<string, ActiveTimer> _active = new();
    private static readonly double TicksToMs = 1000.0 / Stopwatch.Frequency;

    /// <summary>
    /// Begin timing a named section.
    /// </summary>
    public void Begin(string name)
    {
        _active[name] = new ActiveTimer { StartTicks = Stopwatch.GetTimestamp() };
    }

    /// <summary>
    /// End timing a named section and record the measurement.
    /// </summary>
    public void End(string name)
    {
        long endTicks = Stopwatch.GetTimestamp();

        if (!_active.TryGetValue(name, out var timer))
            return;

        double ms = (endTicks - timer.StartTicks) * TicksToMs;
        _active.Remove(name);

        if (_stats.TryGetValue(name, out var stats))
        {
            stats.Count++;
            stats.TotalMs += ms;
            stats.LastMs = ms;
            if (ms < stats.MinMs) stats.MinMs = ms;
            if (ms > stats.MaxMs) stats.MaxMs = ms;
            _stats[name] = stats;
        }
        else
        {
            _stats[name] = new TimingStats
            {
                Name = name,
                Count = 1,
                TotalMs = ms,
                MinMs = ms,
                MaxMs = ms,
                LastMs = ms,
            };
        }
    }

    /// <summary>
    /// Measure a section using an IDisposable scope.
    /// Usage: using (profiler.Scope("Physics")) { ... }
    /// </summary>
    public ProfileScope Scope(string name)
    {
        Begin(name);
        return new ProfileScope(this, name);
    }

    /// <summary>
    /// Get timing statistics for a named section.
    /// </summary>
    public TimingStats GetStats(string name)
    {
        return _stats.TryGetValue(name, out var stats) ? stats : default;
    }

    /// <summary>
    /// Get all timing statistics.
    /// </summary>
    public IReadOnlyDictionary<string, TimingStats> AllStats => _stats;

    /// <summary>
    /// Reset all statistics.
    /// </summary>
    public void Reset()
    {
        _stats.Clear();
        _active.Clear();
    }

    /// <summary>
    /// Format a summary of all measured sections.
    /// </summary>
    public string Summary()
    {
        var sb = new System.Text.StringBuilder();
        sb.AppendLine("=== Profile Summary ===");
        foreach (var kvp in _stats)
            sb.AppendLine(kvp.Value.ToString());
        return sb.ToString();
    }

    /// <summary>
    /// Disposable scope for automatic Begin/End pairing.
    /// </summary>
    public struct ProfileScope : IDisposable
    {
        private readonly ProfileTimer _timer;
        private readonly string _name;

        internal ProfileScope(ProfileTimer timer, string name)
        {
            _timer = timer;
            _name = name;
        }

        public void Dispose()
        {
            _timer.End(_name);
        }
    }
}
