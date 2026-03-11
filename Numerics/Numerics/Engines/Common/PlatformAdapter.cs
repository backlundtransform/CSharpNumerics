using System;

namespace CSharpNumerics.Engines.Common
{
    /// <summary>
    /// Describes the runtime platform.
    /// </summary>
    public enum RuntimePlatform
    {
        Standalone,
        Unity,
        Web,
        Other
    }

    /// <summary>
    /// Abstraction layer for platform-specific behaviour.
    /// Override in host applications to provide platform-specific logging, timing, etc.
    /// </summary>
    public class PlatformAdapter
    {
        /// <summary>Current runtime platform.</summary>
        public virtual RuntimePlatform Platform => RuntimePlatform.Standalone;

        /// <summary>Log an informational message.</summary>
        public virtual void Log(string message) => Console.WriteLine(message);

        /// <summary>Log a warning.</summary>
        public virtual void LogWarning(string message) => Console.WriteLine($"[WARN] {message}");

        /// <summary>Log an error.</summary>
        public virtual void LogError(string message) => Console.Error.WriteLine($"[ERROR] {message}");

        /// <summary>High-resolution timestamp in seconds (monotonic clock).</summary>
        public virtual double GetTimestamp() => Environment.TickCount / 1000.0;

        /// <summary>Default singleton for standalone usage.</summary>
        public static PlatformAdapter Default { get; set; } = new PlatformAdapter();
    }
}
