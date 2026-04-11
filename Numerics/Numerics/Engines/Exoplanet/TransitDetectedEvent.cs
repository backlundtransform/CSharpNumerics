using CSharpNumerics.Engines.Exoplanet.Data;

namespace CSharpNumerics.Engines.Exoplanet;

/// <summary>
/// Event published via <see cref="Common.EventBus"/> when the exoplanet engine detects a transit candidate.
/// </summary>
public class TransitDetectedEvent
{
    /// <summary>The detected transit candidate with fitted parameters and features.</summary>
    public TransitCandidate Candidate { get; }

    /// <summary>Simulation time (seconds) at which the detection occurred.</summary>
    public double Timestamp { get; }

    public TransitDetectedEvent(TransitCandidate candidate, double timestamp)
    {
        Candidate = candidate;
        Timestamp = timestamp;
    }
}
