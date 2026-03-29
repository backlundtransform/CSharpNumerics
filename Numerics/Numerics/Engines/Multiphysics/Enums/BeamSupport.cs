namespace CSharpNumerics.Engines.Multiphysics.Enums;

/// <summary>
/// Beam support conditions for beam stress analysis.
/// </summary>
public enum BeamSupport
{
    /// <summary>Fixed at one end, free at the other.</summary>
    Cantilever,

    /// <summary>Supported at both ends with free rotation.</summary>
    SimplySupported,

    /// <summary>Fixed at both ends (no rotation or displacement).</summary>
    FixedFixed
}
