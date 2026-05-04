using System;

namespace CSharpNumerics.Engines.GIS.Spread.VolumetricContamination;

/// <summary>
/// Concentration-vs-depth profile at a single horizontal (ix, iy) column.
/// </summary>
public class DepthProfile
{
    /// <summary>Depth values (metres from ZMin, ascending).</summary>
    public double[] Depths { get; }

    /// <summary>Concentration (mg/L) at each depth.</summary>
    public double[] Concentrations { get; }

    /// <summary>Depth at which concentration is highest.</summary>
    public double MaxConcentrationDepth
    {
        get
        {
            if (Concentrations.Length == 0) return 0;
            int best = 0;
            for (int i = 1; i < Concentrations.Length; i++)
                if (Concentrations[i] > Concentrations[best]) best = i;
            return Depths[best];
        }
    }

    /// <summary>Peak concentration value across all depths.</summary>
    public double MaxConcentration
    {
        get
        {
            double max = 0;
            for (int i = 0; i < Concentrations.Length; i++)
                if (Concentrations[i] > max) max = Concentrations[i];
            return max;
        }
    }

    /// <summary>Creates a depth profile.</summary>
    public DepthProfile(double[] depths, double[] concentrations)
    {
        if (depths == null) throw new ArgumentNullException(nameof(depths));
        if (concentrations == null) throw new ArgumentNullException(nameof(concentrations));
        if (depths.Length != concentrations.Length)
            throw new ArgumentException("Depths and concentrations must have the same length.");
        Depths = depths;
        Concentrations = concentrations;
    }
}
