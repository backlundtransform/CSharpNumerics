using System.Collections.Generic;

namespace CSharpNumerics.Engines.Exoplanet.Data;

public class TransitFeatureSet
{
    public Dictionary<string, double> Features { get; }

    public TransitFeatureSet()
    {
        Features = new Dictionary<string, double>();
    }

    public TransitFeatureSet(Dictionary<string, double> features)
    {
        Features = new Dictionary<string, double>(features);
    }

    public double this[string key]
    {
        get => Features[key];
        set => Features[key] = value;
    }

    public bool HasFeature(string key) => Features.ContainsKey(key);

    public static class FeatureNames
    {
        public const string Depth = "Depth";
        public const string Duration = "Duration";
        public const string Period = "Period";
        public const string SnrBls = "SnrBls";
        public const string OddEvenRatio = "OddEvenRatio";
        public const string VShapeMetric = "VShapeMetric";
        public const string SecondaryDepth = "SecondaryDepth";
        public const string ScatterInTransit = "ScatterInTransit";
        public const string ScatterOutTransit = "ScatterOutTransit";
        public const string LimbDarkeningU1 = "LimbDarkeningU1";
        public const string LimbDarkeningU2 = "LimbDarkeningU2";
        public const string IngressEgressRatio = "IngressEgressRatio";

        public static readonly string[] All = new[]
        {
            Depth, Duration, Period, SnrBls, OddEvenRatio, VShapeMetric,
            SecondaryDepth, ScatterInTransit, ScatterOutTransit,
            LimbDarkeningU1, LimbDarkeningU2, IngressEgressRatio
        };
    }
}
