using CSharpNumerics.Statistics.Data;
using System;
using System.Linq;

namespace CSharpNumerics.Engines.Exoplanet.Data;

public class LightCurve
{
    public double[] Time { get; }
    public double[] Flux { get; }
    public double[] FluxError { get; }
    public int[] QualityFlags { get; }
    public LightCurveMetadata Metadata { get; }

    public int Length => Time.Length;

    public LightCurve(double[] time, double[] flux, double[] fluxError, int[] qualityFlags, LightCurveMetadata metadata)
    {
        if (time == null) throw new ArgumentNullException(nameof(time));
        if (flux == null) throw new ArgumentNullException(nameof(flux));
        if (fluxError == null) throw new ArgumentNullException(nameof(fluxError));
        if (qualityFlags == null) throw new ArgumentNullException(nameof(qualityFlags));
        if (metadata == null) throw new ArgumentNullException(nameof(metadata));

        if (flux.Length != time.Length)
            throw new ArgumentException("Flux array must match Time length.");
        if (fluxError.Length != time.Length)
            throw new ArgumentException("FluxError array must match Time length.");
        if (qualityFlags.Length != time.Length)
            throw new ArgumentException("QualityFlags array must match Time length.");

        Time = time;
        Flux = flux;
        FluxError = fluxError;
        QualityFlags = qualityFlags;
        Metadata = metadata;
    }

    public static LightCurve FromTimeSeries(TimeSeries ts, string timeCol, string fluxCol, string errCol)
    {
        if (ts == null) throw new ArgumentNullException(nameof(ts));

        int timeIdx = Array.IndexOf(ts.Cols, timeCol);
        int fluxIdx = Array.IndexOf(ts.Cols, fluxCol);
        int errIdx = Array.IndexOf(ts.Cols, errCol);

        if (timeIdx < 0) throw new ArgumentException($"Column '{timeCol}' not found in TimeSeries.");
        if (fluxIdx < 0) throw new ArgumentException($"Column '{fluxCol}' not found in TimeSeries.");
        if (errIdx < 0) throw new ArgumentException($"Column '{errCol}' not found in TimeSeries.");

        double[] time = ts.Data[timeIdx];
        double[] flux = ts.Data[fluxIdx];
        double[] fluxError = ts.Data[errIdx];
        int[] qualityFlags = new int[ts.RowCount];

        return new LightCurve(time, flux, fluxError, qualityFlags, new LightCurveMetadata());
    }

    public static LightCurve FromArrays(double[] time, double[] flux, double[] fluxError)
    {
        int[] qualityFlags = new int[time.Length];
        return new LightCurve(time, flux, fluxError, qualityFlags, new LightCurveMetadata());
    }

    public static LightCurve FromArrays(double[] time, double[] flux)
    {
        double[] fluxError = new double[time.Length];
        int[] qualityFlags = new int[time.Length];
        return new LightCurve(time, flux, fluxError, qualityFlags, new LightCurveMetadata());
    }
}
