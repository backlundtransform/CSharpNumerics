using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Engines.Multiphysics.Snapshots;
using System;
using System.Globalization;
using System.IO;
using System.Text;

namespace CSharpNumerics.Engines.Multiphysics.Export;

/// <summary>
/// Exports multiphysics simulation data as JSON for web visualization.
/// Zero external dependencies — builds JSON strings with <see cref="StringBuilder"/>.
/// <para>
/// Output structure:
/// <code>
/// {
///   "type": "HeatPlate",
///   "dimensions": { "nx": 20, "ny": 20, "dx": 0.005, "dy": 0.005 },
///   "metadata": { ... },
///   "timeSteps": [
///     { "time": 0.0, "stepIndex": 0, "min": ..., "max": ..., "values": [...] },
///     ...
///   ]
/// }
/// </code>
/// </para>
/// </summary>
public static class MultiphysicsJsonExporter
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    // ═══════════════════════════════════════════════════════════════
    //  2D Timeline
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Serialize a <see cref="SimulationTimeline"/> to JSON.
    /// </summary>
    public static string ToJson(SimulationTimeline timeline, ExportMetadata metadata = null)
    {
        if (timeline.Count == 0)
            throw new ArgumentException("Timeline is empty.");

        var first = timeline[0];
        var sb = new StringBuilder(1024);

        sb.Append('{');
        AppendString(sb, "type", timeline.Type.ToString());
        sb.Append(',');

        // Dimensions
        sb.Append("\"dimensions\":{");
        AppendInt(sb, "nx", first.Nx); sb.Append(',');
        AppendInt(sb, "ny", first.Ny); sb.Append(',');
        AppendDouble(sb, "dx", first.Dx); sb.Append(',');
        AppendDouble(sb, "dy", first.Dy);
        sb.Append("},");

        AppendMetadata(sb, metadata);

        // Time steps
        sb.Append("\"timeSteps\":[");
        for (int t = 0; t < timeline.Count; t++)
        {
            if (t > 0) sb.Append(',');
            AppendFieldStep(sb, timeline[t]);
        }
        sb.Append("]}");

        return sb.ToString();
    }

    /// <summary>
    /// Save a <see cref="SimulationTimeline"/> to a JSON file.
    /// </summary>
    public static void Save(SimulationTimeline timeline, string path, ExportMetadata metadata = null)
        => File.WriteAllText(path, ToJson(timeline, metadata), Encoding.UTF8);

    // ═══════════════════════════════════════════════════════════════
    //  Single FieldSnapshot
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Serialize a single <see cref="FieldSnapshot"/> to JSON.
    /// </summary>
    public static string ToJson(FieldSnapshot snapshot, ExportMetadata metadata = null)
    {
        var sb = new StringBuilder(512);
        sb.Append('{');
        AppendString(sb, "type", snapshot.Type.ToString());
        sb.Append(',');

        sb.Append("\"dimensions\":{");
        AppendInt(sb, "nx", snapshot.Nx); sb.Append(',');
        AppendInt(sb, "ny", snapshot.Ny); sb.Append(',');
        AppendDouble(sb, "dx", snapshot.Dx); sb.Append(',');
        AppendDouble(sb, "dy", snapshot.Dy);
        sb.Append("},");

        AppendMetadata(sb, metadata);
        AppendFieldStep(sb, snapshot);
        sb.Append('}');

        return sb.ToString();
    }

    /// <summary>
    /// Save a single <see cref="FieldSnapshot"/> to a JSON file.
    /// </summary>
    public static void Save(FieldSnapshot snapshot, string path, ExportMetadata metadata = null)
        => File.WriteAllText(path, ToJson(snapshot, metadata), Encoding.UTF8);

    // ═══════════════════════════════════════════════════════════════
    //  BeamSnapshot
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Serialize a <see cref="BeamSnapshot"/> to JSON.
    /// </summary>
    public static string ToJson(BeamSnapshot beam, ExportMetadata metadata = null)
    {
        var sb = new StringBuilder(512);
        sb.Append('{');
        AppendString(sb, "type", "BeamStress");
        sb.Append(',');
        AppendString(sb, "support", beam.Support.ToString());
        sb.Append(',');
        AppendDouble(sb, "length", beam.Length);
        sb.Append(',');
        AppendInt(sb, "nodeCount", beam.NodeCount);
        sb.Append(',');
        AppendDouble(sb, "maxDeflection", beam.MaxDeflection);
        sb.Append(',');
        AppendDouble(sb, "maxStress", beam.MaxStress);
        sb.Append(',');

        AppendMetadata(sb, metadata);

        AppendArray(sb, "positions", beam.Positions); sb.Append(',');
        AppendArray(sb, "deflection", beam.Deflection); sb.Append(',');
        AppendArray(sb, "bendingMoment", beam.BendingMoment); sb.Append(',');
        AppendArray(sb, "shearForce", beam.ShearForce); sb.Append(',');
        AppendArray(sb, "stress", beam.Stress);

        sb.Append('}');
        return sb.ToString();
    }

    /// <summary>
    /// Save a <see cref="BeamSnapshot"/> to a JSON file.
    /// </summary>
    public static void Save(BeamSnapshot beam, string path, ExportMetadata metadata = null)
        => File.WriteAllText(path, ToJson(beam, metadata), Encoding.UTF8);

    // ═══════════════════════════════════════════════════════════════
    //  SimulationResult (convenience)
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Serialize a <see cref="SimulationResult"/> to JSON.
    /// For 2D results, exports the final field. For beams, exports all curves.
    /// For pipe flow, exports positions + velocity profile.
    /// </summary>
    public static string ToJson(SimulationResult result, ExportMetadata metadata = null)
    {
        var sb = new StringBuilder(512);
        sb.Append('{');
        AppendString(sb, "type", result.Type.ToString());
        sb.Append(',');
        AppendDouble(sb, "maxValue", result.MaxValue);
        sb.Append(',');
        AppendDouble(sb, "minValue", result.MinValue);
        sb.Append(',');
        AppendInt(sb, "iterations", result.Iterations);
        sb.Append(',');
        AppendDouble(sb, "finalTime", result.FinalTime);
        sb.Append(',');

        AppendMetadata(sb, metadata);

        switch (result.Type)
        {
            case MultiphysicsType.HeatPlate:
            case MultiphysicsType.ElectricField:
                if (result.Field != null)
                    AppendField2D(sb, "field", result.Field);
                if (result.Ex != null)
                {
                    sb.Append(',');
                    AppendField2D(sb, "ex", result.Ex);
                    sb.Append(',');
                    AppendField2D(sb, "ey", result.Ey);
                }
                break;

            case MultiphysicsType.BeamStress:
            case MultiphysicsType.PipeFlow:
                if (result.Positions != null)
                    AppendArray(sb, "positions", result.Positions);
                if (result.Values != null)
                {
                    sb.Append(',');
                    AppendArray(sb, "values", result.Values);
                }
                if (result.BendingMoment != null)
                {
                    sb.Append(',');
                    AppendArray(sb, "bendingMoment", result.BendingMoment);
                    sb.Append(',');
                    AppendArray(sb, "shearForce", result.ShearForce);
                    sb.Append(',');
                    AppendArray(sb, "stress", result.Stress);
                }
                break;
        }

        sb.Append('}');
        return sb.ToString();
    }

    /// <summary>Save a <see cref="SimulationResult"/> to a JSON file.</summary>
    public static void Save(SimulationResult result, string path, ExportMetadata metadata = null)
        => File.WriteAllText(path, ToJson(result, metadata), Encoding.UTF8);

    // ═══════════════════════════════════════════════════════════════
    //  Internal helpers
    // ═══════════════════════════════════════════════════════════════

    private static void AppendFieldStep(StringBuilder sb, FieldSnapshot snap)
    {
        sb.Append('{');
        AppendDouble(sb, "time", snap.Time); sb.Append(',');
        AppendInt(sb, "stepIndex", snap.StepIndex); sb.Append(',');
        AppendDouble(sb, "min", snap.Min()); sb.Append(',');
        AppendDouble(sb, "max", snap.Max()); sb.Append(',');

        var values = snap.GetValues();
        sb.Append("\"values\":[");
        for (int i = 0; i < values.Length; i++)
        {
            if (i > 0) sb.Append(',');
            sb.Append(values[i].ToString("G6", Inv));
        }
        sb.Append(']');
        sb.Append('}');
    }

    private static void AppendField2D(StringBuilder sb, string name, double[,] field)
    {
        int nx = field.GetLength(0), ny = field.GetLength(1);
        sb.Append('"').Append(name).Append("\":[");
        bool first = true;
        for (int iy = 0; iy < ny; iy++)
        {
            for (int ix = 0; ix < nx; ix++)
            {
                if (!first) sb.Append(',');
                first = false;
                sb.Append(field[ix, iy].ToString("G6", Inv));
            }
        }
        sb.Append(']');
    }

    private static void AppendArray(StringBuilder sb, string name, double[] arr)
    {
        sb.Append('"').Append(name).Append("\":[");
        for (int i = 0; i < arr.Length; i++)
        {
            if (i > 0) sb.Append(',');
            sb.Append(arr[i].ToString("G6", Inv));
        }
        sb.Append(']');
    }

    private static void AppendString(StringBuilder sb, string key, string value)
    {
        sb.Append('"').Append(key).Append("\":\"").Append(value).Append('"');
    }

    private static void AppendInt(StringBuilder sb, string key, int value)
    {
        sb.Append('"').Append(key).Append("\":").Append(value);
    }

    private static void AppendDouble(StringBuilder sb, string key, double value)
    {
        sb.Append('"').Append(key).Append("\":").Append(value.ToString("G6", Inv));
    }

    private static void AppendMetadata(StringBuilder sb, ExportMetadata metadata)
    {
        if (metadata != null)
        {
            sb.Append("\"metadata\":{");
            bool first = true;
            if (metadata.Simulation != null)
            {
                AppendString(sb, "simulation", metadata.Simulation);
                first = false;
            }
            if (metadata.Unit != null)
            {
                if (!first) sb.Append(',');
                AppendString(sb, "unit", metadata.Unit);
                first = false;
            }
            if (!double.IsNaN(metadata.MaterialDensity))
            {
                if (!first) sb.Append(',');
                AppendDouble(sb, "density", metadata.MaterialDensity);
            }
            sb.Append("},");
        }
    }
}

/// <summary>Optional metadata attached to JSON exports.</summary>
public class ExportMetadata
{
    /// <summary>Simulation description (e.g. "Steel cantilever beam, 1 kN point load").</summary>
    public string Simulation { get; set; }

    /// <summary>Unit of the primary field (e.g. "K", "V", "m", "m/s").</summary>
    public string Unit { get; set; }

    /// <summary>Material density used (kg/m³). NaN if not set.</summary>
    public double MaterialDensity { get; set; } = double.NaN;
}
