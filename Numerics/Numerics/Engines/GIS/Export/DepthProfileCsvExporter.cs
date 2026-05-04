using CSharpNumerics.Engines.GIS.Spread.VolumetricContamination;
using System.Globalization;
using System.IO;
using System.Text;

namespace CSharpNumerics.Engines.GIS.Export;

/// <summary>
/// Exports volumetric contamination depth profiles to CSV.
/// </summary>
public static class DepthProfileCsvExporter
{
    /// <summary>
    /// Exports a depth profile at horizontal position (ix, iy) to CSV.
    /// Columns: Depth_m, Concentration_mgL.
    /// </summary>
    public static string ToCsv(VolumetricContaminationResult result, int ix, int iy, int timeIndex = -1)
    {
        var profile = result.GetDepthProfile(ix, iy, timeIndex);
        var sb = new StringBuilder();
        sb.AppendLine("Depth_m,Concentration_mgL");
        for (int i = 0; i < profile.Depths.Length; i++)
        {
            sb.Append(profile.Depths[i].ToString("R", CultureInfo.InvariantCulture));
            sb.Append(',');
            sb.AppendLine(profile.Concentrations[i].ToString("R", CultureInfo.InvariantCulture));
        }
        return sb.ToString();
    }

    /// <summary>Saves a depth profile CSV to a file.</summary>
    public static void Save(VolumetricContaminationResult result, int ix, int iy, string path, int timeIndex = -1)
        => File.WriteAllText(path, ToCsv(result, ix, iy, timeIndex), Encoding.UTF8);
}
