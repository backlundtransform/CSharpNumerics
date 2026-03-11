using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;

namespace CSharpNumerics.Engines.Common
{
    /// <summary>
    /// Serializes and deserializes engine state (arrays of doubles) to binary or CSV.
    /// </summary>
    public static class FieldSerializer
    {
        #region Binary

        /// <summary>Write a double array to a binary file.</summary>
        public static void SaveBinary(double[] data, string path)
        {
            using var fs = new FileStream(path, FileMode.Create, FileAccess.Write);
            using var bw = new BinaryWriter(fs);
            bw.Write(data.Length);
            foreach (var d in data)
                bw.Write(d);
        }

        /// <summary>Read a double array from a binary file.</summary>
        public static double[] LoadBinary(string path)
        {
            using var fs = new FileStream(path, FileMode.Open, FileAccess.Read);
            using var br = new BinaryReader(fs);
            int len = br.ReadInt32();
            var data = new double[len];
            for (int i = 0; i < len; i++)
                data[i] = br.ReadDouble();
            return data;
        }

        #endregion

        #region CSV

        /// <summary>Write a double array as one-value-per-line CSV.</summary>
        public static void SaveCsv(double[] data, string path)
        {
            var sb = new StringBuilder();
            sb.AppendLine("Value");
            foreach (var d in data)
                sb.AppendLine(d.ToString(CultureInfo.InvariantCulture));
            File.WriteAllText(path, sb.ToString(), Encoding.UTF8);
        }

        /// <summary>Read a double array from a one-value-per-line CSV (first line = header).</summary>
        public static double[] LoadCsv(string path)
        {
            var lines = File.ReadAllLines(path, Encoding.UTF8);
            var result = new List<double>();
            for (int i = 1; i < lines.Length; i++)
            {
                var line = lines[i].Trim();
                if (line.Length > 0)
                    result.Add(double.Parse(line, CultureInfo.InvariantCulture));
            }
            return result.ToArray();
        }

        #endregion

        #region JSON (minimal, no external deps)

        /// <summary>Serialize a double array to a minimal JSON string: [v0, v1, …].</summary>
        public static string ToJson(double[] data)
        {
            var sb = new StringBuilder("[");
            for (int i = 0; i < data.Length; i++)
            {
                if (i > 0) sb.Append(',');
                sb.Append(data[i].ToString("R", CultureInfo.InvariantCulture));
            }
            sb.Append(']');
            return sb.ToString();
        }

        /// <summary>Deserialize a JSON array string back to double[].</summary>
        public static double[] FromJson(string json)
        {
            json = json.Trim();
            if (json.StartsWith("[")) json = json.Substring(1);
            if (json.EndsWith("]")) json = json.Substring(0, json.Length - 1);
            var parts = json.Split(',');
            var result = new List<double>();
            foreach (var p in parts)
            {
                var trimmed = p.Trim();
                if (trimmed.Length > 0)
                    result.Add(double.Parse(trimmed, CultureInfo.InvariantCulture));
            }
            return result.ToArray();
        }

        #endregion
    }
}
