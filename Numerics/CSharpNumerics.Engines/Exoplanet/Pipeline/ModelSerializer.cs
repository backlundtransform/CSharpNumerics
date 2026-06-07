using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Reflection;
using System.Text;
using CSharpNumerics.Engines.Common;

namespace CSharpNumerics.Engines.Exoplanet.Pipeline;

/// <summary>
/// Serializes and deserializes a <see cref="TrainedTransitModel"/> to/from a byte array or file.
/// Format: JSON metadata (config + metrics) + binary weight data extracted via reflection
/// from the trained pipeline's neural network layers.
/// </summary>
public static class ModelSerializer
{
    private const int MagicHeader = 0x45584F50; // "EXOP"
    private const int FormatVersion = 1;

    /// <summary>
    /// Serializes a <see cref="TrainedTransitModel"/> to a byte array.
    /// The metadata (configuration, model name, metrics) is stored as JSON.
    /// The model weights are extracted via reflection and stored as binary double arrays.
    /// </summary>
    public static byte[] Serialize(TrainedTransitModel model)
    {
        if (model == null) throw new ArgumentNullException(nameof(model));

        using var ms = new MemoryStream();
        using var bw = new BinaryWriter(ms, Encoding.UTF8, leaveOpen: true);

        // Header
        bw.Write(MagicHeader);
        bw.Write(FormatVersion);

        // Metadata JSON
        string metadataJson = SerializeMetadata(model);
        byte[] metadataBytes = Encoding.UTF8.GetBytes(metadataJson);
        bw.Write(metadataBytes.Length);
        bw.Write(metadataBytes);

        // Extract and write all weights from the pipeline's model
        var weights = ExtractWeights(model.Pipeline.Model);
        bw.Write(weights.Count);
        foreach (var (path, data) in weights)
        {
            byte[] pathBytes = Encoding.UTF8.GetBytes(path);
            bw.Write(pathBytes.Length);
            bw.Write(pathBytes);
            bw.Write(data.Length);
            foreach (double d in data)
                bw.Write(d);
        }

        bw.Flush();
        return ms.ToArray();
    }

    /// <summary>
    /// Deserializes a <see cref="TrainedTransitModel"/> from a byte array.
    /// The pipeline is reconstructed by cloning, fitting with dummy data, then restoring weights.
    /// Note: Requires the same model type to be available at runtime.
    /// </summary>
    public static TrainedTransitModel Deserialize(byte[] data, TrainedTransitModel template)
    {
        if (data == null) throw new ArgumentNullException(nameof(data));
        if (template == null) throw new ArgumentNullException(nameof(template));

        using var ms = new MemoryStream(data);
        using var br = new BinaryReader(ms, Encoding.UTF8);

        // Header
        int magic = br.ReadInt32();
        if (magic != MagicHeader) throw new InvalidDataException("Invalid model file format.");
        int version = br.ReadInt32();
        if (version != FormatVersion) throw new InvalidDataException($"Unsupported format version: {version}");

        // Metadata JSON
        int metadataLen = br.ReadInt32();
        byte[] metadataBytes = br.ReadBytes(metadataLen);
        string metadataJson = Encoding.UTF8.GetString(metadataBytes);
        var metadata = DeserializeMetadata(metadataJson);

        // Read weights
        var weights = new List<(string path, double[] data)>();
        int weightCount = br.ReadInt32();
        for (int i = 0; i < weightCount; i++)
        {
            int pathLen = br.ReadInt32();
            byte[] pathBytes = br.ReadBytes(pathLen);
            string path = Encoding.UTF8.GetString(pathBytes);
            int arrayLen = br.ReadInt32();
            double[] arr = new double[arrayLen];
            for (int j = 0; j < arrayLen; j++)
                arr[j] = br.ReadDouble();
            weights.Add((path, arr));
        }

        // Clone the template pipeline and restore weights
        var pipeline = template.Pipeline.Clone();

        // Fit with dummy data to initialize the model's internal structure
        int windowSize = template.TrainerCfg.WindowSize;
        int numFeatures = CSharpNumerics.Engines.Exoplanet.Data.TransitFeatureSet.FeatureNames.All.Length;
        int totalCols = windowSize + numFeatures;
        var dummyX = new CSharpNumerics.Numerics.Objects.Matrix(2, totalCols);
        var dummyY = new CSharpNumerics.Numerics.Objects.VectorN(new double[] { 0.0, 1.0 });
        pipeline.Fit(dummyX, dummyY);

        // Restore weights
        RestoreWeights(pipeline.Model, weights);

        return new TrainedTransitModel(
            pipeline,
            template.DetectionConfig,
            template.TrainerCfg,
            metadata.modelName,
            metadata.metrics);
    }

    /// <summary>
    /// Saves a trained model to a file.
    /// </summary>
    public static void SaveToFile(TrainedTransitModel model, string path)
    {
        if (string.IsNullOrEmpty(path)) throw new ArgumentNullException(nameof(path));
        byte[] bytes = Serialize(model);
        File.WriteAllBytes(path, bytes);
    }

    /// <summary>
    /// Loads a trained model from a file.
    /// </summary>
    public static TrainedTransitModel LoadFromFile(string path, TrainedTransitModel template)
    {
        if (string.IsNullOrEmpty(path)) throw new ArgumentNullException(nameof(path));
        byte[] bytes = File.ReadAllBytes(path);
        return Deserialize(bytes, template);
    }

    #region Weight Extraction / Restoration

    private static List<(string path, double[] data)> ExtractWeights(object model)
    {
        var result = new List<(string, double[])>();
        ExtractWeightsRecursive(model, "", result, new HashSet<object>());
        return result;
    }

    private static void ExtractWeightsRecursive(object obj, string prefix,
        List<(string, double[])> result, HashSet<object> visited)
    {
        if (obj == null || visited.Contains(obj)) return;
        visited.Add(obj);

        var type = obj.GetType();
        var flags = BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic;

        foreach (var field in type.GetFields(flags))
        {
            string fieldPath = string.IsNullOrEmpty(prefix) ? field.Name : $"{prefix}.{field.Name}";
            var value = field.GetValue(obj);
            if (value == null) continue;

            if (value is double[] arr)
            {
                result.Add((fieldPath, (double[])arr.Clone()));
            }
            else if (value is double[,] arr2d)
            {
                int rows = arr2d.GetLength(0);
                int cols = arr2d.GetLength(1);
                double[] flat = new double[rows * cols + 2];
                flat[0] = rows;
                flat[1] = cols;
                for (int r = 0; r < rows; r++)
                    for (int c = 0; c < cols; c++)
                        flat[2 + r * cols + c] = arr2d[r, c];
                result.Add(($"{fieldPath}:2d", flat));
            }
            else if (value is Array objArr && !field.FieldType.IsPrimitive)
            {
                for (int i = 0; i < objArr.Length; i++)
                {
                    var elem = objArr.GetValue(i);
                    if (elem != null && !elem.GetType().IsPrimitive && !elem.GetType().IsEnum)
                        ExtractWeightsRecursive(elem, $"{fieldPath}[{i}]", result, visited);
                }
            }
            else if (!field.FieldType.IsPrimitive && !field.FieldType.IsEnum
                     && !field.FieldType.IsValueType && field.FieldType != typeof(string))
            {
                ExtractWeightsRecursive(value, fieldPath, result, visited);
            }
        }
    }

    private static void RestoreWeights(object model, List<(string path, double[] data)> weights)
    {
        var weightMap = new Dictionary<string, double[]>();
        foreach (var (path, data) in weights)
            weightMap[path] = data;

        RestoreWeightsRecursive(model, "", weightMap, new HashSet<object>());
    }

    private static void RestoreWeightsRecursive(object obj, string prefix,
        Dictionary<string, double[]> weightMap, HashSet<object> visited)
    {
        if (obj == null || visited.Contains(obj)) return;
        visited.Add(obj);

        var type = obj.GetType();
        var flags = BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic;

        foreach (var field in type.GetFields(flags))
        {
            string fieldPath = string.IsNullOrEmpty(prefix) ? field.Name : $"{prefix}.{field.Name}";
            var value = field.GetValue(obj);

            if (value is double[] existingArr && weightMap.TryGetValue(fieldPath, out var storedArr))
            {
                Array.Copy(storedArr, existingArr, Math.Min(existingArr.Length, storedArr.Length));
            }
            else if (value is double[,] existing2d && weightMap.TryGetValue($"{fieldPath}:2d", out var stored2d))
            {
                int rows = (int)stored2d[0];
                int cols = (int)stored2d[1];
                for (int r = 0; r < rows && r < existing2d.GetLength(0); r++)
                    for (int c = 0; c < cols && c < existing2d.GetLength(1); c++)
                        existing2d[r, c] = stored2d[2 + r * cols + c];
            }
            else if (value is Array objArr && !field.FieldType.IsPrimitive && value != null)
            {
                for (int i = 0; i < objArr.Length; i++)
                {
                    var elem = objArr.GetValue(i);
                    if (elem != null && !elem.GetType().IsPrimitive && !elem.GetType().IsEnum)
                        RestoreWeightsRecursive(elem, $"{fieldPath}[{i}]", weightMap, visited);
                }
            }
            else if (value != null && !field.FieldType.IsPrimitive && !field.FieldType.IsEnum
                     && !field.FieldType.IsValueType && field.FieldType != typeof(string))
            {
                RestoreWeightsRecursive(value, fieldPath, weightMap, visited);
            }
        }
    }

    #endregion

    #region Metadata Serialization (minimal JSON, no external deps)

    private static string SerializeMetadata(TrainedTransitModel model)
    {
        var sb = new StringBuilder();
        sb.Append('{');
        sb.Append($"\"modelName\":\"{EscapeJson(model.ModelName)}\"");
        sb.Append($",\"windowSize\":{model.TrainerCfg.WindowSize}");
        sb.Append($",\"phaseBins\":{model.TrainerCfg.PhaseBins}");
        sb.Append($",\"accuracy\":{model.Metrics.Accuracy.ToString("R", CultureInfo.InvariantCulture)}");
        sb.Append($",\"precision\":{model.Metrics.Precision.ToString("R", CultureInfo.InvariantCulture)}");
        sb.Append($",\"recall\":{model.Metrics.Recall.ToString("R", CultureInfo.InvariantCulture)}");
        sb.Append($",\"f1\":{model.Metrics.F1Score.ToString("R", CultureInfo.InvariantCulture)}");
        sb.Append($",\"bestCvScore\":{model.Metrics.BestCvScore.ToString("R", CultureInfo.InvariantCulture)}");
        sb.Append('}');
        return sb.ToString();
    }

    private static (string modelName, TrainingMetrics metrics) DeserializeMetadata(string json)
    {
        var metrics = new TrainingMetrics();
        string modelName = "";

        json = json.Trim().TrimStart('{').TrimEnd('}');
        var pairs = SplitJsonPairs(json);

        foreach (var (key, val) in pairs)
        {
            switch (key)
            {
                case "modelName": modelName = val.Trim('"'); break;
                case "accuracy": metrics.Accuracy = ParseDouble(val); break;
                case "precision": metrics.Precision = ParseDouble(val); break;
                case "recall": metrics.Recall = ParseDouble(val); break;
                case "f1": metrics.F1Score = ParseDouble(val); break;
                case "bestCvScore": metrics.BestCvScore = ParseDouble(val); break;
            }
        }

        return (modelName, metrics);
    }

    private static List<(string key, string value)> SplitJsonPairs(string json)
    {
        var result = new List<(string, string)>();
        int depth = 0;
        bool inString = false;
        int start = 0;

        for (int i = 0; i <= json.Length; i++)
        {
            char c = i < json.Length ? json[i] : ',';
            if (c == '"' && (i == 0 || json[i - 1] != '\\'))
                inString = !inString;
            else if (!inString && (c == '{' || c == '['))
                depth++;
            else if (!inString && (c == '}' || c == ']'))
                depth--;
            else if (!inString && depth == 0 && c == ',')
            {
                string pair = json.Substring(start, i - start).Trim();
                int colonIdx = pair.IndexOf(':');
                if (colonIdx > 0)
                {
                    string key = pair.Substring(0, colonIdx).Trim().Trim('"');
                    string val = pair.Substring(colonIdx + 1).Trim();
                    result.Add((key, val));
                }
                start = i + 1;
            }
        }

        return result;
    }

    private static double ParseDouble(string s)
    {
        s = s.Trim().Trim('"');
        if (double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture, out double d))
            return d;
        return 0.0;
    }

    private static string EscapeJson(string s)
    {
        return s?.Replace("\\", "\\\\").Replace("\"", "\\\"") ?? "";
    }

    #endregion
}
