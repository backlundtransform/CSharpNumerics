
namespace NumericTest.TestData;

public static class CsvTestDataGenerator
{
    private static Random _rand = new Random(42);
    /// <summary>
    /// Generates a simple linear time series CSV for Rolling CV tests.
    /// Timestamp, Feature1, Target
    /// </summary>
    public static void GenerateTimeSeriesCsv(string path, int n = 30)
    {
        using var sw = new StreamWriter(path);
        sw.WriteLine("Timestamp,Feature1,Target");

        for (int i = 0; i < n; i++)
        {
            double feature = i / 10.0;
            double target = 2 * feature + 1; // simple linear relation
            sw.WriteLine($"2026-01-{i + 1:00},{feature},{target}");
        }
    }

    /// <summary>
    /// Generates tabular/grouped CSV for Leave-One-Out or Grouped CV tests.
    /// Columns: Group, Feature1, Target
    /// </summary>
    public static void GenerateGroupedCsv(string path, int groups = 5, int samplesPerGroup = 5)
    {
        using var sw = new StreamWriter(path);
        sw.WriteLine("Group,Feature1,Target");

        var rnd = new Random(42);
        for (int g = 0; g < groups; g++)
        {
            for (int s = 0; s < samplesPerGroup; s++)
            {
                double feature = s + rnd.NextDouble();       // small noise
                double target = g * 10 + feature;           // unique group bias
                sw.WriteLine($"{g},{feature},{target}");
            }
        }
    }

    public static void GenerateClassificationCsv(
           string filePath,
           int nSamples = 100,
           int nFeatures = 3,
           int nClasses = 2,
           double imbalanceRatio = 0.5)
    {
        if (nClasses != 2)
            throw new NotImplementedException("Currently only binary classification supported.");

        using var writer = new StreamWriter(filePath);

        // header
        for (int f = 1; f <= nFeatures; f++)
            writer.Write($"Feature{f},");
        writer.WriteLine("Target");

        int nClass0 = (int)(nSamples * imbalanceRatio);
        int nClass1 = nSamples - nClass0;

        for (int i = 0; i < nClass0; i++)
        {
            string line = "";
            for (int f = 0; f < nFeatures; f++)
                line += $"{_rand.NextDouble() * 10:F2},";
            line += "0";
            writer.WriteLine(line);
        }

        for (int i = 0; i < nClass1; i++)
        {
            string line = "";
            for (int f = 0; f < nFeatures; f++)
                line += $"{_rand.NextDouble() * 10:F2},";
            line += "1";
            writer.WriteLine(line);
        }
    }
}
