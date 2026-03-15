using CSharpNumerics.Engines.GIS.Analysis;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using System;
using System.IO;

namespace CSharpNumerics.Engines.GIS.Export
{
    /// <summary>
    /// Exports grid data to a compact binary format designed for direct loading
    /// in Unity3D via <c>NativeArray&lt;float&gt;</c>.
    /// <para>
    /// Format:
    /// <list type="bullet">
    ///   <item>Header: magic "GPLM", version, grid dims, time count, axis bounds, time range.</item>
    ///   <item>Body: per time step — float[cells] concentration + float[cells] probability.</item>
    /// </list>
    /// </para>
    /// </summary>
    public static class UnityBinaryExporter
    {
        private static readonly byte[] Magic = { (byte)'G', (byte)'P', (byte)'L', (byte)'M' };
        private const int Version = 1;

        // ═══════════════════════════════════════════════════════════════
        //  Save — concentration snapshots only (deterministic)
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Saves a time series of <see cref="GridSnapshot"/> (no probability).
        /// The probability layer is written as zeros.
        /// </summary>
        public static void Save(
            GridSnapshot[] snapshots,
            GeoGrid grid,
            TimeFrame timeFrame,
            string path)
        {
            if (snapshots == null) throw new ArgumentNullException(nameof(snapshots));
            if (grid == null) throw new ArgumentNullException(nameof(grid));
            if (timeFrame == null) throw new ArgumentNullException(nameof(timeFrame));

            using var fs = new FileStream(path, FileMode.Create, FileAccess.Write);
            using var bw = new BinaryWriter(fs);

            WriteHeader(bw, grid, timeFrame, snapshots.Length);

            var zeroes = new float[grid.CellCount];
            foreach (var snap in snapshots)
            {
                WriteFloatLayer(bw, snap.GetValues());
                WriteFloatLayerRaw(bw, zeroes);
            }
        }

        // ═══════════════════════════════════════════════════════════════
        //  Save — concentration + probability (TimeAnimator)
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Saves concentration snapshots paired with probability maps.
        /// </summary>
        public static void Save(
            GridSnapshot[] concentrationSnapshots,
            TimeAnimator animation,
            string path)
        {
            if (concentrationSnapshots == null) throw new ArgumentNullException(nameof(concentrationSnapshots));
            if (animation == null) throw new ArgumentNullException(nameof(animation));

            var grid = animation.Grid;
            var timeFrame = animation.TimeFrame;
            int timeCount = Math.Min(concentrationSnapshots.Length, animation.Count);

            using var fs = new FileStream(path, FileMode.Create, FileAccess.Write);
            using var bw = new BinaryWriter(fs);

            WriteHeader(bw, grid, timeFrame, timeCount);

            for (int t = 0; t < timeCount; t++)
            {
                WriteFloatLayer(bw, concentrationSnapshots[t].GetValues());
                WriteFloatLayer(bw, animation[t].GetValues());
            }
        }

        // ═══════════════════════════════════════════════════════════════
        //  Save — probability only (TimeAnimator)
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Saves probability maps (concentration layer written as zeros).
        /// </summary>
        public static void Save(TimeAnimator animation, string path)
        {
            if (animation == null) throw new ArgumentNullException(nameof(animation));

            var grid = animation.Grid;
            var timeFrame = animation.TimeFrame;

            using var fs = new FileStream(path, FileMode.Create, FileAccess.Write);
            using var bw = new BinaryWriter(fs);

            WriteHeader(bw, grid, timeFrame, animation.Count);

            var zeroes = new float[grid.CellCount];
            for (int t = 0; t < animation.Count; t++)
            {
                WriteFloatLayerRaw(bw, zeroes);
                WriteFloatLayer(bw, animation[t].GetValues());
            }
        }

        // ═══════════════════════════════════════════════════════════════
        //  Read — load header + data
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Reads the binary header from a Unity export file.
        /// </summary>
        public static UnityBinaryHeader ReadHeader(string path)
        {
            using var fs = new FileStream(path, FileMode.Open, FileAccess.Read);
            using var br = new BinaryReader(fs);
            return ReadHeaderInternal(br);
        }

        /// <summary>
        /// Reads both concentration and probability layers from a Unity export file.
        /// </summary>
        public static UnityBinaryData Read(string path)
        {
            using var fs = new FileStream(path, FileMode.Open, FileAccess.Read);
            using var br = new BinaryReader(fs);

            var header = ReadHeaderInternal(br);
            int cellCount = header.Nx * header.Ny * header.Nz;

            var concentration = new float[header.TimeStepCount][];
            var probability = new float[header.TimeStepCount][];

            for (int t = 0; t < header.TimeStepCount; t++)
            {
                concentration[t] = ReadFloatLayer(br, cellCount);
                probability[t] = ReadFloatLayer(br, cellCount);
            }

            return new UnityBinaryData(header, concentration, probability);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Internal — header
        // ═══════════════════════════════════════════════════════════════

        private static void WriteHeader(BinaryWriter bw, GeoGrid grid, TimeFrame timeFrame, int timeStepCount)
        {
            // Magic
            bw.Write(Magic);
            // Version
            bw.Write(Version);
            // Grid dims
            bw.Write(grid.Nx);
            bw.Write(grid.Ny);
            bw.Write(grid.Nz);
            // Time step count
            bw.Write(timeStepCount);
            // Axis bounds (6 doubles)
            bw.Write(grid.XMin);
            bw.Write(grid.XMax);
            bw.Write(grid.YMin);
            bw.Write(grid.YMax);
            bw.Write(grid.ZMin);
            bw.Write(grid.ZMax);
            // Time range (3 doubles)
            bw.Write(timeFrame.Start);
            bw.Write(timeFrame.End);
            bw.Write(timeFrame.StepSeconds);
        }

        private static UnityBinaryHeader ReadHeaderInternal(BinaryReader br)
        {
            var magic = br.ReadBytes(4);
            if (magic.Length < 4 ||
                magic[0] != Magic[0] || magic[1] != Magic[1] ||
                magic[2] != Magic[2] || magic[3] != Magic[3])
                throw new InvalidDataException("Invalid magic bytes — not a GPLM file.");

            int version = br.ReadInt32();
            if (version != Version)
                throw new InvalidDataException($"Unsupported version {version}.");

            int nx = br.ReadInt32();
            int ny = br.ReadInt32();
            int nz = br.ReadInt32();
            int timeStepCount = br.ReadInt32();
            double xMin = br.ReadDouble();
            double xMax = br.ReadDouble();
            double yMin = br.ReadDouble();
            double yMax = br.ReadDouble();
            double zMin = br.ReadDouble();
            double zMax = br.ReadDouble();
            double tStart = br.ReadDouble();
            double tEnd = br.ReadDouble();
            double tStep = br.ReadDouble();

            return new UnityBinaryHeader(
                nx, ny, nz, timeStepCount,
                xMin, xMax, yMin, yMax, zMin, zMax,
                tStart, tEnd, tStep);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Internal — layers
        // ═══════════════════════════════════════════════════════════════

        private static void WriteFloatLayer(BinaryWriter bw, double[] values)
        {
            for (int i = 0; i < values.Length; i++)
                bw.Write((float)values[i]);
        }

        private static void WriteFloatLayerRaw(BinaryWriter bw, float[] values)
        {
            for (int i = 0; i < values.Length; i++)
                bw.Write(values[i]);
        }

        private static float[] ReadFloatLayer(BinaryReader br, int count)
        {
            var data = new float[count];
            for (int i = 0; i < count; i++)
                data[i] = br.ReadSingle();
            return data;
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Header & Data types
    // ═══════════════════════════════════════════════════════════════════

    /// <summary>
    /// Parsed header from a Unity binary export file.
    /// </summary>
    public class UnityBinaryHeader
    {
        public int Nx { get; }
        public int Ny { get; }
        public int Nz { get; }
        public int TimeStepCount { get; }
        public double XMin { get; }
        public double XMax { get; }
        public double YMin { get; }
        public double YMax { get; }
        public double ZMin { get; }
        public double ZMax { get; }
        public double TStart { get; }
        public double TEnd { get; }
        public double TStep { get; }

        internal UnityBinaryHeader(
            int nx, int ny, int nz, int timeStepCount,
            double xMin, double xMax, double yMin, double yMax,
            double zMin, double zMax,
            double tStart, double tEnd, double tStep)
        {
            Nx = nx; Ny = ny; Nz = nz;
            TimeStepCount = timeStepCount;
            XMin = xMin; XMax = xMax;
            YMin = yMin; YMax = yMax;
            ZMin = zMin; ZMax = zMax;
            TStart = tStart; TEnd = tEnd; TStep = tStep;
        }
    }

    /// <summary>
    /// Full contents of a Unity binary export: header + per-timestep layers.
    /// </summary>
    public class UnityBinaryData
    {
        public UnityBinaryHeader Header { get; }

        /// <summary>Concentration layers: [timeStep][cellIndex].</summary>
        public float[][] Concentration { get; }

        /// <summary>Probability layers: [timeStep][cellIndex].</summary>
        public float[][] Probability { get; }

        internal UnityBinaryData(UnityBinaryHeader header, float[][] concentration, float[][] probability)
        {
            Header = header;
            Concentration = concentration;
            Probability = probability;
        }
    }
}
