using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Engines.Multiphysics.Snapshots;
using System;
using System.IO;

namespace CSharpNumerics.Engines.Multiphysics.Export;

/// <summary>
/// Exports multiphysics simulation data in a compact binary format
/// optimised for Unity NativeArray loading.
/// <para>
/// Format (little-endian):
/// <code>
/// [4 bytes] Magic: "MPHY"
/// [int]     Version: 1
/// [int]     SimulationType (enum ordinal)
/// [int]     Nx (or NodeCount for 1D)
/// [int]     Ny (0 for 1D data)
/// [double]  Dx
/// [double]  Dy
/// [int]     TimeStepCount
/// [double]  Dt
/// [double]  TotalTime
/// [int]     LayerCount (number of named layers per step)
/// For 2D: per time step × per layer → float[Nx*Ny]
/// For 1D: float[Nx] positions + per layer → float[Nx] values
/// </code>
/// </para>
/// </summary>
public static class MultiphysicsBinaryExporter
{
    private static readonly byte[] Magic = { (byte)'M', (byte)'P', (byte)'H', (byte)'Y' };
    private const int Version = 1;

    // ═══════════════════════════════════════════════════════════════
    //  Save — 2D timeline (HeatPlate, ElectricField)
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Save a <see cref="SimulationTimeline"/> to binary.
    /// Each time step stores one float layer (the primary field).
    /// </summary>
    public static void Save(SimulationTimeline timeline, string path)
    {
        if (timeline.Count == 0)
            throw new ArgumentException("Timeline is empty.");

        var first = timeline[0];
        int nx = first.Nx, ny = first.Ny;
        int cellCount = nx * ny;

        using var fs = new FileStream(path, FileMode.Create, FileAccess.Write);
        using var w = new BinaryWriter(fs);

        // Header
        w.Write(Magic);
        w.Write(Version);
        w.Write((int)timeline.Type);
        w.Write(nx);
        w.Write(ny);
        w.Write(first.Dx);
        w.Write(first.Dy);
        w.Write(timeline.Count);
        w.Write(timeline.Dt);
        w.Write(timeline.EndTime);
        w.Write(1); // LayerCount = 1 (primary field)

        // Body: per time step → float[cellCount]
        for (int t = 0; t < timeline.Count; t++)
        {
            var snap = timeline[t];
            var values = snap.GetValues();
            for (int i = 0; i < cellCount; i++)
                w.Write((float)values[i]);
        }
    }

    // ═══════════════════════════════════════════════════════════════
    //  Save — 1D beam snapshot
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Save a <see cref="BeamSnapshot"/> to binary.
    /// Stores positions + 4 layers (deflection, moment, shear, stress).
    /// </summary>
    public static void Save(BeamSnapshot beam, string path)
    {
        int n = beam.NodeCount;

        using var fs = new FileStream(path, FileMode.Create, FileAccess.Write);
        using var w = new BinaryWriter(fs);

        // Header
        w.Write(Magic);
        w.Write(Version);
        w.Write((int)MultiphysicsType.BeamStress);
        w.Write(n);
        w.Write(0);             // Ny = 0 → 1D
        w.Write(beam.Length / (n - 1)); // Dx
        w.Write(0.0);           // Dy = 0
        w.Write(1);             // TimeStepCount = 1 (static)
        w.Write(0.0);           // Dt = 0
        w.Write(0.0);           // TotalTime = 0
        w.Write(4);             // LayerCount = 4

        // Positions
        for (int i = 0; i < n; i++)
            w.Write((float)beam.Positions[i]);

        // Layer 0: deflection
        for (int i = 0; i < n; i++)
            w.Write((float)beam.Deflection[i]);

        // Layer 1: bending moment
        for (int i = 0; i < n; i++)
            w.Write((float)beam.BendingMoment[i]);

        // Layer 2: shear force
        for (int i = 0; i < n; i++)
            w.Write((float)beam.ShearForce[i]);

        // Layer 3: stress
        for (int i = 0; i < n; i++)
            w.Write((float)beam.Stress[i]);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Save — CylinderFlow result (vx, vy, pressure, vorticity)
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Save a CylinderFlow <see cref="SimulationResult"/> to binary.
    /// Stores 4 layers: vx, vy, pressure, vorticity.
    /// </summary>
    public static void Save(SimulationResult result, string path)
    {
        if (result.Vx == null)
            throw new ArgumentException("Result does not contain velocity fields.");

        int nx = result.Vx.GetLength(0);
        int ny = result.Vx.GetLength(1);
        int cellCount = nx * ny;

        using var fs = new FileStream(path, FileMode.Create, FileAccess.Write);
        using var w = new BinaryWriter(fs);

        // Header
        w.Write(Magic);
        w.Write(Version);
        w.Write((int)result.Type);
        w.Write(nx);
        w.Write(ny);
        w.Write(0.0); // Dx unknown from result alone
        w.Write(0.0); // Dy
        w.Write(1);   // TimeStepCount = 1 (final frame)
        w.Write(0.0); // Dt
        w.Write(result.FinalTime);
        w.Write(4);   // LayerCount = 4 (vx, vy, p, vorticity)

        // Layer 0: vx
        WriteField(w, result.Vx, nx, ny);
        // Layer 1: vy
        WriteField(w, result.Vy, nx, ny);
        // Layer 2: pressure
        WriteField(w, result.Pressure, nx, ny);
        // Layer 3: vorticity
        WriteField(w, result.Vorticity, nx, ny);
    }

    private static void WriteField(BinaryWriter w, double[,] field, int nx, int ny)
    {
        for (int iy = 0; iy < ny; iy++)
            for (int ix = 0; ix < nx; ix++)
                w.Write((float)field[ix, iy]);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Read
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Read the header from a binary export file.
    /// </summary>
    public static MultiphysicsBinaryHeader ReadHeader(string path)
    {
        using var fs = new FileStream(path, FileMode.Open, FileAccess.Read);
        using var r = new BinaryReader(fs);
        return ParseHeader(r);
    }

    /// <summary>
    /// Read the complete binary export (header + all data) from a file.
    /// </summary>
    public static MultiphysicsBinaryData Read(string path)
    {
        using var fs = new FileStream(path, FileMode.Open, FileAccess.Read);
        using var r = new BinaryReader(fs);

        var header = ParseHeader(r);
        int cellCount = header.Ny > 0 ? header.Nx * header.Ny : header.Nx;

        float[] positions = null;
        if (header.Ny == 0) // 1D: read positions
        {
            positions = new float[header.Nx];
            for (int i = 0; i < header.Nx; i++)
                positions[i] = r.ReadSingle();
        }

        // Read layers
        var layers = new float[header.TimeStepCount * header.LayerCount][];
        int idx = 0;
        for (int t = 0; t < header.TimeStepCount; t++)
        {
            for (int layer = 0; layer < header.LayerCount; layer++)
            {
                var data = new float[cellCount];
                for (int i = 0; i < cellCount; i++)
                    data[i] = r.ReadSingle();
                layers[idx++] = data;
            }
        }

        return new MultiphysicsBinaryData
        {
            Header = header,
            Positions = positions,
            Layers = layers
        };
    }

    private static MultiphysicsBinaryHeader ParseHeader(BinaryReader r)
    {
        var magic = r.ReadBytes(4);
        if (magic[0] != Magic[0] || magic[1] != Magic[1] ||
            magic[2] != Magic[2] || magic[3] != Magic[3])
            throw new InvalidDataException("Invalid MPHY magic bytes.");

        return new MultiphysicsBinaryHeader
        {
            Version = r.ReadInt32(),
            Type = (MultiphysicsType)r.ReadInt32(),
            Nx = r.ReadInt32(),
            Ny = r.ReadInt32(),
            Dx = r.ReadDouble(),
            Dy = r.ReadDouble(),
            TimeStepCount = r.ReadInt32(),
            Dt = r.ReadDouble(),
            TotalTime = r.ReadDouble(),
            LayerCount = r.ReadInt32()
        };
    }
}

/// <summary>Parsed header from a multiphysics binary export.</summary>
public class MultiphysicsBinaryHeader
{
    public int Version { get; internal set; }
    public MultiphysicsType Type { get; internal set; }
    public int Nx { get; internal set; }
    public int Ny { get; internal set; }
    public double Dx { get; internal set; }
    public double Dy { get; internal set; }
    public int TimeStepCount { get; internal set; }
    public double Dt { get; internal set; }
    public double TotalTime { get; internal set; }
    public int LayerCount { get; internal set; }
    public bool Is2D => Ny > 0;
}

/// <summary>Complete data read from a multiphysics binary export.</summary>
public class MultiphysicsBinaryData
{
    public MultiphysicsBinaryHeader Header { get; internal set; }

    /// <summary>Node positions (1D only, null for 2D).</summary>
    public float[] Positions { get; internal set; }

    /// <summary>
    /// All layers flattened: index = timeStep * layerCount + layerIndex.
    /// Each entry is a float[] of length (Nx*Ny) for 2D or Nx for 1D.
    /// </summary>
    public float[][] Layers { get; internal set; }
}
