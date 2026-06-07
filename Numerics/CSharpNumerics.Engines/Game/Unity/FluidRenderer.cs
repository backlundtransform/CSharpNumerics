using CSharpNumerics.Engines.Game.Fluids;
using System;

namespace CSharpNumerics.Engines.Game.Unity;

/// <summary>
/// Extracts density and velocity data from a <see cref="GameFluidSolver3D"/> 
/// in formats suitable for Unity VFX Graph or shader consumption.
///
/// Provides flat float arrays that can be uploaded to Unity 3D textures:
///   - Density texture (single channel)
///   - Velocity texture (3 channels: R=vx, G=vy, B=vz)
///
/// No Unity dependency — produces raw arrays that the Unity shim uploads.
/// </summary>
public class FluidRenderer
{
    private readonly GameFluidSolver3D _solver;
    private float[] _densityTexture;
    private float[] _velocityTexture;

    /// <summary>Grid dimensions.</summary>
    public int NX => _solver.NX;

    /// <summary>Grid dimensions.</summary>
    public int NY => _solver.NY;

    /// <summary>Grid dimensions.</summary>
    public int NZ => _solver.NZ;

    /// <summary>Total voxel count.</summary>
    public int VoxelCount => NX * NY * NZ;

    /// <summary>Density scale factor for visualization.</summary>
    public double DensityScale { get; set; } = 1.0;

    /// <summary>Velocity scale factor for visualization.</summary>
    public double VelocityScale { get; set; } = 0.1;

    /// <summary>Minimum density threshold for rendering (below this = transparent).</summary>
    public double DensityThreshold { get; set; } = 0.01;

    /// <summary>
    /// Creates a fluid renderer for the given 3D solver.
    /// </summary>
    public FluidRenderer(GameFluidSolver3D solver)
    {
        _solver = solver ?? throw new ArgumentNullException(nameof(solver));
        int size = NX * NY * NZ;
        _densityTexture = new float[size];
        _velocityTexture = new float[size * 3];
    }

    /// <summary>
    /// Update the density texture data from the current solver state.
    /// Returns a float array suitable for uploading to a Unity Texture3D (R float).
    /// Values are clamped to [0,1].
    /// </summary>
    public float[] UpdateDensityTexture()
    {
        var density = _solver.Density;
        int size = VoxelCount;

        if (_densityTexture.Length != size)
            _densityTexture = new float[size];

        double scale = DensityScale;
        double threshold = DensityThreshold;

        for (int i = 0; i < size; i++)
        {
            double d = density[i] * scale;
            _densityTexture[i] = d < threshold ? 0f : (float)Math.Min(d, 1.0);
        }

        return _densityTexture;
    }

    /// <summary>
    /// Update the velocity texture data from the current solver state.
    /// Returns a float array of [vx,vy,vz, vx,vy,vz, ...] suitable for a Texture3D (RGB float).
    /// Values are scaled and clamped to [-1,1].
    /// </summary>
    public float[] UpdateVelocityTexture()
    {
        var vx = _solver.VelocityX;
        var vy = _solver.VelocityY;
        var vz = _solver.VelocityZ;
        int size = VoxelCount;

        if (_velocityTexture.Length != size * 3)
            _velocityTexture = new float[size * 3];

        double scale = VelocityScale;
        for (int i = 0; i < size; i++)
        {
            int idx = i * 3;
            _velocityTexture[idx] = ClampFloat(vx[i] * scale);
            _velocityTexture[idx + 1] = ClampFloat(vy[i] * scale);
            _velocityTexture[idx + 2] = ClampFloat(vz[i] * scale);
        }

        return _velocityTexture;
    }

    /// <summary>
    /// Extract a 2D density slice at a given Z index.
    /// Returns a float array of NX*NY for a 2D texture.
    /// </summary>
    public float[] GetDensitySlice(int zIndex)
    {
        if (zIndex < 0 || zIndex >= NZ)
            throw new ArgumentOutOfRangeException(nameof(zIndex));

        var density = _solver.Density;
        var slice = new float[NX * NY];
        int nxy = NX * NY;

        for (int j = 0; j < NY; j++)
            for (int i = 0; i < NX; i++)
                slice[i + j * NX] = (float)(density[i + j * NX + zIndex * nxy] * DensityScale);

        return slice;
    }

    /// <summary>
    /// Compute maximum density in the field (for auto-scaling).
    /// </summary>
    public double MaxDensity()
    {
        var density = _solver.Density;
        double max = 0;
        for (int i = 0; i < density.Length; i++)
            if (density[i] > max) max = density[i];
        return max;
    }

    /// <summary>
    /// Check if the renderer has valid data at the current solver state.
    /// Returns true if density texture has been populated and has non-zero values.
    /// </summary>
    public bool HasValidData()
    {
        if (_densityTexture == null || _densityTexture.Length == 0) return false;
        for (int i = 0; i < _densityTexture.Length; i++)
            if (_densityTexture[i] > 0) return true;
        return false;
    }

    private static float ClampFloat(double v)
    {
        if (v < -1.0) return -1f;
        if (v > 1.0) return 1f;
        return (float)v;
    }
}
