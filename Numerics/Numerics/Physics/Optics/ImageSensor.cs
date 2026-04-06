using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Physics.Optics;

/// <summary>
/// A rectangular image sensor (film plane / CCD array) that collects ray hits.
/// Each hit records the position on the sensor and the ray properties.
/// </summary>
public class ImageSensor : IOpticalSurface
{
    /// <summary>Centre of the sensor.</summary>
    public Vector Center { get; }

    /// <summary>Normal to the sensor plane.</summary>
    public Vector Normal { get; }

    /// <summary>Local "right" direction on the sensor plane.</summary>
    public Vector Right { get; }

    /// <summary>Local "up" direction on the sensor plane.</summary>
    public Vector Up { get; }

    /// <summary>Half-width of the sensor.</summary>
    public double HalfWidth { get; }

    /// <summary>Half-height of the sensor.</summary>
    public double HalfHeight { get; }

    /// <summary>Horizontal resolution (number of pixels).</summary>
    public int ResolutionX { get; }

    /// <summary>Vertical resolution (number of pixels).</summary>
    public int ResolutionY { get; }

    private readonly List<SensorHit> _hits = new();
    private readonly double[,] _accumulator;

    /// <summary>All recorded hits.</summary>
    public IReadOnlyList<SensorHit> Hits => _hits;

    /// <summary>
    /// Creates a rectangular sensor.
    /// </summary>
    public ImageSensor(Vector center, Vector normal, Vector right,
        double width, double height, int resolutionX = 256, int resolutionY = 256)
    {
        Center = center;
        Normal = normal.GetUnitVector();
        Right = right.GetUnitVector();
        Up = Normal.Cross(Right).GetUnitVector();
        HalfWidth = width / 2.0;
        HalfHeight = height / 2.0;
        ResolutionX = resolutionX;
        ResolutionY = resolutionY;
        _accumulator = new double[resolutionX, resolutionY];
    }

    /// <inheritdoc />
    public RayHit? Intersect(Ray ray, double tMin = 1e-6, double tMax = double.MaxValue)
    {
        double denom = ray.Direction.Dot(Normal);
        if (Math.Abs(denom) < 1e-12) return null;

        double t = (Center - ray.Origin).Dot(Normal) / denom;
        if (t < tMin || t > tMax) return null;

        Vector hitPoint = ray.PointAt(t);
        Vector offset = hitPoint - Center;
        double u = offset.Dot(Right);
        double v = offset.Dot(Up);

        if (Math.Abs(u) > HalfWidth || Math.Abs(v) > HalfHeight) return null;

        // Record hit
        int px = (int)((u / HalfWidth * 0.5 + 0.5) * (ResolutionX - 1));
        int py = (int)((v / HalfHeight * 0.5 + 0.5) * (ResolutionY - 1));
        px = Math.Max(0, Math.Min(ResolutionX - 1, px));
        py = Math.Max(0, Math.Min(ResolutionY - 1, py));

        _hits.Add(new SensorHit(hitPoint, u, v, px, py, ray.Intensity, ray.WavelengthNm));
        _accumulator[px, py] += ray.Intensity;

        Vector outNormal = denom < 0 ? Normal : (-1.0) * Normal;
        return new RayHit(hitPoint, outNormal, t,
            new OpticalMedium("Sensor", 0), this);
    }

    /// <summary>Returns the accumulated intensity at pixel (x, y).</summary>
    public double GetPixelIntensity(int x, int y) => _accumulator[x, y];

    /// <summary>Returns the full intensity image as a 2D array [x, y].</summary>
    public double[,] GetImage() => (double[,])_accumulator.Clone();

    /// <summary>Clears all recorded hits and resets the accumulator.</summary>
    public void Clear()
    {
        _hits.Clear();
        Array.Clear(_accumulator, 0, _accumulator.Length);
    }

    /// <summary>Total number of hits recorded.</summary>
    public int HitCount => _hits.Count;
}

/// <summary>
/// A single hit on the image sensor.
/// </summary>
public readonly struct SensorHit
{
    /// <summary>World-space hit point.</summary>
    public Vector Point { get; }

    /// <summary>Local U coordinate on the sensor.</summary>
    public double U { get; }

    /// <summary>Local V coordinate on the sensor.</summary>
    public double V { get; }

    /// <summary>Pixel X index.</summary>
    public int PixelX { get; }

    /// <summary>Pixel Y index.</summary>
    public int PixelY { get; }

    /// <summary>Ray intensity at this hit.</summary>
    public double Intensity { get; }

    /// <summary>Wavelength of the ray in nm.</summary>
    public double WavelengthNm { get; }

    public SensorHit(Vector point, double u, double v, int pixelX, int pixelY,
        double intensity, double wavelengthNm)
    {
        Point = point;
        U = u;
        V = v;
        PixelX = pixelX;
        PixelY = pixelY;
        Intensity = intensity;
        WavelengthNm = wavelengthNm;
    }
}
