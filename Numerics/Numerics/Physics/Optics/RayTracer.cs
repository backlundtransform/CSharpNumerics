using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Physics.Optics;

/// <summary>
/// Result of tracing a single ray through an optical scene.
/// Contains the full path of the ray including all bounces and refractions.
/// </summary>
public class TraceResult
{
    /// <summary>The original ray.</summary>
    public Ray OriginalRay { get; }

    /// <summary>Sequential segments of the traced path.</summary>
    public List<TraceSegment> Segments { get; } = new();

    /// <summary>Whether the ray was terminated early (exceeded max bounces or intensity cutoff).</summary>
    public bool Terminated { get; internal set; }

    public TraceResult(Ray originalRay) => OriginalRay = originalRay;
}

/// <summary>
/// A single segment of a traced ray path: from one intersection to the next.
/// </summary>
public readonly struct TraceSegment
{
    /// <summary>Start point of this segment.</summary>
    public Vector Start { get; }

    /// <summary>End point (intersection with a surface).</summary>
    public Vector End { get; }

    /// <summary>Ray intensity at the start of this segment.</summary>
    public double Intensity { get; }

    /// <summary>The medium through which this segment travels.</summary>
    public OpticalMedium Medium { get; }

    /// <summary>The surface hit at the end (null for the final segment that goes to infinity).</summary>
    public IOpticalSurface HitSurface { get; }

    public TraceSegment(Vector start, Vector end, double intensity,
        OpticalMedium medium, IOpticalSurface hitSurface)
    {
        Start = start;
        End = end;
        Intensity = intensity;
        Medium = medium;
        HitSurface = hitSurface;
    }
}

/// <summary>
/// A basic recursive ray tracer for geometric optics.
/// Traces rays through an <see cref="OpticalScene"/>, handling reflection,
/// refraction, and absorption at each surface interaction.
/// </summary>
public class RayTracer
{
    /// <summary>Maximum number of reflections/refractions (default 16).</summary>
    public int MaxBounces { get; set; } = 16;

    /// <summary>Rays with intensity below this threshold are terminated (default 0.001).</summary>
    public double IntensityCutoff { get; set; } = 0.001;

    /// <summary>The optical scene to trace through.</summary>
    public OpticalScene Scene { get; }

    public RayTracer(OpticalScene scene)
    {
        Scene = scene;
    }

    /// <summary>
    /// Traces a single ray through the scene, returning all reflected/refracted paths.
    /// </summary>
    public TraceResult Trace(Ray ray)
    {
        var result = new TraceResult(ray);
        TraceRecursive(ray, Scene.AmbientMedium, 0, result);
        return result;
    }

    /// <summary>
    /// Traces multiple rays, returning all results.
    /// </summary>
    public List<TraceResult> TraceAll(IEnumerable<Ray> rays)
    {
        var results = new List<TraceResult>();
        foreach (var ray in rays)
            results.Add(Trace(ray));
        return results;
    }

    private void TraceRecursive(Ray ray, OpticalMedium currentMedium,
        int depth, TraceResult result)
    {
        if (depth >= MaxBounces || ray.Intensity < IntensityCutoff)
        {
            result.Terminated = true;
            return;
        }

        var hit = Scene.ClosestHit(ray);
        if (!hit.HasValue) return; // ray escapes

        RayHit h = hit.Value;

        // Beer–Lambert absorption through the current medium
        double intensity = OpticsExtensions.AttenuateIntensity(
            ray.Intensity, currentMedium.AbsorptionCoefficient, h.Distance);

        result.Segments.Add(new TraceSegment(
            ray.Origin, h.Point, intensity, currentMedium, h.Surface));

        // Determine interaction based on surface type
        if (h.Surface is PlaneMirror mirror)
        {
            // Pure reflection
            Vector reflected = OpticsExtensions.Reflect(ray.Direction, h.Normal);
            var newRay = new Ray(h.Point, reflected, ray.WavelengthNm,
                intensity * mirror.Reflectivity);
            TraceRecursive(newRay, currentMedium, depth + 1, result);
        }
        else if (h.Surface is SphericalMirror sMirror)
        {
            // Pure reflection
            Vector reflected = OpticsExtensions.Reflect(ray.Direction, h.Normal);
            var newRay = new Ray(h.Point, reflected, ray.WavelengthNm,
                intensity * sMirror.Reflectivity);
            TraceRecursive(newRay, currentMedium, depth + 1, result);
        }
        else if (h.Surface is ThinLens lens)
        {
            // Thin lens refracts using paraxial model
            var refractedRay = lens.RefractRay(
                new Ray(ray.Origin, ray.Direction, ray.WavelengthNm, intensity));
            TraceRecursive(refractedRay, currentMedium, depth + 1, result);
        }
        else if (h.Surface is CircularAperture || h.Surface is RectangularAperture)
        {
            // Aperture: ray passes through unchanged
            var continuedRay = new Ray(h.Point, ray.Direction, ray.WavelengthNm, intensity);
            TraceRecursive(continuedRay, currentMedium, depth + 1, result);
        }
        else if (h.Surface is ImageSensor)
        {
            // Sensor absorbs the ray — no further propagation
        }
        else
        {
            // Generic refracting surface: compute Fresnel split
            double n1 = currentMedium.RefractiveIndex;
            double n2 = h.MediumBeyond.RefractiveIndex;
            if (n2 <= 0) n2 = n1; // fallback

            double cosI = Math.Abs(ray.Direction.Dot(h.Normal));
            double R = OpticsExtensions.FresnelReflectance(n1, n2, Math.Acos(cosI));

            // Reflected ray
            if (R > IntensityCutoff)
            {
                Vector reflected = OpticsExtensions.Reflect(ray.Direction, h.Normal);
                var reflRay = new Ray(h.Point, reflected, ray.WavelengthNm, intensity * R);
                TraceRecursive(reflRay, currentMedium, depth + 1, result);
            }

            // Refracted ray
            double T = 1.0 - R;
            if (T > IntensityCutoff)
            {
                var refracted = OpticsExtensions.Refract(ray.Direction, h.Normal, n1, n2);
                if (refracted.HasValue)
                {
                    var refrRay = new Ray(h.Point, refracted.Value,
                        ray.WavelengthNm, intensity * T);
                    TraceRecursive(refrRay, h.MediumBeyond, depth + 1, result);
                }
                // else TIR — all energy already in reflected ray
            }
        }
    }
}
