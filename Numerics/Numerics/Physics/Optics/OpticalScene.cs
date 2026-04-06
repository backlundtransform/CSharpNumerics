using System.Collections.Generic;

namespace CSharpNumerics.Physics.Optics;

/// <summary>
/// An optical scene containing surfaces and a current medium.
/// Surfaces are tested in insertion order; the closest intersection wins.
/// </summary>
public class OpticalScene
{
    private readonly List<IOpticalSurface> _surfaces = new();
    private readonly List<OpticalMedium> _media = new();

    /// <summary>
    /// The ambient medium that fills the scene (default: air).
    /// </summary>
    public OpticalMedium AmbientMedium { get; set; } = OpticalMaterialLibrary.Air;

    /// <summary>All surfaces in the scene.</summary>
    public IReadOnlyList<IOpticalSurface> Surfaces => _surfaces;

    /// <summary>
    /// Adds a surface to the scene.
    /// </summary>
    public void Add(IOpticalSurface surface) => _surfaces.Add(surface);

    /// <summary>
    /// Adds a surface paired with the medium on its far side.
    /// Used for refracting surfaces (e.g. glass/air boundaries).
    /// </summary>
    public void Add(IOpticalSurface surface, OpticalMedium mediumBeyond)
    {
        _surfaces.Add(surface);
        _media.Add(mediumBeyond);
    }

    /// <summary>
    /// Finds the closest intersection of the ray with any surface in the scene.
    /// Returns null when no surface is hit.
    /// </summary>
    public RayHit? ClosestHit(Ray ray, double tMin = 1e-6, double tMax = double.MaxValue)
    {
        RayHit? closest = null;
        double bestT = tMax;

        for (int i = 0; i < _surfaces.Count; i++)
        {
            var hit = _surfaces[i].Intersect(ray, tMin, bestT);
            if (hit.HasValue && hit.Value.Distance < bestT)
            {
                closest = hit;
                bestT = hit.Value.Distance;
            }
        }

        return closest;
    }
}
