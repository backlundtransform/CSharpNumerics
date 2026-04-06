namespace CSharpNumerics.Physics.Optics;

/// <summary>
/// Common interface for any surface that can interact with light rays.
/// </summary>
public interface IOpticalSurface
{
    /// <summary>
    /// Tests whether the given ray intersects this surface.
    /// Returns null when there is no intersection.
    /// </summary>
    /// <param name="ray">The incident ray.</param>
    /// <param name="tMin">Minimum ray parameter (avoids self-intersection).</param>
    /// <param name="tMax">Maximum ray parameter.</param>
    RayHit? Intersect(Ray ray, double tMin = 1e-6, double tMax = double.MaxValue);
}
