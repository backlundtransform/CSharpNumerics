# Optics Roadmap

## Phase 1 — Core Types & Ray Primitives
- [x] `Ray` struct (origin, direction, wavelength, intensity)
- [x] `OpticalMedium` readonly struct (refractive index, absorption coefficient)
- [x] `OpticalMaterialLibrary` static class with common media (air, glass, water, diamond, etc.)
- [x] `RayHit` struct (hit point, normal, distance, medium)
- [x] `IOpticalSurface` interface

## Phase 2 — Reflection & Refraction (Snell's Law / Fresnel)
- [x] `OpticsExtensions` — Snell's law (refraction angle)
- [x] `OpticsExtensions` — Fresnel reflectance (s-pol, p-pol, unpolarized)
- [x] `OpticsExtensions` — Total internal reflection detection
- [x] `OpticsExtensions` — Beer–Lambert absorption

## Phase 3 — Optical Elements (Mirrors, Lenses, Prisms)
- [x] `PlaneMirror` — flat reflection surface
- [x] `SphericalMirror` — concave/convex with focal length
- [x] `ThinLens` — converging/diverging with lensmaker's equation
- [x] `Prism` — triangular prism with dispersion
- [x] `MirrorType` / `LensType` enums

## Phase 4 — Apertures & Sensors
- [x] `CircularAperture` — blocks rays outside radius
- [x] `RectangularAperture` — rectangular opening
- [x] `ImageSensor` — collects hits on a rectangular sensor plane

## Phase 5 — Basic Ray Tracer
- [x] `OpticalScene` — holds surfaces, media, light sources
- [x] `RayTracer` — traces ray through scene (recursive reflection/refraction)
- [x] Max bounce depth & intensity cutoff

## Phase 6 — Materials Integration
- [x] Add `OpticalMedium` to `Materials/Optical/` folder
- [x] Add factory method `Materials.Optical("CrownGlass")`

## Phase 7 — Game Engine Integration
- [x] `RaycastExtensions` for game-engine ray intersection with AABB/BoundingSphere

## Phase 8 — Tests & Documentation
- [x] Unit tests for all optics classes (49 tests)
- [x] Update Physics README
