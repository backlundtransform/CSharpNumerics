# Engine extraction — status and remaining steps

This document tracks the split of the simulation engines out of the open-source
`CSharpNumerics` core into a separate `CSharpNumerics.Engines` package intended
to live in a **private** repository and depend on the public `CSharpNumerics`
NuGet package.

## Rationale

`CSharpNumerics` is marketed as an open-source (MIT) scientific numerics library.
The simulation engines (Game physics, GIS, Multiphysics, Quantum, Audio,
Exoplanet) were built to sell Unity assets and as commercial productizations.
Splitting them keeps the core focused and open, and keeps the engine source
closed/proprietary. The sole copyright holder may relicense the engine code in
the private repo; MIT on the core permits the engine package to depend on it
without restriction.

## Done (in this branch / PR)

- **Step 0** — removed all reverse dependencies; core no longer references the
  engine layer.
- **Step 1a** — `CSharpNumerics.Engines` project created
  (`net10.0;net8.0;netstandard2.1`, `netstandard2.1` retained for Unity),
  engine source moved out of the core project. Core builds independently.
- **Step 1b** — tests split into `NumericTest` (core) and `NumericTest.Engines`.
- **Step 2** — core package metadata cleaned (engine-only tags removed:
  GamePhysics, Audio, Gis, Exoplanets; README engine section removed) and
  version bumped to **4.0.0** (breaking: anyone consuming engine types via the
  `CSharpNumerics` NuGet package must switch to `CSharpNumerics.Engines`).

Neither core nor engines has any external NuGet dependency, so the engine
package will depend only on `CSharpNumerics`.

## Remaining — gated behind review (do NOT run before the PR is approved)

These steps publish artifacts or create external repositories. They are
intentionally not part of the reviewed diff.

1. **Publish core `CSharpNumerics` 4.0.0 to NuGet.**
   ```
   dotnet pack Numerics/Numerics/CSharpNumerics.csproj -c Release
   dotnet nuget push <CSharpNumerics.4.0.0.nupkg> --source nuget.org --api-key <KEY>
   ```

2. **Create the private repository** (e.g. `backlundtransform/CSharpNumerics.Engines`)
   and move these projects into it:
   - `CSharpNumerics.Engines/`
   - `NumericTest.Engines/`
   - a solution + the `LICENSE` of your choice (may be proprietary).

3. **Swap the engine project's reference** from project to package:
   ```xml
   <!-- remove -->
   <ProjectReference Include="..\Numerics\CSharpNumerics.csproj" />
   <!-- add -->
   <PackageReference Include="CSharpNumerics" Version="4.0.0" />
   ```
   (This cannot build green until step 1 has published the package, which is why
   it is not done in this PR.)

4. **Remove the engine projects from the public repo** (follow-up PR): delete
   `CSharpNumerics.Engines/` and `NumericTest.Engines/` and drop them from
   `Numerics.sln`. Until this is done the engine source remains public in this
   repo; the split so far is structural only.

## Known pre-existing test failure (unrelated to this work)

`NumericTest.Engines/CollisionTests.FullScenario_TwoSpheres_BounceApart` fails on
`master` as well (the two spheres are positioned 3 apart with radii 1 each, so
they do not overlap, yet the test asserts a contact). Left untouched here.
