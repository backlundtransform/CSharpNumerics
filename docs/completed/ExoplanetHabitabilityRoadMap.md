# Exoplanet Habitability Roadmap

Formulas needed to derive exoplanet catalogue fields (as in the PHL Habitable Worlds Catalog)
from NASA Exoplanet Archive data, so ExoplanetHunter no longer depends on the PHL catalogue.
All results were cross-checked against PHL values.

## Phase 1 — Research & Design
- [x] Identify missing formulas (luminosity, insolation, equilibrium temperature, flux ESI, optimistic HZ)
- [x] Verify Kopparapu et al. (2014) coefficients against the official `HZ_coefficients.dat`
- [x] Confirm PHL uses the 2014 optimistic limits (Kepler-1556 reproduced to 1e-7 AU)

## Phase 2 — Core Implementation
- [x] `HabitableZoneLimit` and `HabitableZoneModel` enums
- [x] `PhysicsConstants.SolarEffectiveTemperature` (IAU 2015, 5772 K)
- [x] `StellarLuminosity`, `InsolationFlux`, `InsolationFluxRange`, `EquilibriumTemperature`
- [x] `CalculateEsiFromFlux`
- [x] `HabitableZoneFlux` and `CalculateGoldilocksZone(L, Teff, HabitableZoneModel)`
- [x] Unit tests in `AstronomyTests`
- [x] Tests pass in CI (local test runs are blocked by Windows App Control on the dev machine)

## Phase 3 — Docs
- [x] Update Physics README (Exoplanet Classification section)

## Notes
- The existing two-argument `CalculateGoldilocksZone(L, Teff)` is unchanged. Its inner-edge
  coefficients (S_eff,☉ = 1.0146) correspond to the Kopparapu et al. (2013) *moist* greenhouse
  limit, although the code comment calls it runaway greenhouse. Aligning it with the 2014
  conservative model would change results and is left for a decision.
