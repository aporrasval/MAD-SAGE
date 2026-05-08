# CURRENT_STATE.md

Current progress, active goals, and known blockers for the MAD_SAGE project.

---

## Current Goal

Validate the MAD black hole accretion module (`model_mad.c`) on IllustrisTNG100-1 merger trees and calibrate free parameters via autocalibration.

---

## Recent Changes

### [2026-05-06] Bug fix: missing `RadiativeEfficiency` in `minimillennium.par`

**Problem:** `RadiativeEfficiency` was commented out in `input/minimillennium.par`, causing an assertion failure at startup. The parameter is still required by the parameter file reader even though it is unused at runtime when `SpinEvolutionOn=1` (both code paths that use it are guarded by `if(SpinEvolutionOn == 0)`).

**Fix:** Uncommented with value `0.10` and added an inline comment noting it is unused when `SpinEvolutionOn=1`.

**Files changed:** `input/minimillennium.par`

---

## Known Issues / Blockers

- **Uninitialized `quasar_energy`** (`model_mergers.c:175`): Compiler warns this could be used uninitialized when `SpinEvolutionOn=0`. Harmless with current settings but should be fixed.
- **Unused variables** (`model_mad.c`): `r_Edd_final` and `underRadical` — leftover from development, no functional impact.

---

## Environment

- Local machine: macOS, Clang, GSL 2.6 at `/usr/local/Cellar/gsl/2.6`
- Supercomputer: Bouchet (Yale Center for Research Computing)
- Test simulation: Mini-Millennium (local, 8 tree files)
- Production simulation: IllustrisTNG100-1 (Bouchet)
