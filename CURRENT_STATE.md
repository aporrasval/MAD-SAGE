# CURRENT_STATE.md

Current progress, active goals, and known blockers for the MAD_SAGE project.

---

## Current Goal

Validate the MAD black hole accretion module (`model_mad.c`) on IllustrisTNG100-1 merger trees and calibrate free parameters via autocalibration.

---

## Recent Changes

### [2026-05-06] Bug fix: RK45 integrator stall on stiff ODEs (`model_mad.c`)

**Problem:** The inner `while(!step_accepted)` loop in `integrate_spin_RK45()` had no minimum substep floor. On stiff ODEs (triggered by massive halos in TNG100 with extreme accretion conditions), `substep` would shrink toward zero. This caused k-values → 0, truncation error → 0, `stepsize_factor → ∞`, and a "step accepted" with `time += 0` — locking the outer loop indefinitely. The TNG100 run stalled for >3 hours instead of the expected 1–1.5 hours.

**Fix:** Added a substep floor of `dt_years * 1e-10` in `model_mad.c` (lines ~465–490). If `substep` falls below this threshold, the current proposed spin/mass state is force-accepted and integration continues.

**Impact on outputs:** Behavior change for galaxies whose spin ODE hits the floor. In the stalling case the outputs were undefined; now they reflect the best available proposed state at the minimum step. If the floor is hit repeatedly for the same galaxy, `MAX_SUBSTEPS = 500` may cause only a small fraction of the timestep to be integrated — those galaxies will show near-static spin evolution for that snapshot. Monitor for a population of galaxies with unexpectedly constant spin in TNG100 output.

**Files changed:** `code/model_mad.c`

---

### [2026-05-06] Bug fix: missing `RadiativeEfficiency` in `minimillennium.par`

**Problem:** `RadiativeEfficiency` was commented out in `input/minimillennium.par`, causing an assertion failure at startup. The parameter is still required by the parameter file reader even though it is unused at runtime when `SpinEvolutionOn=1` (both code paths that use it are guarded by `if(SpinEvolutionOn == 0)`).

**Fix:** Uncommented with value `0.10` and added an inline comment noting it is unused when `SpinEvolutionOn=1`.

**Files changed:** `input/minimillennium.par`

---

## Known Issues / Blockers

- **TNG100 stall**: Fix applied above. Needs validation — re-run TNG100 on Bouchet to confirm the stall is resolved and runtime returns to expected 1–1.5 hours.
- **Uninitialized `quasar_energy`** (`model_mergers.c:175`): Compiler warns this could be used uninitialized when `SpinEvolutionOn=0`. Harmless with current settings but should be fixed.
- **Unused variables** (`model_mad.c:584,842`): `r_Edd_final` and `underRadical` — leftover from development, no functional impact.

---

## Environment

- Local machine: macOS, Clang, GSL 2.6 at `/usr/local/Cellar/gsl/2.6`
- Supercomputer: Bouchet (Yale Center for Research Computing)
- Test simulation: Mini-Millennium (local, 8 tree files)
- Production simulation: IllustrisTNG100-1 (Bouchet)
