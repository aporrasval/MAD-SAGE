# DECISIONS.md

Key design and implementation decisions with rationale. Append new entries chronologically; do not delete old ones.

---

## [2026-05-06] RK45 substep floor: force-accept at `dt_years * 1e-10`

**Context:** The `integrate_spin_RK45()` function in `model_mad.c` uses an adaptive RK45 integrator to evolve BH spin and mass over a SAM timestep. The inner `while(!step_accepted)` loop shrinks `substep` whenever the truncation error exceeds tolerance. No minimum floor existed, so on stiff ODEs the loop could run indefinitely — confirmed as the cause of a >3 hour stall when running on IllustrisTNG100-1 trees on the Bouchet supercomputer.

**Decision:** Add a hard substep floor of `dt_years * 1e-10`. If `substep` falls below this threshold after the error-controller adjustment, force-accept the current proposed spin/mass state and advance time by the floor amount.

**Alternatives considered:**

1. *Maximum inner iteration count* — cap the number of rejected steps and break. Rejected because it leaves the state undefined (break on a rejected step means proposed values may be far from converged), and the resulting spin/mass would be unreliable.

2. *Force-accept at floor (chosen)* — physically motivated: if the adaptive controller drives `substep` to near machine precision, the ODE is locally saturated (spin is near a limit, or the system is quasi-stationary). Accepting the proposed state at the floor is equivalent to taking the best available estimate at the smallest meaningful timescale.

**Known limitation:** If the floor is triggered repeatedly, each forced step advances time by only `dt_years * 1e-10`. With `MAX_SUBSTEPS = 500`, the integrator may cover at most `500 * dt_years * 1e-10 = dt_years * 5e-8` before the outer loop exits. Affected galaxies will show near-static spin for that timestep. This is conservative and detectable in post-analysis — look for a population of galaxies with constant spin across a snapshot in TNG100 output.

**Parameter choice — `1e-10`:** Small enough to be below any physically meaningful timescale (even at `dt_years ~ 500 Myr`, the floor is `~50 years`), large enough to avoid true floating-point underflow (`double` underflows below ~`1e-308`). The value can be tightened or loosened without affecting normal runs (the floor is never reached under non-stiff conditions).

**Files changed:** `code/model_mad.c` — `integrate_spin_RK45()`, inner `while(!step_accepted)` loop.
