# CURRENT_STATE.md

Current progress, active goals, and known blockers for the MAD_SAGE project.

---

## Current Goal

Validate and diagnose the MAD black hole accretion module (`model_mad.c`) using individual BH history plots on the Mini-Millennium run. Current focus: understanding the decomposition of BH growth into radio (hot-gas) and instability-driven channels.

---

## Recent Changes

### [2026-06-13] New fields: RadioEddrate and InstaEddrate

Added supply-based Eddington ratio fields for each accretion channel, computed in the same way as `Eddrate` (= supply / mdot_Edd). Their sum equals `Eddrate` exactly.

**Why:** The previous approach split `Eddrate` proportionally using *net* accreted masses (`RadioBHaccreted`, `InstaBHaccreted`), which was an approximation. These new fields are derived directly from the supply fractions at the time `Eddrate` is set, giving an exact decomposition.

**Physics:** `RadioEddrate = Eddrate * radio_fraction`, `InstaEddrate = Eddrate * (1 - radio_fraction)`, where `radio_fraction = radio_supply / total_supply` computed in `evolve_black_hole()`.

**Files changed:**
- `code/core_allvars.h` — added `float RadioEddrate`, `float InstaEddrate` to all three struct definitions (GALAXY_OUTPUT small, large, and full GALAXY)
- `code/model_misc.c` — initialized both to 0.0 in `init_galaxy()`
- `code/model_mad.c` — computed both after `radio_fraction` is known, before the SpinEvolutionOn branch (~line 633)
- `code/core_save.c` — output both in both `prepare_galaxy_for_output` functions (~lines 313, 539)
- `output/routines.py` — added both to numpy dtype after `Eddrate`
- `output/plot_bh_history.py` — updated history builder to read fields directly instead of proportional-split approximation

### [2026-06-11] Bug fix: InstaBHaccreted always zero

**Problem:** `InstaBHaccreted` was written with `=` (overwrite) inside the sub-step loop. Instability accretion fires in early sub-steps; the last sub-step usually has none, so the field was always 0 at output time.

**Fix:** Reset both `InstaBHaccreted` and `RadioBHaccreted` to 0 before the sub-step loop in `core_build_model.c`, use `+=` inside the loop, and divide by `deltaT` after the loop to get snapshot-averaged rates. Same `+=` fix applied to `RadioBHaccreted` in `model_mad.c`.

**Files changed:** `code/core_build_model.c`, `code/model_mad.c`

### [2026-06-11] Par file fix: EtaFloor not in C code

`input/minimillennium_mad_mixed_etafloor_z0.par` had `EtaFloor 0.1` but the C code has no such parameter, causing an abort at startup. Commented out the line.

**File changed:** `input/minimillennium_mad_mixed_etafloor_z0.par`

### [2026-06-11] BH history plot: two-column comparison (quiescent vs active)

`output/plot_bh_history.py` now produces a 5-row × 2-column figure:
- Row 1: BH mass vs time
- Row 2: BH spin vs time
- Row 3: Total Eddington ratio (`Eddrate`)
- Row 4: Cumulative BH mass split by channel (RadioBlackHoleMass solid, InstaBlackHoleMass dashed)
- Row 5: Per-channel Eddington ratio (RadioEddrate solid, InstaEddrate dashed)

Left column: 10 randomly drawn quiescent BHs (Eddrate < 0.01 for ALL z < 4)
Right column: 10 randomly drawn active BHs (Eddrate > 0.01 for ALL z < 4)
Selection pool: log M_star > 9 AND log M_BH > 6 at z=0 (~1570 galaxies)
Output: `output/bh_history_comparison.pdf` and `.png`

### [2026-05-06] Bug fix: missing `RadiativeEfficiency` in `minimillennium.par`

**Fix:** Uncommented with value `0.10`; unused when `SpinEvolutionOn=1`.

**File changed:** `input/minimillennium.par`

---

## Output Directories

| Par file | Output dir | NumOutputs | Purpose |
|---|---|---|---|
| `input/minimillennium.par` | `output/history_results/` | -1 (all 64 snaps) | BH history plots |
| `input/minimillennium_mad_mixed_z0.par` | `output/mad_mixed_results/` | 1 (z=0 only) | Autocalibration |

**IMPORTANT:** Always run `/bin/rm -f output/history_results/model_z*` (or the relevant dir) before rerunning DarkSAGE, as stale files with the old struct layout will cause dtype mismatches in Python.

---

## Known Issues / Blockers

- **Uninitialized `quasar_energy`** (`model_mergers.c:175`): Compiler warns this could be used uninitialized when `SpinEvolutionOn=0`. Harmless with current settings but should be fixed.
- **Unused variables** (`model_mad.c`): `r_Edd_final` and `underRadical` — leftover from development, no functional impact.
- **`minimillennium.par` SpinEvolutionOn=1 only:** The history run par file always uses MAD spin physics. There is no all-snapshot par file for the `SpinEvolutionOn=0` path.

---

## Environment

- Local machine: macOS, Clang, GSL 2.6 at `/usr/local/Cellar/gsl/2.6`
- Supercomputer: Bouchet (Yale Center for Research Computing)
- Test simulation: Mini-Millennium (local, 8 tree files, 64 snapshots)
- Production simulation: IllustrisTNG100-1 (Bouchet)
- `rm` is aliased to `rm -i` — always use `/bin/rm -f` for non-interactive deletion
