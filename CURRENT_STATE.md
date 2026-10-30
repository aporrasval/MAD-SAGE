# CURRENT_STATE.md

Current progress, active goals, and known blockers for the MAD_SAGE project.

---

## Current Goal

Calibrate MAD-SAGE on Mini-Millennium with the PSO autocalibration (`autocalibration/`), then push the current version to GitHub (`origin` = github.com/aporrasval/MAD-SAGE).

Previous goal (June 2026): validate the MAD module with BH history plots and the radio/instability decomposition of BH growth.

---

## Recent Changes

### [2026-10-30] Renamed darksage → madsage
Executable is now `madsage` (Makefile `EXEC`). `build_coherent.sh` deleted (it compiled from a missing `coherentonlylcode/` folder); run `./madsage input/minimillennium.par` and pass `-b ../madsage` to the autocalibration. Python renames: `execution.run_madsage`, `common.read_madsage_par`, `routines.madsage_snap` / `madsage_out_single` / `galdtype_madsage` (old `darksage_*` names kept as aliases in `routines.py`). Per-particle autocalibration output folder renamed `DS_output_<pid>` → `MS_output_<pid>`. `main.py -b` (`--madsage-binary`, alias `--shark-binary`) now defaults to `../madsage` instead of Shark's `shark`/`../build/shark`. `-S` in `main.py`, `diagnostics.py`, `post_analysis.py` and `best_fit_params.py` now defaults to `autocalibration/space_BHspin.txt` (was Shark's `space.txt`). PSO log renamed `shark_pso.log` → `madsage_pso.log`. `SMF_z1` constraint removed (missing data file). `input/minimillennium.par` now uses repo-relative paths (`./output/results`, `./mini_millennium/treedata`, `./mini_millennium/millennium.a_list`); `constraints.configure()` resolves relative paths against the repo root. `mini_millennium/` (152 MB trees) is git-ignored, not pushed. LICENCE: added MAD-SAGE copyright and a note on the GPL/LGPL Shark-derived autocalibration files. Docs, comments and plot labels say MAD-SAGE; references to upstream DARK SAGE (papers, licence) unchanged. `autocalibration/README.md` added. **Refactor only** — no physics or output change. Entries below keep the old names as written at the time.

### [2026-10-29] Autocalibration fixed and tested on Mini-Millennium

**Classification:** bug fixes + refactor in `autocalibration/` only. No DarkSAGE physics changed. PSO outcomes differ from earlier runs (student-t fix, weights, log sampling); chi2 values themselves are unchanged.

- `constraints.configure(par_file)` reads `Hubble_h`, `Omega`, `BoxLen`, `FirstFile`/`LastFile`, `FileWithSnapList`, `FileNameGalaxies` from the `-c` .par file (replaces hard-coded `sim` switch). Volume = (BoxLen/h)^3 × N_files_run / N_files_total (total counted as `SimulationDir/TreeName.N`).
- `analysis.studentT` now returns a positive negative-log-likelihood (was +log-likelihood → every particle NaN); uses err², `gammaln`; Gaussian limit when var ≤ 1.
- Objective in `execution.run_darksage` = ∏ stat_i^w_i (weights were a no-op constant factor before). Non-finite/≤0 statistic → inf + warning.
- `IsLog=1` parameters now sampled in log10 (was linear). Tracks and log stay in physical units (`pso.py` `dump_transform`).
- `--seed` option (logged; auto-drawn if omitted) → runs reproducible (verified: same seed gives identical fitness).
- Temp .par writing: exact key match, whole file scanned, error if a parameter or `OutputDir` is missing. Stale-file cleanup fixed; temp .par deleted unless `-k`.
- Spawned workers (macOS) configure constraints and logging themselves, so `Particle ... evaluated to` lines reach stdout for `post_analysis.py`.
- `post_analysis.py` `-o` defaults to the dump file's directory.
- Directory cleanup: removed old run outputs, temp .par files, unused obs files (GAMA_SMF.dat, GAMAII_BBD_GSMFs.dat, HIMF_Jones18.dat, JWST_QLFz6.csv, Shen20_QLFz6.csv, weibel24_SMFz6.csv), `space_GUREFT.txt`, `space_fiducial_z0.txt`, `space_mad_etafloor_z0.txt`. Only search space now: `space_BHspin.txt` (GasSinkRate, EnergySN; both log).
- `.gitignore`: `output/results*/`, `autocalibration/DS_output_*/`, `autocalibration/*_temp.par`.

**Test:** `python main.py -c ../input/minimillennium.par -b ../darksage -S space_BHspin.txt -s 2 -m 2 --seed 12345 -t student-t -x "SMF_z0,HIMF,CSFRDH"` → ran in 460 s, finite fitness, tracks/plots OK.

### [2026-10-29] Mini-Millennium baseline run (`input/minimillennium.par`, z=0)

- Output: `output/results/`, 36,043 galaxies; deterministic (two runs identical). Record size 13,440 bytes/galaxy; `output/routines.py` (user-supplied version, adds `MergerBlackHoleMass`, `CurrentRadioBlackHoleMass`) matches the C struct.
- Params: GasSinkRate 0.71393, EnergySN 8.36414e50 (from old PSO run 05, older code).
- Fit vs observations (reduced chi2, constraints.py definitions): SMF_z0 (8–13) 101; SMF_z0(9–12) 131; HIMF 206; CSFRDH 16.8.
- Physics flags: SMF too low at log M* 10.2–10.9, too high above 11.1 (up to +1.4 dex at 11.6). BHs barely grow in M* < 1e11 galaxies (median log M_BH ≈ 2.8 vs 137 Msun seed; 65% below 1e4 Msun at log M* 9–10).

### [2026-10-29] Uncommitted C changes in working tree (status unresolved)

Relative to commit `8b0fd50`, the working tree removes `RadioEddrate`/`InstaEddrate` (June 13 entry below), reverts `InstaBHaccreted`/`RadioBHaccreted` to per-sub-step `= x/dt` (undoes June 11 fix), stops writing Primary/Secondary BH mass/spin in the large output (always 0), changes `eta_floor` 0.1 → 0.02 in `bh_feedback()` (physics), and converts BHaccreted to Msun/yr before computing `Lbol`. The June 11/13 entries below describe the committed state, not the working tree.

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
| `input/minimillennium.par` | `output/results/` | 1 (z=0 only) | Baseline run; template for autocalibration (`-c`) |

Only `input/minimillennium.par` remains in `input/` (other .par files were removed in the Oct 2026 cleanup). Autocalibration writes per-particle output to `<run dir>/MS_output_<pid>/` (e.g. `autocalibration/minimillennium/`) and deletes it after each evaluation.

**IMPORTANT:** Always run `/bin/rm -f output/results/model_z*` (or the relevant dir) before rerunning MAD-SAGE, as stale files with the old struct layout will cause dtype mismatches in Python.

---

## Known Issues / Blockers

- **Uncommitted C changes need a keep/revert decision** before pushing to GitHub (see 2026-10-29 entry above).
- **Autocalibration leftovers:** `-H` (HPC) and `-v` (subvolumes) are Shark code that do nothing for MAD-SAGE; `chi2` hard-codes 2 free parameters in its reduced-chi2 denominator; parallel workers can interleave stdout in `full_dump.txt`; worker log timestamps are local time vs UTC in `main.py`.
- **TNG50-4 volume:** old constraints.py used the full 35 Mpc/h box with 1 tree file; `configure()` now scales by files run / files total — check on the cluster.

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
