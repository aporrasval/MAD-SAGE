# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Guidelines

### Code Style
- Use Python 3
- Prioritize modularity, readability, and editability
- Do not remove existing comments
- Add clear, verbose comments when modifying code
- Ask before introducing new dependencies
- Ask if Cython or C might be more efficient for expensive components

### Physics Integrity (Critical)
- Never modify, simplify, or refactor physics logic without explicit approval
- If a change could affect scientific results, stop and ask
- Clearly distinguish: **refactor** (no behavior change) vs. **behavior change** (affects results); when uncertain, assume behavior change
- Always surface hidden assumptions, implicit constants, and unit conventions
- When editing SAM logic, explain what physical process is affected and whether outputs will change

### Reproducibility
- All parameter choices must be traceable — never introduce magic numbers without explanation
- When modifying parameters, state the previous value, new value, and reason
- Ensure outputs can be tied to a parameter set, code version, and random seed (if applicable)

### Change Classification
Before making changes, classify as one of:
1. **Refactor** — no change in outputs
2. **Performance optimization**
3. **Bug fix**
4. **Scientific/behavior change**

For (3) or (4): explain expected impact on results and ask for confirmation before proceeding.

### Interaction Style
- Explain reasoning briefly before major changes
- Point out potential bugs or edge cases proactively
- Ask for clarification when ambiguity exists

---

## Project Memory Protocol

### Files to Maintain
- `CURRENT_STATE.md` — current progress, goals, blockers
- `ARCHITECTURE.md` — system design and invariants
- `DECISIONS.md` — key decisions and rationale
- `ROADMAP.md` — upcoming tasks and milestones

`CLAUDE.md` contains stable rules and protocols only. `ARCHITECTURE.md` contains evolving system design — do not duplicate content between them. If architecture changes, update `ARCHITECTURE.md`; only update `CLAUDE.md` if rules change.

### Session Start Protocol
At the start of every session:
1. Read `CLAUDE.md` and `CURRENT_STATE.md`
2. Summarize: current goal, key constraints, recent decisions
3. Do not proceed with work until this is done

### Required Workflow
At the end of any significant step:
1. Summarize what changed
2. Identify whether physics, parameters, or architecture were affected
3. Propose updates to the relevant memory files

After any non-trivial task, propose memory file updates and ask for confirmation before overwriting. If the session may end or context may reset, perform a full memory checkpoint.

### Experiment Tracking
When running parameter studies, record:
- Parameter set used
- Purpose of the run
- Output location

Encourage use of a structured log or file for this.

### Rules
- Keep entries concise and structured
- Do not duplicate information across files

---

## Project Overview

DARK SAGE is a semi-analytic model of galaxy formation written in C, simulating disc structure and evolution of galaxies. This fork (`Angelo_MAD`) adds the MAD (Magnetically-Arrested Disc) black hole accretion physics module (`model_mad.c`). Reference: Stevens et al. (2024) PASA 41:53.

## Build and Run

**Requirements:** GCC (Linux) or Clang (macOS, auto-detected), GSL (GNU Scientific Library)

```bash
# Build
make

# Alternative: build via Python wrapper (allows parameter file customization)
python compile.py

# Run
./madsage input/minimillennium.par

# Run with MPI (experimental)
mpirun -np 8 ./madsage input/minimillennium.par
```

There are no automated tests. Manual validation is done by running on a small simulation (Mini-Millennium) and checking galaxy property distributions against observations.

## Code Architecture

The model reads N-body merger trees, then evolves galaxies forward through simulation snapshots by applying physics modules in sequence. Disk structure is tracked on 30 radial annuli.

**Entry point:** `main.c` → `core_build_model.c` (main evolution loop per merger tree)

**Key files:**
- `core_allvars.c/.h` — Global variables and the central `GALAXY` struct (all galaxy properties)
- `core_proto.h` — All function prototypes
- `core_read_parameter_file.c` — Parses `.par` configuration files
- `core_io_tree.c` — Reads halo merger tree data
- `core_save.c` — Writes binary galaxy output
- `model_mad.c` — **MAD black hole physics** (this fork's main contribution): spin evolution via RK45 integration, Eddington luminosity, radiative efficiency, jet power, magnetization dynamics, merger spin calculations (Rezzolla et al. 2008)

**Physics modules** (applied per snapshot per galaxy):
- `model_cooling.c` — Gas cooling
- `model_infall.c` — Gas infall into halos
- `model_starformation_and_feedback.c` — Star formation and SN feedback
- `model_disk_instability.c` — Bar/instability physics
- `model_mergers.c` — Galaxy mergers
- `model_reincorporation.c` — Hot gas reincorporation
- `model_misc.c` — Miscellaneous helpers

## Parameter Files

The `.par` files in `input/` control all model behavior. Key sections:
- Output paths and file naming
- Simulation tree files and cosmology
- Physics toggles: AGN recipe, disk instability, ram pressure, H2 prescription
- Black hole parameters: spin evolution mode, fEdd threshold, seed mass
- Free parameters: `GasSinkRate`, `EnergySN`, `EpsilonR`, `EpsilonQ`

## Output

Binary galaxy catalogs. Use `output/routines.py` to read output — it defines the full numpy dtype matching the `GALAXY` struct in `core_allvars.h`. Galaxy properties include masses, gas reservoirs, disk annulus data (HI, H2, SFR per ring), black hole spin/mass/magnetization/jet power, and metallicity.

## MAD Module Details (`model_mad.c`)

The central contribution of this fork. Key physics implemented:
- **Spin evolution**: RK45 numerical integration of spin parameter `a` during accretion
- **Accretion modes**: coherent (prograde/retrograde) vs. chaotic, super- vs. sub-Eddington
- **ISCO radius**: spin-dependent innermost stable circular orbit
- **Radiative efficiency**: `η(a)` from spin-dependent ISCO
- **Jet power**: magnetically-driven jets from MAD state
- **Magnetization**: tracks magnetic flux threading the black hole
- **Merger spin**: post-merger spin using Rezzolla et al. (2008) fitting formulas

When modifying `model_mad.c`, the `GALAXY` struct fields for black holes are defined in `core_allvars.h` and must stay in sync with the output dtype in `output/routines.py`.

## Autocalibration

`autocalibration/` contains a Particle Swarm Optimization (PSO) system to tune MAD-SAGE free parameters against observational constraints.

**Entry point:** `main.py`

**Typical run (Mini-Millennium):**
```bash
cd autocalibration
mkdir -p minimillennium
python main.py \
  -c ../input/minimillennium.par \
  -b ../madsage \
  -S space_BHspin.txt \
  -o minimillennium \
  -m 20 \
  -s 12 \
  -x "SMF_z0,HIMF,CSFRDH" \
  -t chi2 \
  --seed 12345 > minimillennium/full_dump.txt
```

**Choosing the swarm size (`-s`):** set the number of particles equal to the number of tasks (CPU cores / SLURM `--ntasks`) available, so every particle in an iteration runs MAD-SAGE at the same time. More particles than tasks makes each iteration run in several sequential batches. On a laptop, use 2 particles. For a quick test, also use 1 iteration (`-m 1`), which only evaluates the random starting swarm (no PSO update) and takes about 4 minutes on Mini-Millennium:
```bash
cd autocalibration
mkdir -p minimillennium
python main.py -c ../input/minimillennium.par -b ../madsage -S space_BHspin.txt \
  -o minimillennium -s 2 -m 1 -x "SMF_z0,HIMF,CSFRDH" -t chi2 --seed 12345 \
  > minimillennium/full_dump.txt
```

**Key arguments:**
- `-c` — MAD-SAGE `.par` template (required); h, Omega, BoxLen, tree files, a_list and output prefix are read from it by `constraints.configure()`
- `-b` — MAD-SAGE binary (default: `../madsage`)
- `-S` — Search space file defining parameter bounds (default: `space_BHspin.txt` in `autocalibration/`)
- `-x` — Comma-separated constraints; supports domain `(low-high)` and weight `*N` suffixes
- `-t` — Goodness-of-fit test: `chi2` (reduced chi2; default and preferred) or `student-t` (negative log-likelihood, more tolerant of outlier bins)
- `-m` — Max PSO iterations (default: 10)
- `-s` — Swarm size (default: 10 + 2√D where D = number of free parameters); set equal to the number of tasks/cores, 2 on a laptop
- `--seed` — RNG seed; drawn from the OS and logged if omitted
- `-o` — Run directory: `madsage_pso.log`, `tracks/`, `temp/` (per-particle `.par` files `<par>_<pid>_<NNNN>_temp.par`, always kept), `MS_output_<pid>/` (MAD-SAGE output, deleted after each evaluation), `plots/` (from `post_analysis.py` and `best_fit_params.py --plot`)
- `-H`, `-v` — inherited from Shark; not functional for MAD-SAGE

**Objective:** ∏ stat_i^w_i over constraints (both statistics are positive, smaller is better).

**Search space file format** (`space_BHspin.txt`):
```
ParamName, PlotLabel, IsLog, LowerBound, UpperBound
GasSinkRate,fmove,1,0.01,0.999
EnergySN,eSN,1,5e49,3e51
```
Bounds are physical. `IsLog=1` means the PSO samples the parameter log-uniformly (search runs in log10); the temp `.par`, log and track files are always in physical units.

**Current free parameters:** `GasSinkRate` (gas inflow rate to BH), `EnergySN` (SN feedback energy).

**Available constraints** (defined in `constraints.py`, all z = 0):
- `SMF_z0` — Stellar mass function (Driver et al. 2022, GAMA); default domain log M* 8–13
- `HIMF` — HI mass function (Jones et al. 2018); default domain log M_HI 7–12
- `CSFRDH` — Cosmic star formation rate density history from stellar age bins (D'Silva et al. 2023)

**Script roles:**
- `pso.py` — PSO algorithm (adapted from PySwarm); saves per-iteration tracks to `output_dir/tracks/` (positions in physical units)
- `execution.py` — Writes a temp `.par` file per particle, runs MAD-SAGE, evaluates constraints, returns fitness
- `constraints.py` — Loads observation files and model snapshots, computes goodness-of-fit per constraint; `configure()` sets simulation info from the `.par` file
- `analysis.py` — Statistical tests (`chi2`, `studentT`), search space loader, log/physical transforms
- `common.py` — I/O utilities: reads MAD-SAGE `.par` files, counts tree files, loads observation data
- `best_fit_params.py` — Best-fit parameters and convergence plot from tracks
- `diagnostics.py` — Visualizes PSO convergence from tracks (`python diagnostics.py -S space_BHspin.txt minimillennium/tracks`)
- `post_analysis.py` — Reconstructs and plots model vs. observation fits from PSO log

**Python dependencies:** `numpy`, `scipy`, `h5py`, `matplotlib`; optional `pandas`/`seaborn` for pairplots.
