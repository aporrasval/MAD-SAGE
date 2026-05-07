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
./darksage input/minimillennium.par

# Run with MPI (experimental)
mpirun -np 8 ./darksage input/tng100.par
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

`autocalibration/` contains a Particle Swarm Optimization (PSO) system to tune DarkSAGE free parameters against observational constraints.

**Entry point:** `main.py`

**Typical run:**
```bash
cd autocalibration
python main.py \
  -c /path/to/parameter.par \
  -S space_GUREFT.txt \
  -o output_dir \
  -m 50 \
  -s 20 \
  -v 0-31 \
  -x "SMF_z6,QLF_z6,CSFRDH*2" \
  -t student-t
```

**Key arguments:**
- `-c` — DarkSAGE `.par` config file (required)
- `-S` — Search space file defining parameter bounds (default: `space.txt`)
- `-x` — Comma-separated constraints; supports domain `(low-high)` and weight `*N` suffixes
- `-t` — Goodness-of-fit test: `student-t` or `chi2`
- `-m` — Max PSO iterations (default: 10)
- `-s` — Swarm size (default: 10 + 2√D where D = number of free parameters)
- `-H` — HPC/SLURM mode for cluster submission

**Search space file format** (`space_GUREFT.txt`):
```
ParamName, PlotLabel, IsLog, LowerBound, UpperBound
RadiativeEfficiency, eBH, 0, 0.0378, 0.4226
GasSinkRate, fmove, 1, 0.01, 0.999
EnergySN, eSN, 1, 5e49, 3e51
```
`IsLog=1` means the parameter is sampled log-uniformly.

**Current free parameters:** `RadiativeEfficiency` (BH accretion efficiency), `GasSinkRate` (gas inflow rate to BH), `EnergySN` (SN feedback energy).

**Available constraints** (defined in `constraints.py`):
- `SMF_z6` — Stellar mass function at z=6 (Weibel et al. 2024)
- `QLF_z6` — AGN bolometric luminosity function at z=6 (JWST + Shen+20)
- `CSFRDH` — Cosmic star formation rate density history (D'Silva et al. 2023)

**Script roles:**
- `pso.py` — PSO algorithm (adapted from PySwarm); saves per-iteration tracks to `output_dir/tracks/`
- `execution.py` — Writes a temp `.par` file per particle, runs DarkSAGE, evaluates constraints, returns fitness
- `constraints.py` — Loads observation CSVs and model snapshots, computes goodness-of-fit per constraint
- `analysis.py` — Statistical tests (`chi2`, `studentT`) and search space loader
- `common.py` — I/O utilities: reads `.par` files, observation data, HDF5 galaxy catalogs
- `diagnostics.py` — Visualizes PSO convergence from tracks (`python diagnostics.py -S space.txt tracks/`)
- `post_analysis.py` — Reconstructs and plots model vs. observation fits from PSO log

**Python dependencies:** `numpy`, `scipy`, `h5py`, `matplotlib`; optional `pandas`/`seaborn` for pairplots.
