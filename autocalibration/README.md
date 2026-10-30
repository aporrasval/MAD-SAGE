# MAD-SAGE autocalibration (PSO)

Particle Swarm Optimisation (PSO) of MAD-SAGE free parameters against z = 0
observational constraints. Adapted from the DarkSage 2024 and Shark PSO code (ICRAR/UWA) and PySwarm.

Each PSO particle is one set of free-parameter values. For every particle the
code writes a copy of your MAD-SAGE `.par` file with those values, runs
MAD-SAGE, compares the output with observations, and returns a single number
(smaller is better). The swarm then moves towards the best particles and the
process repeats for the requested number of iterations.

---

## 1. Before you run

1. **Build MAD-SAGE** in the repository root (`make`), so that `../madsage`
   exists. `main.py` uses this binary by default, so `-b` is optional.
2. **Python 3** with `numpy`, `scipy`, `h5py`, `matplotlib`
   (`pandas`/`seaborn` are optional, only for pair plots in `diagnostics.py`).
3. **Check that `output/routines.py` matches the C output struct** (it is used
   to read the galaxy files; see section 2).
4. **Edit the files listed in section 2** for your simulation and parameters.

---

## 2. Files you need to edit, and where

### Always check

| File | What to edit | Notes |
|---|---|---|
| `../input/<sim>.par` (e.g. `../input/minimillennium.par`) | `SimulationDir`, `TreeName`, `FileWithSnapList`, `FirstFile`/`LastFile`, `BoxLen`, `Hubble_h`, `Omega`, `FileNameGalaxies`, plus all physics settings | This is the **template** for every particle. The autocalibration reads h, Omega, BoxLen, the tree-file range, the a_list and the output prefix from it (`constraints.configure()`), so nothing about the simulation is hard-coded in Python. Relative paths are taken from the repository root, where `madsage` runs. |
| same `.par` file | An **uncommented** line for `OutputDir` and for **every** free parameter in the space file (e.g. `GasSinkRate`, `EnergySN`) | The parameter name must be the first word on the line. The values are overwritten for each particle; `OutputDir` is replaced by `<run dir>/MS_output_<pid>/`. If any of them is missing, the run stops with an error instead of silently using the template value. |
| `space_BHspin.txt` (or your own space file, passed with `-S`) | Which parameters are free and their bounds | Format below. `space_BHspin.txt` is the default for `main.py` and the post-run scripts, so `-S` is only needed for a different file. |

**Search space file** (one line per free parameter, no header line):
```
ParamName, PlotLabel, IsLog, LowerBound, UpperBound
GasSinkRate,fmove,1,0.01,0.999
EnergySN,eSN,1,5e49,3e51
```
- `ParamName` must match the name in the `.par` file exactly.
- Bounds are always in **physical** units.
- `IsLog=1`: the swarm searches in log10 (log-uniform sampling; use for
  parameters spanning more than ~1 dex). Bounds must then be > 0.
- The temp `.par` files, log and track files are always in physical units.

### Only if you change the simulation output or the C code

| File | Where | When |
|---|---|---|
| `../output/routines.py` | the numpy dtype in `galdtype_madsage()` | Whenever the `GALAXY_OUTPUT_LARGE` struct (written when `AgeStructOut 1`, as in `minimillennium.par`) in `code/core_allvars.h` / `code/core_save.c` changes. If it does not match byte-for-byte, every galaxy is read wrongly. |
| `constraints.py` | `_load_model_data()` (≈ line 158–162): `fields`, `Nage = 30`, `Nannuli=30`, and the snapshot name `'_z0.000'` | If the number of annuli or age bins in the `.par` file changes, or the z = 0 output file is named differently (e.g. another simulation's redshift list). |

### Only if you change or add a constraint

| File | Where | What |
|---|---|---|
| `constraints.py` | Binning constants at the top (≈ lines 36–66): `mbins`, `dm`, `Nmin` | Mass-function bins and the minimum number of galaxies per bin. |
| `constraints.py` | Default fit range `domain = (...)` in each class: `HIMF` (≈ line 278), `SMF_z0` (≈ line 308), `CSFRDH` (≈ line 342) | A narrower range can also be given on the command line, e.g. `SMF_z0(9-12)`, but it must lie inside the class default. |
| `constraints.py` | `get_obs_x_y_err()` of each class | Observation file name and columns. Observation files are read from **this directory** (`autocalibration/`). |
| `constraints.py` | New constraint: write a class with `domain`, `get_obs_x_y_err()` and `get_model_x_y()`, then register it in the `_constraints` dict inside `parse()` (≈ line 388) | The name in that dict is what you pass to `-x`. |
| `post_analysis.py` | `obs_files` dict (≈ line 85), `constraint_names` (≈ line 90), `plot_settings` (≈ line 169) | The plotting script identifies each constraint by its observation file name, so update these whenever an observation file is renamed or a constraint is added. |

### Usually not needed

| File | Where | When |
|---|---|---|
| `execution.py` | `run_madsage()`, command line (≈ lines 261–264) | Only to change how MAD-SAGE is launched (e.g. the commented-out `mpirun` line). Each particle currently runs one serial MAD-SAGE. |
| `analysis.py` | `chi2`, `studentT` | Only to change the goodness-of-fit statistic. |

Notes:
- `Driver_SFRD.dat` is an alternative CSFRD data set; it is only referenced in a commented-out line in `CSFRDH`.

---

## 3. Running

All commands are run from inside `autocalibration/`. The run directory (`-o`)
must exist before starting, because stdout is redirected into it.

### Recommended settings

- **Particles (`-s`) = number of tasks (CPU cores / SLURM `--ntasks`).**
  Every particle in an iteration then runs MAD-SAGE at the same time. With
  more particles than tasks, each iteration runs in several sequential batches
  and takes proportionally longer. With fewer, cores sit idle.
- **Iterations (`-m`) = 10** (the default). Wall time is roughly
  iterations × one MAD-SAGE run (≈ 4 min on Mini-Millennium), as long as
  particles ≤ tasks.
- **Statistic (`-t`) = `chi2`** (the default).
- **Always pass `--seed`** so the run can be reproduced exactly. (If omitted, a
  seed is drawn and written to the log.)

### Full run (10 iterations)

Replace `-s 16` with your number of tasks:
```bash
cd autocalibration
mkdir -p minimillennium
python main.py \
  -c ../input/minimillennium.par \
  -b ../madsage \
  -S space_BHspin.txt \
  -o minimillennium \
  -s 16 \
  -m 10 \
  -x "SMF_z0,HIMF,CSFRDH" \
  -t chi2 \
  --seed 12345 > minimillennium/full_dump.txt
```

### Quick test on a laptop (2 particles, 1 iteration)

One iteration only evaluates the random starting swarm (no PSO update). It is
for checking that everything runs and writes to the right place, and takes
about 4 minutes on Mini-Millennium (make sure you ran MAD SAGE on Mini-Millennium prior to this):
```bash
cd autocalibration
mkdir -p minimillennium
python main.py -c ../input/minimillennium.par -b ../madsage -S space_BHspin.txt \
  -o minimillennium -s 2 -m 1 -x "SMF_z0,HIMF,CSFRDH" -t chi2 --seed 12345 \
  > minimillennium/full_dump.txt
```
Clear the run directory (`/bin/rm -rf minimillennium`) before starting a new
run in it. Otherwise `madsage_pso.log` and `temp/` keep the old run's entries,
and old track files from later iterations are left next to the new ones.

### On a cluster (SLURM)

The particles are run with Python `multiprocessing`, so the job must be on a
**single node**, and `-s` should equal `--ntasks`:
```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00

cd /path/to/MAD_SAGE/autocalibration
mkdir -p minimillennium
python main.py -c ../input/minimillennium.par -b ../madsage -S space_BHspin.txt \
  -o minimillennium -s ${SLURM_NTASKS} -m 10 -x "SMF_z0,HIMF,CSFRDH" -t chi2 --seed 12345 \
  > minimillennium/full_dump.txt
```
The number of worker processes is `min(cores seen on the node, -s)`. Python
sees every core on the node, not only the ones SLURM allocated, so keep
`-s` equal to `--ntasks` to avoid oversubscribing.

### Command-line options

| Option | Meaning | Default |
|---|---|---|
| `-c` | MAD-SAGE `.par` template (required) | none |
| `-b`, `--madsage-binary` | MAD-SAGE executable. Optional: if omitted, `main.py` uses `madsage` in the repository root | `../madsage` (built by `make`) |
| `-S` | Search space file. Optional: if omitted, `main.py` uses `space_BHspin.txt` in `autocalibration/` | `space_BHspin.txt` |
| `-o` | Run directory | `.` |
| `-s` | Swarm size (number of particles) | 10 + 2√D (12 for 2 parameters) |
| `-m` | Number of PSO iterations | 10 |
| `-x` | Constraints, comma-separated. Optional fit range `(low-high)` and weight `*N`, e.g. `"SMF_z0(9-12)*2,HIMF,CSFRDH"` | `HIMF,SMF_z0,CSFRDH` |
| `-t` | Statistic: `chi2` (reduced χ²) or `student-t` (negative log-likelihood, more tolerant of outlier bins) | `chi2` |
| `--seed` | Random seed for the PSO | drawn from the OS and logged |
| `-k`, `-v`, `-H` and the HPC options | Inherited from Shark; not used by the MAD-SAGE runner | |

**Objective:** for each particle, ∏ stat_i^w_i over the constraints (the
product of the per-constraint statistics, each raised to its weight). Lower is
better.

**Available constraints** (all z = 0):

| Name | Observation | File | Default fit range |
|---|---|---|---|
| `SMF_z0` | Stellar mass function, Driver et al. (2022), GAMA | `GAMA_SMF_highres.csv` | log M* 8–13 |
| `HIMF` | HI mass function, Jones et al. (2018), ALFALFA | `Jones2018_HIMF.dat` | log M_HI 7–12 |
| `CSFRDH` | Cosmic SFR density history from stellar age bins, D'Silva et al. (2023) | `CSFH_DSILVA+23_Ver_Final.csv` | look-back time 0–14 Gyr |

---

## 4. Where the output goes

Everything from one run is inside the run directory (`-o`):
```
minimillennium/
├── full_dump.txt        stdout: model/observation arrays + "Particle ... evaluated to" lines
├── madsage_pso.log      run settings (seed, bounds, constraints) and results
├── tracks/              track_NNN_pos.npy (positions, physical units) and
│                        track_NNN_fx.npy (objective) for each iteration
├── temp/                <par>_<pid>_<NNNN>_temp.par: the exact .par file each
│                        particle ran with (always kept)
├── MS_output_<pid>/     MAD-SAGE output while a particle runs (deleted afterwards)
└── plots/               created by the post-run scripts below
```
Each `evaluated to` line in the log names its `.par` file in `temp/`, so any
result can be rerun by hand with `../madsage minimillennium/temp/<file>`.
(Its `OutputDir` points to the deleted `MS_output_<pid>/`; change it first.)

---

## 5. After the run

```bash
# Best-fit parameters and convergence plot  -> minimillennium/plots/pso_convergence.pdf
python best_fit_params.py -t minimillennium/tracks -S space_BHspin.txt --plot

# Model vs observations for every particle, best fit highlighted
#   -> minimillennium/plots/reconstructed_{SMF,HIMF,CSFRD}.pdf
python post_analysis.py -d minimillennium/full_dump.txt -t minimillennium/tracks -S space_BHspin.txt

# Swarm positions per iteration
python diagnostics.py -S space_BHspin.txt minimillennium/tracks
```
`-S` is optional in all three scripts (default `space_BHspin.txt`); pass it if the
run used a different space file.
To keep the best fit, copy its values into `../input/<sim>.par` and rerun
MAD-SAGE normally.

---

## 6. Scripts

| Script | Role |
|---|---|
| `main.py` | Entry point: reads options, sets up logging and the run directory, starts the PSO |
| `pso.py` | PSO algorithm (adapted from PySwarm); writes `tracks/` |
| `execution.py` | `run_madsage()`: writes the temp `.par`, runs MAD-SAGE, evaluates the constraints, returns the objective |
| `constraints.py` | Observations, model histograms and fit ranges for each constraint; `configure()` reads the simulation from the `.par` file |
| `analysis.py` | `chi2`, `studentT`, search-space loader, log/physical transforms |
| `common.py` | `.par` reader, tree-file counter, observation loader, helpers |
| `best_fit_params.py` | Best fit and convergence plot from `tracks/` |
| `post_analysis.py` | Model-vs-observation plots from `full_dump.txt` |
| `diagnostics.py` | Swarm diagnostics from `tracks/` |

## Known limitations

- Particles run in parallel write to the same stdout, so `full_dump.txt` blocks
  can interleave. `post_analysis.py` warns when this happens.
- Timestamps from worker processes are local time; those from `main.py` are UTC.
- The `evaluated to` line gives only the combined objective, not the
  statistic of each constraint.
