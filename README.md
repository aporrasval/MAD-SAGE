# MAD-SAGE

**MAD-SAGE** is a fork of [DARK SAGE](https://github.com/arhstevens/DarkSage) — a semi-analytic model (SAM) of galaxy formation — extended with a Magnetically Arrested Disc (MAD) black hole accretion and spin evolution module.

## Overview

Semi-analytic models of galaxy formation evolve dark matter merger trees (from N-body simulations) by applying analytic prescriptions for baryonic physics: gas cooling, star formation, supernova feedback, AGN feedback, and galaxy mergers. DARK SAGE distinguishes itself by resolving disk structure across 30 radial annuli, tracking HI, H₂, and star formation radially rather than as single integrated quantities.

This fork adds a self-consistent model for supermassive black hole spin evolution under MAD accretion, enabling jet power to be computed from first principles and fed back into the host galaxy.

### What is a MAD disk?

In the Magnetically Arrested Disk regime, magnetic flux threading the black hole becomes dynamically important, halting infall and launching powerful relativistic jets. The MAD state produces the highest radiative and mechanical efficiencies achievable in accretion physics, making it a key ingredient for AGN feedback in massive galaxies.

## New Physics (`model_mad.c`)

The central contribution of this fork. Implements:

- **Spin evolution** — RK45 adaptive numerical integration of the dimensionless BH spin parameter `a` during accretion, accounting for torques from accreted angular momentum and the Bardeen-Petterson effect
- **Accretion modes** — coherent (prograde/retrograde) vs. chaotic accretion; sub- and super-Eddington regimes handled separately
- **ISCO radius** — spin-dependent innermost stable circular orbit `r_ISCO(a)`, used to compute radiative efficiency
- **Radiative efficiency** — `η(a)` derived from the ISCO, replacing a fixed free parameter during spin-active runs
- **Jet power** — magnetically driven jets computed from the MAD magnetization state and spin
- **Magnetization dynamics** — tracks the magnetic flux threading the black hole and its evolution through accretion and spin-down
- **Post-merger spin** — spin of BH remnants after binary mergers computed using the Rezzolla et al. (2008) fitting formulae

These fields are stored in the `GALAXY` struct (`core_allvars.h`) and written to the binary output.

## Build and Run

**Requirements:** GCC (Linux) or Clang (macOS, auto-detected), [GSL](https://www.gnu.org/software/gsl/) ≥ 2.0

```bash
# Build: produces the executable ./madsage in the repository root
make

# Run in serial
./madsage input/minimillennium.par

# Run with MPI (experimental — each rank processes a separate merger tree file)
mpirun -np 8 ./madsage input/minimillennium.par
```

Before rerunning into the same output directory, delete the old galaxy files (e.g. `/bin/rm -f output/results/model_z*`). Stale files from a build with a different output struct cannot be read by `output/routines.py`.

There is no automated test suite. Validation is done by running on the Mini-Millennium simulation and comparing galaxy property distributions (stellar mass function, HI mass function, cosmic star formation rate density) against observations. A Mini-Millennium run takes about 4 minutes on a laptop.

## Parameter Files

Configuration files in `input/` control all model behavior. `input/minimillennium.par` is the reference file. Its paths are relative to the repository root (where `madsage` is run): it expects the Mini-Millennium trees in `mini_millennium/treedata/` and the snapshot list in `mini_millennium/millennium.a_list`, and writes to `output/results/`. The trees are not included in the repository (152 MB); download them separately, or change `SimulationDir` and `FileWithSnapList` to where yours are. Key MAD and calibration parameters:

| Parameter | Description |
|---|---|
| `SpinEvolutionOn` | `1` = BH spin evolves (MAD module on); `0` = BHs have no spin |
| `SpinSeedValue` | Initial BH spin |
| `fEdd_ADAF_threshold` | Eddington ratio of the transition between the thin-disk (MAD) and ADAF regimes |
| `BH_SeedMass` | Initial BH mass, in code units of 10¹⁰ M☉/h (1e-8 = 100 M☉/h) |
| `BH_SeedMinHaloMass` | Minimum halo mass for BH seeding (M☉) |
| `GasSinkRate` | Fraction of disc gas driven inward per dynamical time (calibrated) |
| `EnergySN` | Supernova feedback energy per event in erg (calibrated) |
| `RadiativeEfficiency` | Fixed η; only used when `SpinEvolutionOn=0` |

## Output

Binary galaxy catalogs written per snapshot, named `<FileNameGalaxies>_z<redshift>_<file number>`. Use `output/routines.py` to read them into numpy structured arrays (`routines.madsage_snap(...)`); its dtype must match the `GALAXY_OUTPUT_LARGE` struct in `code/core_allvars.h` byte for byte (the struct written when `AgeStructOut 1`). Masses are in 10¹⁰ M☉/h. Key black hole output fields:

- `BlackHoleMass` — BH mass
- `BlackHoleSpin` — dimensionless spin parameter `a ∈ [−0.998, 0.998]`
- `BH_Magnetization` — magnetic flux parameter
- `BH_JetPower`, `BH_JetEfficiency` — jet power and jet efficiency
- `BH_RadiativeEfficiency` — spin-dependent radiative efficiency η(a)

## Autocalibration

`autocalibration/` contains a Particle Swarm Optimization (PSO) system that tunes MAD-SAGE free parameters (currently `GasSinkRate` and `EnergySN`) against the z = 0 stellar mass function, HI mass function and cosmic star formation rate density history. **See [`autocalibration/README.md`](autocalibration/README.md)** for the full guide: which files to edit and where, all options, output layout and post-run analysis.

Recommended settings:
- **Particles (`-s`) = number of tasks** (CPU cores / SLURM `--ntasks`, single node), so every particle in an iteration runs at the same time. On a laptop, use 2.
- **Iterations (`-m`) = 10.** Wall time ≈ iterations × one MAD-SAGE run (≈ 4 min on Mini-Millennium).
- **Statistic `-t chi2`** (default), and always pass **`--seed`** for reproducibility.

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

Quick test on a laptop (2 particles, 1 iteration, ≈ 4 minutes): use `-s 2 -m 1`.

Everything for a run is written inside the `-o` directory: `madsage_pso.log`, `full_dump.txt`, `tracks/`, `temp/` (the `.par` file of every particle), `plots/`, and `MS_output_<pid>/` while a particle runs.

## Code Architecture

```
MAD-SAGE/
├── code/               # C source: physics modules and core infrastructure
│   ├── main.c                          # Entry point
│   ├── core_build_model.c              # Main evolution loop per merger tree
│   ├── core_allvars.c/.h               # Global variables and GALAXY struct
│   ├── core_proto.h                    # All function prototypes
│   ├── core_read_parameter_file.c      # .par file parser
│   ├── core_io_tree.c                  # Merger tree reader
│   ├── core_save.c                     # Binary output writer
│   ├── model_mad.c                     # MAD BH spin/jet physics (this fork)
│   ├── model_cooling.c                 # Gas cooling
│   ├── model_infall.c                  # Gas infall
│   ├── model_starformation_and_feedback.c  # Star formation and SN feedback
│   ├── model_disk_instability.c        # Bar and disc instability
│   ├── model_mergers.c                 # Galaxy mergers
│   ├── model_reincorporation.c         # Hot gas reincorporation
│   └── model_misc.c                    # Helpers (incl. BH seeding)
├── input/              # Parameter files (.par); minimillennium.par is the reference
├── autocalibration/    # PSO parameter calibration (see autocalibration/README.md)
├── output/
│   └── routines.py    # Python reader for binary galaxy catalogs
├── extra/CoolFunctions/  # Cooling function tables
└── Makefile            # Builds ./madsage
```

The model reads N-body merger trees, then evolves galaxies forward through simulation snapshots. Each snapshot applies the physics modules in sequence: cooling → infall → star formation → disc instability → mergers → BH accretion (MAD) → feedback → reincorporation.

## References

- **DARK SAGE (current version):** Stevens et al. (2024), PASA 41:53 — [ADS](https://ui.adsabs.harvard.edu/abs/2024PASA...41...53S/abstract)
- **DARK SAGE (original):** Stevens, Croton & Mutch (2016), MNRAS 461:859 — [ADS](http://adsabs.harvard.edu/abs/2016MNRAS.461..859S)
- **SAGE:** Croton et al. (2016) — [GitHub](https://github.com/darrencroton/sage/)
- **Post-merger spin:** Rezzolla et al. (2008), ApJL 679:L37
- **SMF at z=0:** Driver et al. (2022), MNRAS 513:439
- **HIMF at z=0:** Jones et al. (2018)
- **CSFRD history:** D'Silva et al. (2023)

## License

See `LICENCE` file. This fork inherits the DARK SAGE MIT license. The MAD physics additions (`model_mad.c`) are released under the same terms. Scripts in `autocalibration/` adapted from the Shark PSO code and pyswarm keep their original licences (GPL-3.0 / LGPL-2.1 for the Shark files), as listed in `LICENCE` and in each file's header.
