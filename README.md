# MAD-SAGE

**MAD-SAGE** is a fork of [DARK SAGE](https://github.com/arhstevens/DarkSage) — a semi-analytic model (SAM) of galaxy formation — extended with a Magnetically Arrested Disc (MAD) black hole accretion and spin evolution module.

## Overview

Semi-analytic models of galaxy formation evolve dark matter merger trees (from N-body simulations) by applying analytic prescriptions for baryonic physics: gas cooling, star formation, supernova feedback, AGN feedback, and galaxy mergers. DARK SAGE distinguishes itself by resolving disc structure across 30 radial annuli, tracking HI, H₂, and star formation radially rather than as single integrated quantities.

This fork adds a self-consistent model for supermassive black hole spin evolution under MAD accretion, enabling jet power to be computed from first principles and fed back into the host galaxy.

### What is a MAD disc?

In the Magnetically Arrested Disc regime, magnetic flux threading the black hole becomes dynamically important, halting infall and launching powerful relativistic jets. The MAD state produces the highest radiative and mechanical efficiencies achievable in accretion physics, making it a key ingredient for AGN feedback in massive galaxies.

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
# Build with make
make

# Or via Python wrapper (detects compiler, optionally patches parameter-dependent globals)
python compile.py

# Run in serial
./darksage input/minimillennium.par

# Run with MPI (experimental — each rank processes a separate merger tree file)
mpirun -np 8 ./darksage input/tng100.par
```

There is no automated test suite. Validation is done by running on the Mini-Millennium simulation and comparing galaxy property distributions (stellar mass function, HI mass function, cosmic star formation rate density) against observations.

## Parameter Files

Configuration files in `input/` control all model behavior. Key parameters:

| Parameter | Description |
|---|---|
| `SpinEvolutionOn` | Enable MAD spin evolution module (`1` = on) |
| `AccretionMode` | `0` = coherent, `1` = chaotic |
| `fEddThreshold` | Eddington ratio threshold separating thin-disc from thick-disc/ADAF regime |
| `BlackHoleSeedMass` | Initial BH seed mass in M☉/h |
| `GasSinkRate` | Fraction of disc gas driven inward per dynamical time |
| `EnergySN` | Supernova feedback energy per event (erg) |
| `RadiativeEfficiency` | Fixed η; only used when `SpinEvolutionOn=0` |

Parameter files are tracked in version control so all runs are reproducible.

## Output

Binary galaxy catalogs written per snapshot. Use `output/routines.py` to read output into numpy structured arrays — it defines the full dtype matching the `GALAXY` struct. Key black hole output fields:

- `BlackHoleMass` — BH mass (M☉/h)
- `BlackHoleSpin` — dimensionless spin parameter `a ∈ [−1, 1]`
- `BlackHoleMagnetization` — magnetic flux parameter
- `BlackHoleJetPower` — instantaneous jet mechanical power (erg/s)

## Autocalibration

`autocalibration/` contains a Particle Swarm Optimization (PSO) system to tune DarkSAGE free parameters against observational constraints.

```bash
cd autocalibration
python main.py \
  -c /path/to/parameter.par \
  -S space_GUREFT.txt \
  -o output_dir \
  -m 50 \
  -s 20 \
  -x "SMF_z6,QLF_z6,CSFRDH*2" \
  -t student-t
```

**Key arguments:**

| Flag | Description |
|---|---|
| `-c` | DarkSAGE `.par` config file (required) |
| `-S` | Search space file defining parameter bounds |
| `-x` | Constraints to fit; supports domain `(low-high)` and weight `*N` suffixes |
| `-t` | Goodness-of-fit statistic: `student-t` or `chi2` |
| `-m` | Maximum PSO iterations |
| `-s` | Swarm size (default: 10 + 2√D, D = free parameters) |
| `-H` | HPC/SLURM mode for cluster submission |

**Available constraints:**

| Constraint | Observation |
|---|---|
| `SMF_z6` | Stellar mass function at z=6 (Weibel et al. 2024) |
| `QLF_z6` | AGN bolometric luminosity function at z=6 (JWST + Shen et al. 2020) |
| `CSFRDH` | Cosmic star formation rate density history (D'Silva et al. 2023) |

**Search space file format** (`space_GUREFT.txt`):
```
ParamName, PlotLabel, IsLog, LowerBound, UpperBound
RadiativeEfficiency, eBH, 0, 0.0378, 0.4226
GasSinkRate, fmove, 1, 0.01, 0.999
EnergySN, eSN, 1, 5e49, 3e51
```
`IsLog=1` samples the parameter log-uniformly.

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
│   └── model_reincorporation.c         # Hot gas reincorporation
├── input/              # Parameter files (.par) for supported simulations
├── autocalibration/    # PSO-based parameter calibration system
├── output/
│   └── routines.py    # Python reader for binary galaxy catalogs
└── Makefile
```

The model reads N-body merger trees, then evolves galaxies forward through simulation snapshots. Each snapshot applies the physics modules in sequence: cooling → infall → star formation → disc instability → mergers → BH accretion (MAD) → feedback → reincorporation.

## References

- **DARK SAGE (current version):** Stevens et al. (2024), PASA 41:53 — [ADS](https://ui.adsabs.harvard.edu/abs/2024PASA...41...53S/abstract)
- **DARK SAGE (original):** Stevens, Croton & Mutch (2016), MNRAS 461:859 — [ADS](http://adsabs.harvard.edu/abs/2016MNRAS.461..859S)
- **SAGE:** Croton et al. (2016) — [GitHub](https://github.com/darrencroton/sage/)
- **Post-merger spin:** Rezzolla et al. (2008), ApJL 679:L37
- **SMF at z=6:** Weibel et al. (2024)
- **QLF at z=6:** JWST observations + Shen et al. (2020)
- **CSFRD history:** D'Silva et al. (2023)

## License

See `LICENCE` file. This fork inherits the DARK SAGE license. The MAD physics additions (`model_mad.c`) are released under the same terms.
