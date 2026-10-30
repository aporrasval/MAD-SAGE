# DECISIONS.md

Key design and implementation decisions with rationale. Append new entries chronologically; do not delete old ones.

---

### [2026-10-29] Autocalibration reads simulation info from the .par file

**Decision:** `constraints.configure()` takes h, Omega, BoxLen, tree-file range, a_list and output prefix from the DarkSAGE .par file given with `-c`, instead of a hard-coded `sim` switch.
**Why:** switching between Mini-Millennium (local) and TNG50-4 (cluster) required editing code, and the stale Mini-Millennium a_list path made every constraint crash. The .par file is already the source of truth for the run.
**Note:** BoxLen is in Mpc/h; volume = (BoxLen/h)^3 × files run / files total. Omega_Lambda is taken as 1 − Omega (flat), as before.

### [2026-10-29] PSO objective = ∏ stat_i^w_i

**Decision:** combine constraints as 10^(Σ w_i log10 stat_i).
**Why:** the previous 10^(Σ log10(stat_i w_i)) = ∏ stat_i × ∏ w_i, so weights never changed which particle was best. With all weights 1 the two forms agree, so unweighted chi2 runs are unchanged.

### [2026-10-29] Student-t statistic is a negative log-likelihood with a Gaussian fallback

**Decision:** `studentT` returns −Σ log t(σ_i; ν) with ν = 2 var/(var − 1) (variance matching, var = mean σ²); for var ≤ 1 it uses the Gaussian limit.
**Why:** the PSO minimises and the objective takes log10, so the statistic must be positive and smaller-is-better. A t-distribution cannot have variance ≤ 1, so the Gaussian (ν → ∞) is the consistent limit there.

### [2026-10-29] Log-uniform sampling for IsLog=1, physical units on disk

**Decision:** the swarm moves in log10 coordinates for `IsLog=1` parameters; the temporary .par file, track files and log are always in physical units.
**Why:** `IsLog` was ignored, so EnergySN (5e49–3e51) was sampled linearly and the low end was almost never explored. Keeping tracks physical leaves `best_fit_params.py`, `post_analysis.py`, `diagnostics.py` and older track files compatible.

### [2026-10-29] PSO runs are seeded

**Decision:** `--seed` sets `np.random.seed`; if omitted, a seed is drawn from the OS and logged.
**Why:** reproducibility requirement (CLAUDE.md): every calibration must be traceable to parameter set, code version and seed.

### [2026-10-29] chi2 is the default PSO statistic

**Decision:** `-t` defaults to `chi2` (reduced χ²); `student-t` remains available.
**Why:** user preference — χ² is simpler to interpret (≈1 is a good fit) and rewards matching every bin. Student-t down-weights outlier bins and its values have no absolute scale.
**Note:** `analysis.chi2` divides by (N_points − 2) per constraint. Because the objective is a product, this constant does not change which particle is best, only the reported value.
