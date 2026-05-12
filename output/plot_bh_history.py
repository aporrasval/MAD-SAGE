"""
plot_bh_history.py

For each target log10(M_BH/Msun) in [10, 9, 8, 7, 6], finds the z=0 black hole
whose mass is closest to that target, traces its main progenitor across all
snapshots via RootGalaxyIndex, and plots BH mass, spin, and Eddington ratio
vs cosmic time on a 3-panel figure.

Progenitor linking: all progenitors of the same z=0 galaxy share the same
RootGalaxyIndex. At each snapshot we pick the Type=0 (central) galaxy with
that RootGalaxyIndex and the largest BlackHoleMass.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.integrate import quad
import os, sys

sys.path.insert(0, os.path.dirname(__file__))
from routines import galdtype_darksage

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
OUTPUT_DIR  = os.path.join(os.path.dirname(__file__), 'history_results')
SNAPLIST    = '/Users/porrasaj/Documents/Astronomy_research/DarkSAGE/running_test_minimillennium/mini_millennium/millennium.a_list'
N_FILES     = 8
H_h         = 0.73          # Hubble_h
OMEGA_M     = 0.25
OMEGA_L     = 0.75
H0_KMS_MPC  = 73.0          # km/s/Mpc
fEdd_ADAF   = 0.01          # ADAF threshold in our model

# Target log10(M_BH/Msun) values to trace
LOG_TARGETS = [10, 9, 8, 7, 6]

# ---------------------------------------------------------------------------
# Cosmology: age of universe at scale factor a  [Gyr]
# ---------------------------------------------------------------------------
_H0_inv_gyr = 1.0 / (H0_KMS_MPC * 1e3 / 3.0857e22) / 3.1536e16

def _integrand(ap):
    return 1.0 / (ap * np.sqrt(OMEGA_M / ap**3 + OMEGA_L))

def age_gyr(a):
    val, _ = quad(_integrand, 1e-10, a)
    return _H0_inv_gyr * val

# ---------------------------------------------------------------------------
# Read a single output file (AgeStructOut=1 → GALAXY_OUTPUT_LARGE dtype)
# ---------------------------------------------------------------------------
dtype = galdtype_darksage(Nannuli=30, Nage=30)

def read_file(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    with open(path, 'rb') as f:
        tmp = np.fromfile(f, dtype=np.int32, count=1)
        if len(tmp) == 0:
            return None
        n_trees = tmp[0]
        n_gals  = np.fromfile(f, dtype=np.int32, count=1)[0]
        _       = np.fromfile(f, dtype=np.int32, count=n_trees)
        if n_gals > 0:
            return np.fromfile(f, dtype=dtype, count=n_gals)
    return None

# ---------------------------------------------------------------------------
# Load scale factors and compute ages
# ---------------------------------------------------------------------------
scale_factors = np.loadtxt(SNAPLIST)
n_snaps       = len(scale_factors)
redshifts     = 1.0 / scale_factors - 1.0
ages_gyr      = np.array([age_gyr(a) for a in scale_factors])

print(f"Loaded {n_snaps} snapshots  (z=0 at {ages_gyr[-1]:.2f} Gyr)")

# ---------------------------------------------------------------------------
# Step 1: load all z=0 galaxies and find one BH per target mass
# ---------------------------------------------------------------------------
z0_z   = redshifts[-1]
z0_str = f'{z0_z:.3f}'
z0_gals = []

for fi in range(N_FILES):
    g = read_file(os.path.join(OUTPUT_DIR, f'model_z{z0_str}_{fi}'))
    if g is not None:
        z0_gals.append(g)

z0_all = np.concatenate(z0_gals)
bh_mass_msun = z0_all['BlackHoleMass'] * 1e10 / H_h

# For each target, pick the z=0 BH whose log mass is closest to the target
targets = []
for log_target in LOG_TARGETS:
    target_msun = 10.0**log_target
    valid = bh_mass_msun > 0
    if not np.any(valid):
        print(f"  WARNING: no BHs with positive mass for log target {log_target}")
        continue
    idx = np.argmin(np.abs(np.log10(bh_mass_msun[valid]) - log_target))
    full_idx = np.where(valid)[0][idx]
    gal = z0_all[full_idx]
    actual_log = np.log10(bh_mass_msun[full_idx])
    print(f"  log10 target {log_target}: found M_BH={bh_mass_msun[full_idx]:.3e} Msun "
          f"(log={actual_log:.2f})  spin={gal['BlackHoleSpin']:.3f}  "
          f"RootGalaxyIndex={gal['RootGalaxyIndex']}")
    targets.append({
        'log_target': log_target,
        'root_id':    gal['RootGalaxyIndex'],
        'label':      rf'$\log M_{{\bullet}}={log_target}$',
    })

# ---------------------------------------------------------------------------
# Step 2: load all snapshots once, then trace each target
# ---------------------------------------------------------------------------
print("\nLoading all snapshots...")
snap_data = []  # list of (t, z, all_galaxies_array)

for snap_i, (a, z, t) in enumerate(zip(scale_factors, redshifts, ages_gyr)):
    z_str     = f'{z:.3f}'
    snap_gals = []
    for fi in range(N_FILES):
        g = read_file(os.path.join(OUTPUT_DIR, f'model_z{z_str}_{fi}'))
        if g is not None:
            snap_gals.append(g)
    if snap_gals:
        snap_data.append((t, z, np.concatenate(snap_gals)))

print(f"Loaded data for {len(snap_data)} snapshots")

# ---------------------------------------------------------------------------
# Step 3: trace progenitors for each target
# ---------------------------------------------------------------------------
colors  = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
histories = []

for tgt in targets:
    root_id = tgt['root_id']
    hist    = []

    for (t, z, all_g) in snap_data:
        # Prefer central (Type=0) with matching root
        mask = (all_g['RootGalaxyIndex'] == root_id) & (all_g['Type'] == 0)
        prog = all_g[mask]
        if len(prog) == 0:
            mask = all_g['RootGalaxyIndex'] == root_id
            prog = all_g[mask]
        if len(prog) == 0:
            continue

        best = prog[np.argmax(prog['BlackHoleMass'])]
        hist.append({
            't':    t,
            'M_BH': best['BlackHoleMass'] * 1e10 / H_h,
            'spin': best['BlackHoleSpin'],
            'Edd':  best['Eddrate'],
        })

    hist.sort(key=lambda x: x['t'])
    print(f"  log target {tgt['log_target']}: traced {len(hist)} snapshots")
    histories.append(hist)

# ---------------------------------------------------------------------------
# Step 4: plot
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
fig.subplots_adjust(hspace=0.05)

for hist, tgt, color in zip(histories, targets, colors):
    if not hist:
        continue
    t_arr    = np.array([h['t']    for h in hist])
    M_arr    = np.array([h['M_BH'] for h in hist])
    spin_arr = np.array([h['spin'] for h in hist])
    Edd_arr  = np.array([h['Edd']  for h in hist])
    Edd_arr  = np.where(Edd_arr > 0, Edd_arr, np.nan)

    axes[0].semilogy(t_arr, M_arr,    color=color, lw=1.5, label=tgt['label'])
    axes[1].plot    (t_arr, spin_arr, color=color, lw=1.5)
    axes[2].semilogy(t_arr, Edd_arr,  color=color, lw=1.0)

# --- Panel 1: BH mass ---
axes[0].set_ylabel(r'$M_\bullet\ [M_\odot]$', fontsize=12)
axes[0].set_ylim(1e4, 1e11)
axes[0].legend(fontsize=9, loc='upper left', ncol=2)
axes[0].yaxis.set_major_formatter(ticker.LogFormatterMathtext())

# --- Panel 2: spin ---
axes[1].axhline(0, color='k', lw=0.6, ls=':')
axes[1].set_ylabel(r'$a_\bullet$', fontsize=12)
axes[1].set_ylim(-1.05, 1.05)

# --- Panel 3: Eddington ratio ---
axes[2].axhspan(0.3 * fEdd_ADAF, fEdd_ADAF, color='black', alpha=0.4, zorder=0,
                label='ADAF transition')
axes[2].set_ylabel(r'$\dot{m}$', fontsize=12)
axes[2].set_xlabel(r'$t_{\rm ABB}\ [\rm Gyr]$', fontsize=12)
axes[2].set_ylim(1e-6, 1e2)

# --- Shared redshift axis on top ---
ax_top  = axes[0].twiny()
z_ticks = [10, 6, 4, 2, 1, 0.5, 0]
t_ticks = [age_gyr(1.0/(1+z)) if z < 100 else ages_gyr[0] for z in z_ticks]
ax_top.set_xlim(axes[0].get_xlim())
ax_top.set_xticks(t_ticks)
ax_top.set_xticklabels([str(z) for z in z_ticks])
ax_top.set_xlabel('z', fontsize=12)

# --- x-axis limits ---
for ax in axes:
    ax.set_xlim(ages_gyr[0], ages_gyr[-1])

out_base = os.path.join(os.path.dirname(__file__), 'bh_history')
plt.savefig(out_base + '.pdf', bbox_inches='tight', dpi=150)
plt.savefig(out_base + '.png', bbox_inches='tight', dpi=150)
print(f"\nSaved bh_history.pdf and bh_history.png")
plt.show()
