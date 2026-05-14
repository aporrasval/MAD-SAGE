"""
plot_bh_history.py

Traces the main progenitor of each target BH using the halo merger tree
FirstProgenitor chain — avoiding the RootGalaxyIndex fallback that was
jumping to unrelated galaxies sharing the same root.

Algorithm:
  1. At z=0 find the BH closest to each target log mass.
  2. Decode its GalaxyIndex → filenr, tree, galaxy_nr.
  3. Read the tree file; walk FirstProgenitor backward to get
     (SnapNum, local_halo_idx) for every ancestor.
  4. At each snapshot load model_z{z}_{filenr}, match by
     TreeIndex == tree AND HaloIndex == local_halo_idx.

Plots BH mass, spin, and Eddington ratio vs cosmic time for
log10(M_BH/Msun) targets = 10, 9, 8, 7, 6.
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
TREE_DIR    = '/Users/porrasaj/Documents/Astronomy_research/DarkSAGE/running_test_minimillennium/mini_millennium/treedata'
TREE_NAME   = 'trees_063'
N_FILES     = 8
H_h         = 0.73
OMEGA_M     = 0.25
OMEGA_L     = 0.75
H0_KMS_MPC  = 73.0
fEdd_ADAF   = 0.01

LOG_TARGETS = [10, 9, 8, 7, 6]

# GalaxyIndex encoding constants (from core_save.c)
TREE_MUL_FAC   = 1000000000
FILENR_MUL_FAC = 10000000000000000

# ---------------------------------------------------------------------------
# Cosmology
# ---------------------------------------------------------------------------
_H0_inv_gyr = 1.0 / (H0_KMS_MPC * 1e3 / 3.0857e22) / 3.1536e16

def _integrand(ap):
    return 1.0 / (ap * np.sqrt(OMEGA_M / ap**3 + OMEGA_L))

def age_gyr(a):
    val, _ = quad(_integrand, 1e-10, a)
    return _H0_inv_gyr * val

# ---------------------------------------------------------------------------
# Galaxy dtype (AgeStructOut=1)
# ---------------------------------------------------------------------------
gal_dtype = galdtype_darksage(Nannuli=30, Nage=30)

def read_gal_file(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    with open(path, 'rb') as f:
        tmp = np.fromfile(f, dtype=np.int32, count=1)
        if len(tmp) == 0:
            return None
        n_trees = tmp[0]
        n_gals  = np.fromfile(f, dtype=np.int32, count=1)[0]
        np.fromfile(f, dtype=np.int32, count=n_trees)
        if n_gals > 0:
            return np.fromfile(f, dtype=gal_dtype, count=n_gals)
    return None

# ---------------------------------------------------------------------------
# Halo tree dtype  (104 bytes, matches core_simulation.h)
# ---------------------------------------------------------------------------
halo_dtype = np.dtype([
    ('Descendant',          np.int32),
    ('FirstProgenitor',     np.int32),
    ('NextProgenitor',      np.int32),
    ('FirstHaloInFOFgroup', np.int32),
    ('NextHaloInFOFgroup',  np.int32),
    ('Len',                 np.int32),
    ('Mvir1',  np.float32), ('Mvir2', np.float32), ('Mvir3', np.float32),
    ('Pos',    np.float32, 3),
    ('Vel',    np.float32, 3),
    ('VelDisp', np.float32),
    ('Vmax',    np.float32),
    ('Spin',    np.float32, 3),
    ('MostBoundID', np.int64),
    ('SnapNum',     np.int32),
    ('FileNr',      np.int32),
    ('SubhaloIndex', np.int32),
    ('SubHalfMass',  np.float32),
])
assert halo_dtype.itemsize == 104, f"halo_dtype itemsize={halo_dtype.itemsize}, expected 104"

def read_tree_file(filenr):
    path = os.path.join(TREE_DIR, f'{TREE_NAME}.{filenr}')
    with open(path, 'rb') as f:
        n_trees   = np.fromfile(f, dtype=np.int32, count=1)[0]
        tot_halos = np.fromfile(f, dtype=np.int32, count=1)[0]
        tree_nhalo = np.fromfile(f, dtype=np.int32, count=n_trees)
        halos = np.fromfile(f, dtype=halo_dtype, count=tot_halos)
    # Cumulative offset of each tree's first halo in the flat halos array
    tree_first = np.concatenate([[0], np.cumsum(tree_nhalo[:-1])])
    return halos, tree_first

def decode_galaxy_index(gal_idx):
    """Decode GalaxyIndex → (filenr, tree_idx, galaxy_nr)."""
    filenr_enc = int(gal_idx // FILENR_MUL_FAC)
    tree_enc   = int((gal_idx % FILENR_MUL_FAC) // TREE_MUL_FAC)
    galaxy_nr  = int(gal_idx % TREE_MUL_FAC)
    if   filenr_enc < 100: filenr, tree = filenr_enc,       tree_enc
    elif filenr_enc < 200: filenr, tree = filenr_enc - 100, tree_enc + 10000000
    elif filenr_enc < 300: filenr, tree = filenr_enc - 200, tree_enc + 20000000
    elif filenr_enc < 400: filenr, tree = filenr_enc - 300, tree_enc + 30000000
    else:                  filenr, tree = filenr_enc - 400, tree_enc + 40000000
    return filenr, tree, galaxy_nr

def get_progenitor_chain(halos, tree_first, tree_idx, halo_nr):
    """Walk FirstProgenitor backward from (tree_idx, halo_nr).
    Returns list of (snap_num, local_halo_idx) ordered z=0 → high-z."""
    chain = []
    local_idx = int(halo_nr)
    offset    = int(tree_first[tree_idx])
    while local_idx >= 0:
        h = halos[offset + local_idx]
        chain.append((int(h['SnapNum']), local_idx))
        nxt = int(h['FirstProgenitor'])
        if nxt < 0:
            break
        local_idx = nxt
    return chain

# ---------------------------------------------------------------------------
# Load snaplist
# ---------------------------------------------------------------------------
scale_factors = np.loadtxt(SNAPLIST)
n_snaps       = len(scale_factors)
redshifts     = 1.0 / scale_factors - 1.0
ages_gyr      = np.array([age_gyr(a) for a in scale_factors])
print(f"Loaded {n_snaps} snapshots  (z=0 at {ages_gyr[-1]:.2f} Gyr)")

# snap_num → redshift string for output filename lookup
snap_to_zstr = {i: f'{redshifts[i]:.3f}' for i in range(n_snaps)}

# ---------------------------------------------------------------------------
# Step 1: find z=0 targets
# ---------------------------------------------------------------------------
z0_str  = f'{redshifts[-1]:.3f}'
z0_gals = []
z0_fi   = []  # which file each galaxy came from

for fi in range(N_FILES):
    g = read_gal_file(os.path.join(OUTPUT_DIR, f'model_z{z0_str}_{fi}'))
    if g is not None:
        z0_gals.append(g)
        z0_fi.extend([fi] * len(g))

z0_all  = np.concatenate(z0_gals)
z0_fi   = np.array(z0_fi, dtype=np.int32)
bh_msun = z0_all['BlackHoleMass'] * 1e10 / H_h

targets = []
for log_target in LOG_TARGETS:
    valid = bh_msun > 0
    idx   = np.argmin(np.abs(np.log10(bh_msun[valid]) - log_target))
    full_idx = np.where(valid)[0][idx]
    gal  = z0_all[full_idx]
    filenr, tree_idx, galaxy_nr = decode_galaxy_index(int(gal['GalaxyIndex']))
    halo_nr = int(gal['HaloIndex'])
    print(f"  log target {log_target}: M_BH={bh_msun[full_idx]:.3e} Msun "
          f"(log={np.log10(bh_msun[full_idx]):.2f})  spin={gal['BlackHoleSpin']:.3f}  "
          f"filenr={filenr}  tree={tree_idx}  HaloIndex={halo_nr}")
    targets.append({
        'log_target': log_target,
        'label':      rf'$\log M_{{\bullet}}={log_target}$',
        'filenr':     filenr,
        'tree_idx':   tree_idx,
        'halo_nr':    halo_nr,
    })

# ---------------------------------------------------------------------------
# Step 2: for each target, walk the halo tree and build the history
# ---------------------------------------------------------------------------
# Group targets by filenr to load each tree file only once
from collections import defaultdict
by_filenr = defaultdict(list)
for i, tgt in enumerate(targets):
    by_filenr[tgt['filenr']].append(i)

# Pre-load per-filenr snapshot caches: snap_num → galaxy array
# (only the files we actually need)
needed_filenrs = set(tgt['filenr'] for tgt in targets)

# For each filenr, pre-load all snapshot galaxy files
snap_cache = {}  # (filenr, snap_num) → galaxy array or None
print("\nPre-loading snapshot files for needed tree files...")
for filenr in needed_filenrs:
    for snap_num in range(n_snaps):
        z_str = snap_to_zstr[snap_num]
        path  = os.path.join(OUTPUT_DIR, f'model_z{z_str}_{filenr}')
        snap_cache[(filenr, snap_num)] = read_gal_file(path)
print(f"  Loaded {len(snap_cache)} snapshot-file slots "
      f"({sum(v is not None for v in snap_cache.values())} non-empty)")

histories = []

for i, tgt in enumerate(targets):
    filenr  = tgt['filenr']
    tree_idx = tgt['tree_idx']
    halo_nr  = tgt['halo_nr']

    halos, tree_first = read_tree_file(filenr)
    chain = get_progenitor_chain(halos, tree_first, tree_idx, halo_nr)
    print(f"  log target {tgt['log_target']}: progenitor chain has {len(chain)} halos")

    hist = []
    for snap_num, local_halo_idx in chain:
        gals = snap_cache.get((filenr, snap_num))
        if gals is None:
            continue
        # Match by TreeIndex AND HaloIndex — uniquely identifies the galaxy
        mask = (gals['TreeIndex'] == tree_idx) & (gals['HaloIndex'] == local_halo_idx)
        matches = gals[mask]
        if len(matches) == 0:
            continue
        if len(matches) > 1:
            # Multiple galaxies in same halo (central + orphan satellites); prefer Type=0
            central = matches[matches['Type'] == 0]
            best = central[0] if len(central) > 0 else matches[0]
        else:
            best = matches[0]

        t = ages_gyr[snap_num]
        hist.append({
            't':    t,
            'M_BH': best['BlackHoleMass'] * 1e10 / H_h,
            'spin': best['BlackHoleSpin'],
            'Edd':  best['Eddrate'],
        })

    hist.sort(key=lambda x: x['t'])
    print(f"            → matched galaxies in {len(hist)} snapshots")
    histories.append(hist)

# ---------------------------------------------------------------------------
# Step 3: plot
# ---------------------------------------------------------------------------
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

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

# Panel 1: BH mass
axes[0].set_ylabel(r'$M_\bullet\ [M_\odot]$', fontsize=12)
axes[0].set_ylim(1e4, 1e11)
axes[0].legend(fontsize=9, loc='upper left', ncol=2)
axes[0].yaxis.set_major_formatter(ticker.LogFormatterMathtext())

# Panel 2: spin
axes[1].axhline(0, color='k', lw=0.6, ls=':')
axes[1].set_ylabel(r'$a_\bullet$', fontsize=12)
axes[1].set_ylim(-1.05, 1.05)

# Panel 3: Eddington ratio
axes[2].axhspan(0.3 * fEdd_ADAF, fEdd_ADAF, color='black', alpha=0.4, zorder=0,
                label='ADAF transition')
axes[2].set_ylabel(r'$\dot{m}$', fontsize=12)
axes[2].set_xlabel(r'$t_{\rm ABB}\ [\rm Gyr]$', fontsize=12)
axes[2].set_ylim(1e-6, 1e2)

# Shared redshift axis on top
ax_top  = axes[0].twiny()
z_ticks = [10, 6, 4, 2, 1, 0.5, 0]
t_ticks = [age_gyr(1.0 / (1 + z)) if z < 100 else ages_gyr[0] for z in z_ticks]
ax_top.set_xlim(axes[0].get_xlim())
ax_top.set_xticks(t_ticks)
ax_top.set_xticklabels([str(z) for z in z_ticks])
ax_top.set_xlabel('z', fontsize=12)

for ax in axes:
    ax.set_xlim(ages_gyr[0], ages_gyr[-1])

out_base = os.path.join(os.path.dirname(__file__), 'bh_history')
plt.savefig(out_base + '.pdf', bbox_inches='tight', dpi=150)
plt.savefig(out_base + '.png', bbox_inches='tight', dpi=150)
print(f"\nSaved bh_history.pdf and bh_history.png")
plt.show()
