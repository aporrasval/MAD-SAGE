# Adam Stevens, 2022
# Post-analysis of output of PSO with Dark Sage to check fidelity etc
# Updated for local use: reads dump_file from --dump argument, no galprops dependency

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import argparse

import warnings
warnings.filterwarnings("ignore")

fsize = 20
matplotlib.rcParams.update({
    'font.size': fsize, 'xtick.major.size': 10, 'ytick.major.size': 10,
    'xtick.major.width': 1, 'ytick.major.width': 1,
    'ytick.minor.size': 5, 'xtick.minor.size': 5,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'axes.linewidth': 1, 'font.family': 'serif',
    'legend.numpoints': 1, 'legend.columnspacing': 1,
    'legend.fontsize': fsize - 4, 'lines.markeredgewidth': 1.0,
    'errorbar.capsize': 4.0, 'xtick.top': True, 'ytick.right': True
})

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dump', required=True, help='Path to full_dump.txt (stdout from PSO run)')
parser.add_argument('-t', '--tracks-dir', required=True, help='Path to PSO tracks directory')
parser.add_argument('-S', '--space-file', required=True, help='Space file used for the PSO run')
opts = parser.parse_args()

base_dir   = os.path.dirname(opts.dump) + '/'
track_dir  = opts.tracks_dir
dump_file  = opts.dump

# Load space to know parameter names
sys.path.insert(0, os.path.normpath(os.path.join(__file__, '..', '..', 'output')))
import analysis
space       = analysis.load_space(opts.space_file)
param_names = [n for n in space['name']]

# ---- Find best-fit from track files ----------------------------------------
import glob
fx_files  = sorted(glob.glob(os.path.join(track_dir, 'track_*_fx.npy')))
pos_files = sorted(glob.glob(os.path.join(track_dir, 'track_*_pos.npy')))

if not fx_files:
    sys.exit('No track files found in %s' % track_dir)

all_fx  = np.concatenate([np.load(f) for f in fx_files])
all_pos = np.vstack([np.load(f) for f in pos_files])

valid    = np.isfinite(all_fx)
best_idx = np.argmin(all_fx[valid])
best_fx  = all_fx[valid][best_idx]
best_pos = all_pos[valid][best_idx]

print('\n=== Best-fit parameters ===')
for name, val in zip(param_names, best_pos):
    print('  %-20s = %g' % (name, val))
print('  chi2              = %.6f' % best_fx)

# ---- Parse dump file for obs/model data ------------------------------------
SMF_x_obs, SMF_y_obs, SMF_y_mod = [], [], []

f = open(dump_file, 'r')
lines = f.readlines()
on_a_data_line = False
last_line = False
obs = ''
line_start = ''

for l, line in enumerate(lines):

    # identify which constraint the following data belong to
    if 'GAMA_SMF_highres.csv' in line:
        obs = 'SMF'

    # catch the start of a printed data array
    if line[:5] == 'obs x' or line[:5] == 'obs y' or line[:5] == 'mod y':
        line_start = line[:5]
        on_a_data_line = True
        line = line[8:]
        data = np.array([], dtype=np.float32)

    if on_a_data_line:
        if line[-2:] == ']\n':
            last_line = True
            line = line[:-2]
        else:
            last_line = False
            line = line[:-1]

        line_data = np.array(line.split())
        data = np.append(data, line_data.astype(np.float32))

        if last_line:
            if obs == 'SMF' and line_start == 'obs x':
                SMF_x_obs += [data]
            if obs == 'SMF' and line_start == 'obs y':
                SMF_y_obs += [data]
            if obs == 'SMF' and line_start == 'mod y':
                SMF_y_mod += [data]
            on_a_data_line = False
    else:
        last_line = False

print('\nParsed %d SMF evaluations from dump file' % len(SMF_y_mod))

if not SMF_y_mod:
    print('No SMF data found in dump file. Check that stdout was fully captured.')
    sys.exit(0)

# ---- Plot ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 5))

ax.plot(SMF_x_obs[0], SMF_y_obs[0], 'ko', ms=4, label='Driver+2022 (GAMA)', zorder=5)

for i in range(len(SMF_y_mod)):
    alpha = 0.3 if i < len(SMF_y_mod) - 1 else 1.0
    lw    = 0.8 if i < len(SMF_y_mod) - 1 else 1.8
    label = 'DarkSAGE (best fit)' if i == len(SMF_y_mod) - 1 else None
    ax.plot(SMF_x_obs[i], SMF_y_mod[i], '-', alpha=alpha, lw=lw, color='steelblue', label=label)

ax.set_xlim(9, 12.1)
ax.set_ylim(-6, 0)
ax.set_xlabel(r'$\log_{10}(m_*\;[\mathrm{M}_\odot])$', fontsize=fsize)
ax.set_ylabel(r'$\log_{10}(\Phi_*\;[\mathrm{Mpc}^{-3}\,\mathrm{dex}^{-1}])$', fontsize=fsize)
ax.legend()

plt.tight_layout()
out = os.path.join(base_dir, 'reconstructed_SMF.pdf')
plt.savefig(out, bbox_inches='tight')
print('Plot saved to: %s' % out)
