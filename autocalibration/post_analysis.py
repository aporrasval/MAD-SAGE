# Adam Stevens, 2022
# Post-analysis of output of PSO with MAD-SAGE to check fidelity etc
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
# -S used to be required; it now defaults to space_BHspin.txt next to this file,
# the same default as main.py and diagnostics.py. Pass the space file of the run if it was different.
parser.add_argument('-S', '--space-file', help='Space file used for the PSO run (default: space_BHspin.txt in the autocalibration directory)',
                    default=os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'space_BHspin.txt')))
# Plots are written to the home directory (/home/ajp79/) by default; previously they went
# next to the dump file. Pass -o to save them elsewhere.
parser.add_argument('-o', '--out-dir', default=None, help='Directory to save plots in (default: plots/ next to the dump file, e.g. minimillennium/plots/)')
opts = parser.parse_args()

# Default to the dump file's own directory (the default used to be the cluster
# home directory /home/ajp79/, which does not exist on other machines)
#
# Default is now the plots/ folder of the run directory (the dump file's
# directory), e.g. minimillennium/plots/, created if needed.
if opts.out_dir is None:
    opts.out_dir = os.path.join(os.path.dirname(os.path.abspath(opts.dump)), 'plots')
    os.makedirs(opts.out_dir, exist_ok=True)
base_dir   = os.path.join(opts.out_dir, '')   # ensure trailing slash
if not os.path.isdir(base_dir):
    sys.exit('Output directory does not exist: %s' % base_dir)
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
# Each constraint is identified by the observation file it loads (printed as
# "Loading observations from <file>" just before its data). File names must match
# those used in constraints.py.
obs_files = {
    'GAMA_SMF_highres.csv'       : 'SMF',
    'Jones2018_HIMF.dat'         : 'HIMF',
    'CSFH_DSILVA+23_Ver_Final.csv': 'CSFRD',
}
constraint_names = ['SMF', 'HIMF', 'CSFRD']

# Per constraint: lists of obs x, obs y, mod y (one entry per evaluation), and
# fx = chi2 of the particle each 'mod y' block belongs to.
# Each evaluation in the dump prints SMF, HIMF and CSFRD blocks, then a
# "Particle ... evaluated to <fx>" log line, so the fx that follows a block is its chi2.
blocks = {c: {'x_obs': [], 'y_obs': [], 'y_mod': [], 'fx': []} for c in constraint_names}
n_pending = {c: 0 for c in constraint_names}   # blocks seen since the last "evaluated to" line (should only ever be 0 or 1)

f = open(dump_file, 'r')
lines = f.readlines()
on_a_data_line = False
last_line = False
obs = ''
line_start = ''

for l, line in enumerate(lines):

    # identify which constraint the following data belong to.
    # Every constraint prints a "Loading observations from <file>" line before its data,
    # so reset the label on each one; otherwise HIMF and CSFRD blocks are mislabelled as SMF.
    if 'Loading observations from' in line:
        obs = 'other'
        for fname, cname in obs_files.items():
            if fname in line:
                obs = cname

    # attach this particle's chi2 to the blocks printed during its evaluation
    if 'evaluated to' in line and 'Particle ' in line:
        fx = float(line.rsplit('evaluated to', 1)[1].split()[0])
        for c in constraint_names:
            if n_pending[c] == 1:
                blocks[c]['fx'][-1] = fx
            elif n_pending[c] > 1:
                print('Warning: %d %s blocks before one "evaluated to" line (line %d); '
                      'parallel output interleaved, chi2 matching may be wrong here' % (n_pending[c], c, l + 1))
            n_pending[c] = 0

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
            if obs in blocks:
                if line_start == 'obs x':
                    blocks[obs]['x_obs'] += [data]
                if line_start == 'obs y':
                    blocks[obs]['y_obs'] += [data]
                if line_start == 'mod y':
                    blocks[obs]['y_mod'] += [data]
                    blocks[obs]['fx'] += [np.nan]   # filled in when this particle's "evaluated to" line is reached
                    n_pending[obs] += 1
            on_a_data_line = False
    else:
        last_line = False

print()
for c in constraint_names:
    print('Parsed %d %s evaluations from dump file' % (len(blocks[c]['y_mod']), c))

# ---- Plot settings per constraint ------------------------------------------
# Axis limits are for display only (they do not affect the fit). SMF limits are the
# original ones from this script; HIMF x-range spans the HIMF(9-11) fit range, and the
# CSFRD x-range spans the CSFRDH(0-13.8) look-back-time range. Units follow constraints.py:
# masses in Msun (h removed), number densities per Mpc^3 (h removed), SFRD in Msun/yr/Mpc^3.
plot_settings = {
    'SMF': {
        'obs_label': 'Driver+2022 (GAMA)',
        'xlim': (9, 12.1), 'ylim': (-6, 0),
        'xlabel': r'$\log_{10}(m_*\;[\mathrm{M}_\odot])$',
        'ylabel': r'$\log_{10}(\Phi_*\;[\mathrm{Mpc}^{-3}\,\mathrm{dex}^{-1}])$',
        'outfile': 'reconstructed_SMF.pdf',
    },
    'HIMF': {
        'obs_label': 'Jones+2018 (ALFALFA)',
        'xlim': (9, 11), 'ylim': (-6, 0),
        'xlabel': r'$\log_{10}(m_{\rm HI}\;[\mathrm{M}_\odot])$',
        'ylabel': r'$\log_{10}(\Phi_{\rm HI}\;[\mathrm{Mpc}^{-3}\,\mathrm{dex}^{-1}])$',
        'outfile': 'reconstructed_HIMF.pdf',
    },
    'CSFRD': {
        'obs_label': "D'Silva+2023",
        'xlim': (0, 14), 'ylim': (-3.5, -0.5),
        'xlabel': r'Look-back time [Gyr]',
        'ylabel': r'$\log_{10}(\rho_{\rm SFR}\;[\mathrm{M}_\odot\,\mathrm{yr}^{-1}\,\mathrm{Mpc}^{-3}])$',
        'outfile': 'reconstructed_CSFRD.pdf',
    },
}

# ---- Plot ------------------------------------------------------------------
for c in constraint_names:
    x_obs, y_obs, y_mod = blocks[c]['x_obs'], blocks[c]['y_obs'], blocks[c]['y_mod']
    ps = plot_settings[c]

    if not y_mod:
        print('No %s data found in dump file. Check that stdout was fully captured.' % c)
        continue

    # Pick the evaluation with the lowest chi2 (rather than simply the last one printed)
    fx = np.array(blocks[c]['fx'])
    if np.all(np.isnan(fx)):
        print('Warning: no chi2 values matched to %s blocks; highlighting the last evaluation instead' % c)
        i_best = len(y_mod) - 1
    else:
        i_best = int(np.nanargmin(fx))
        print('Best %s evaluation: #%d of %d, chi2 = %.6f (track-file best = %.6f)'
              % (c, i_best + 1, len(y_mod), fx[i_best], best_fx))

    fig, ax = plt.subplots(figsize=(6, 5))

    ax.plot(x_obs[0], y_obs[0], 'ko', ms=4, label=ps['obs_label'], zorder=5)

    # all evaluations faint, best-chi2 evaluation drawn on top
    for i in range(len(y_mod)):
        if i == i_best:
            continue
        ax.plot(x_obs[i], y_mod[i], '-', alpha=0.3, lw=0.8, color='steelblue')
    ax.plot(x_obs[i_best], y_mod[i_best], '-', alpha=1.0, lw=1.8, color='steelblue',
            label='MAD-SAGE (best fit)', zorder=6)

    ax.set_xlim(*ps['xlim'])
    ax.set_ylim(*ps['ylim'])
    ax.set_xlabel(ps['xlabel'], fontsize=fsize)
    ax.set_ylabel(ps['ylabel'], fontsize=fsize)
    ax.legend()

    plt.tight_layout()
    out = os.path.join(base_dir, ps['outfile'])
    plt.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print('Plot saved to: %s' % out)
