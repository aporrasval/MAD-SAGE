"""
Reads PSO track files and reports the best-fit parameters.
Usage: python best_fit_params.py -t <tracks_dir> -S <space_file>
"""

import argparse
import glob
import os
import sys
import numpy as np

sys.path.insert(0, os.path.normpath(os.path.join(__file__, '..', '..', 'output')))
import analysis


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tracks-dir', required=True, help='Path to PSO tracks directory')
    parser.add_argument('-S', '--space-file', required=True, help='Space file used for the PSO run')
    parser.add_argument('--plot', action='store_true', help='Also plot convergence')
    opts = parser.parse_args()

    space = analysis.load_space(opts.space_file)
    param_names  = [n for n in space['name']]
    param_labels = [n for n in space['plot_label']]
    is_log       = space['is_log']

    # Collect all track files
    fx_files  = sorted(glob.glob(os.path.join(opts.tracks_dir, 'track_*_fx.npy')))
    pos_files = sorted(glob.glob(os.path.join(opts.tracks_dir, 'track_*_pos.npy')))

    if not fx_files:
        sys.exit('No track files found in %s' % opts.tracks_dir)

    all_fx  = np.concatenate([np.load(f) for f in fx_files])   # (N_iter * S,)
    all_pos = np.vstack([np.load(f) for f in pos_files])        # (N_iter * S, D)

    # Filter out NaN / inf
    valid = np.isfinite(all_fx)
    if not np.any(valid):
        sys.exit('All fitness values are NaN/inf — PSO did not converge.')

    best_idx = np.argmin(all_fx[valid])
    best_fx  = all_fx[valid][best_idx]
    best_pos = all_pos[valid][best_idx]

    # Per-iteration best for convergence
    n_iter = len(fx_files)
    S = np.load(fx_files[0]).shape[0]
    running_best = np.full(n_iter, np.inf)
    for i, f in enumerate(fx_files):
        fx_i = np.load(f)
        valid_i = np.isfinite(fx_i)
        if np.any(valid_i):
            running_best[i] = np.min(fx_i[valid_i])
        if i > 0:
            running_best[i] = min(running_best[i], running_best[i-1])

    print('\n=== PSO Results ===')
    print('Track files found : %d iterations x %d particles = %d evaluations' % (n_iter, S, n_iter * S))
    print('Valid evaluations : %d / %d' % (valid.sum(), len(valid)))
    print('Best chi2         : %.6f' % best_fx)
    print('\nBest-fit parameters:')
    print('-' * 40)
    for name, label, val in zip(param_names, param_labels, best_pos):
        print('  %-20s = %g  (%s)' % (name, val, label))
    print('-' * 40)

    print('\nPer-iteration best chi2:')
    for i, bv in enumerate(running_best):
        marker = ' <-- best' if (bv == best_fx and i == np.argmin(running_best)) else ''
        print('  iter %02d: %.6f%s' % (i, bv, marker))

    if opts.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, len(param_names) + 1, figsize=(5*(len(param_names)+1), 4))

        iters = np.arange(n_iter)
        ax = axes[0]
        ax.semilogy(iters, running_best, 'k-o')
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Best chi2 (running)')
        ax.set_title('PSO convergence')

        for d, (name, label) in enumerate(zip(param_names, param_labels)):
            ax = axes[d + 1]
            for i, pf in enumerate(pos_files):
                pos_i = np.load(pf)
                fx_i  = np.load(fx_files[i])
                sc = ax.scatter([i]*len(pos_i), pos_i[:, d], c=np.log10(np.where(np.isfinite(fx_i), fx_i, np.nanmax(fx_i[np.isfinite(fx_i)]))), cmap='viridis', s=40)
            ax.axhline(best_pos[d], color='r', lw=1.5, linestyle='--', label='best fit')
            ax.set_xlabel('Iteration')
            ax.set_ylabel(label)
            ax.set_title(name)
            ax.legend(fontsize=8)

        plt.tight_layout()
        out = os.path.join(opts.tracks_dir, 'pso_convergence.pdf')
        plt.savefig(out, bbox_inches='tight')
        print('\nConvergence plot saved to: %s' % out)


if __name__ == '__main__':
    main()
