#
#    ICRAR - International Centre for Radio Astronomy Research
#    (c) UWA - The University of Western Australia, 2019
#    Copyright by UWA (in the framework of the ICRAR)
#    All rights reserved
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#    MA 02111-1307  USA
#

import glob
import logging
import multiprocessing
import os
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np

import analysis
import common


logger = logging.getLogger(__name__)

def count_jobs(job_name):
    """Returns how many jobs with self.jobs_name are currently queued or running"""

    try:
        out, err, code = common.exec_command("squeue")
    except OSError:
        raise RuntimeError("Couldn't run squeue, is it installed?")

    if code:
        raise RuntimeError("squeue failed with code %d: stdout: %s, stderr: %s" % (code, out, err))

    lines_with_jobname = [l for l in out.splitlines() if job_name in l]
    return len(lines_with_jobname)

def _exec_shark(msg, cmdline, cwd=None):
    logger.info('%s with command line: %s', msg, subprocess.list2cmdline(cmdline))
    out, err, code = common.exec_command(cmdline, cwd=cwd)
    if code != 0:
        logger.error('Error while executing %s (exit code %d):\n' +
                     'stdout:\n%s\nstderr:\n%s', cmdline[0], code,
                     common.b2s(out), common.b2s(err))
        raise RuntimeError('%s error' % cmdline[0])


def _to_shark_options(particle, space):
    """Given `particle` in `space` return an iterable with the corresponding
    shark options settings corresponding to that particle"""
    for value, name, is_log in zip(particle, space['name'], space['is_log']):
        if is_log:
            value = 10 ** value
        yield '%s=%s' % (name, value)

def _evaluate(constraint, stat_test, modeldir, subvols):
    y_obs, y_mod, err = constraint.get_data(modeldir, subvols)
    return stat_test(y_obs, y_mod, err)

count = 0
def run_shark_hpc(particles, *args):
    """
    - Handler function for running PSO on Shark on a SLURM based cluster.
    - Swarm size and number of iterations need to be set within the script for now
    - Function needs the relative path to a Shark config file under the -c option
    - For now the subprocess call within must be altered if you are changing shark submit options
    - To find appropriate memory allocations peruse the initial output lines of each particle
    """

    global count

    opts, space, subvols, statTest = args

    # Prepare the file that will be used by the shark submission scripts
    # to determine which values shark will be run for. We put a final \n so the
    # final line gets properly counted by wc (used by shark-submit)
    shark_options = [
        ' '.join(['-o "%s"' % option for option in _to_shark_options(particle, space)])
        for particle in particles
    ]
    positions_fname = tempfile.mktemp('particle_positions.txt')
    logger.info('Creating particle positions file at %s', positions_fname)
    with open(positions_fname, 'wt') as f:
        f.write('\n'.join(shark_options) + '\n')

    # Submit the execution of multiple shark instances, one for each particle
    job_name = 'PSOSMF_%d' % count
    shark_output_base = os.path.join(opts.outdir, job_name)
    cmdline = ['./shark-submit', '-S', opts.shark_binary, '-w', opts.walltime,
               '-n', job_name, '-O', shark_output_base, '-E', positions_fname,
               '-V', ' '.join(map(str, subvols))]
    if opts.account:
        cmdline += ['-a', opts.account]
    if opts.queue:
        cmdline += ['-Q', opts.queue]
    if opts.nodes:
        cmdline += ['-N', str(opts.nodes)]
    else:
        cmdline += ['-m', opts.memory, '-c', str(opts.cpus)]
    cmdline.append(opts.config)
    _exec_shark('Queueing PSO particles', cmdline)

    # Actually wait for the jobs to finish...
    while count_jobs(job_name) > 0:
        time.sleep(10)

    ss = len(particles)
    fx = np.zeros([ss, 3])
    for i in range(ss):
        _, simu, model, _ = common.read_configuration(opts.config)
        particle_outdir = os.path.join(shark_output_base, str(i))
        modeldir = common.get_shark_output_dir(particle_outdir, simu, model)
        fx[i] = [_evaluate(c, statTest, modeldir, subvols) for c in opts.constraints]
        if not opts.keep:
            shutil.rmtree(particle_outdir)

    fx = np.sum(fx, 1)
    logger.info('Particles %r evaluated to %r', particles, fx)

    # this global count just tracks the number of iterations so they can be saved to different files
    count += 1

    return fx

# this is the 'func' that is passed into the pso routine in pso.py
# this doesn't require the new pso.py file, can just run on the original pyswarm function
# only the hpc module above needs the modified pso function, as it has parallelization updates
# the vanilla version only does one particle at a time, whereas the hpc version takes multiple particles at once
def run_shark(particle, *args):

# space is the thing containing the parameter values for the model
    opts, space, subvols, statTest = args

    pid = multiprocessing.current_process().pid
    shark_output_base = os.path.join(opts.outdir, 'output_%d' % pid)
    _, simu, model, _ = common.read_configuration(opts.config)
    modeldir = common.get_shark_output_dir(shark_output_base, simu, model)

# here is where things are executed.  I will need to modify my parameter file here and then run MAD-SAGE
    cmdline = [opts.shark_binary, opts.config,
               '-o', 'execution.output_directory=%s' % shark_output_base,
               '-o', 'execution.simulation_batches=%s' % ' '.join(map(str, subvols))]
    for option in _to_shark_options(particle, space):
        cmdline += ['-o', option]
    _exec_shark('Executing shark instance', cmdline)

    total = sum(_evaluate(c, statTest, modeldir, subvols) for c in opts.constraints)
    logger.info('Particle %r evaluated to %f', particle, total)

    if not opts.keep:
        shutil.rmtree(shark_output_base)

    return total



# number of particles evaluated by this (worker) process; used to give each
# temporary .par file a unique name
_n_evaluations = 0

def run_madsage(particle, *args):

    # space is the thing containing the parameter values for the model
    opts, space, subvols, statTest = args

    # Set simulation info (h, Omega, volume, file list, a_list) in this process.
    # Pool workers are started with 'spawn' on macOS and so do not inherit the
    # configure() call made in main.py.
    import constraints
    constraints.configure(opts.config)

    # Spawned workers also start with no logging handlers, so logger.info/warning
    # calls below (including the 'Particle ... evaluated to' lines that
    # post_analysis.py parses from the captured stdout) would be silently dropped.
    # Send them to stdout in the same format main.py uses. (Under 'fork' the
    # handlers are inherited and this does nothing.)
    if not logging.root.handlers:
        logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                            format='%(asctime)-15s %(name)s#%(funcName)s:%(lineno)s %(message)s')

    # The PSO works in search coordinates (log10 for is_log parameters);
    # convert to the physical values MAD-SAGE expects.
    values = analysis.to_physical(particle, space)

    # create/clear directory for temporary MAD-SAGE output
    spid = str(multiprocessing.current_process().pid)
    _autocalib_dir = os.path.normpath(os.path.abspath(os.path.dirname(__file__)))
    # Per-particle MAD-SAGE output lives inside the run's output directory (-o),
    # e.g. minimillennium/MS_output_<pid>/, so everything from one run stays together
    # (previously autocalibration/DS_output_<pid>/). It is still deleted after evaluation.
    # Renamed DS_output_<pid> -> MS_output_<pid> (MAD-SAGE) on 2026-10-30.
    modeldir = os.path.join(opts.outdir, 'MS_output_' + spid) + '/'
    if not os.path.exists(modeldir): os.makedirs(modeldir)
    # Remove stale galaxy files left by a crashed run with the same PID. (The old
    # subprocess.call(['rm', modeldir + 'model*']) never expanded the wildcard,
    # because no shell was involved, so it silently removed nothing.)
    for stale in glob.glob(modeldir + '*_z*_*'):
        os.remove(stale)

    # copy template parameter file and edit accordingly
    # The temporary .par file is written to the temp/ folder of the run's output
    # directory (-o), e.g. minimillennium/temp/, and kept, so every evaluated particle
    # can be traced to the exact parameter file MAD-SAGE ran with. The name carries
    # the worker PID and a per-worker evaluation counter so files from later
    # iterations do not overwrite earlier ones: <par name>_<pid>_<NNNN>_temp.par
    global _n_evaluations
    _n_evaluations += 1
    config_base = os.path.splitext(os.path.basename(opts.config))[0]
    temp_dir = os.path.join(opts.outdir, 'temp')
    os.makedirs(temp_dir, exist_ok=True)   # also created by main.py; exist_ok avoids a race between workers
    temp_filename = os.path.join(temp_dir, '%s_%s_%04d_temp.par' % (config_base, spid, _n_evaluations))
    #
    f = open(opts.config).readlines()
    #
    # Replace OutputDir and each free parameter. Matching is on the exact first
    # token of the line, so e.g. 'EnergySN' cannot match a longer key that merely
    # starts with it, and commented-out lines ('%Key', ';Key') are never touched.
    # The whole file is scanned (no early break), so OutputDir is always replaced
    # wherever it sits.
    names = list(space['name'])
    found = set()
    for l, line in enumerate(f):
        tokens = line.split()
        if not tokens:
            continue
        key = tokens[0]
        if key == 'OutputDir':
            f[l] = 'OutputDir              ' + modeldir + '\n'
            found.add(key)
        elif key in names:
            # 6 significant figures, e.g. 0.713934 or 8.36414e+50
            f[l] = '%-28s %.6g\n' % (key, values[names.index(key)])
            found.add(key)
    missing = [n for n in ['OutputDir'] + names if n not in found]
    if missing:
        raise RuntimeError('Parameter(s) %s not found in %s; refusing to run with the '
                           'template values instead' % (', '.join(missing), opts.config))
    #
    with open(temp_filename, 'w') as s:
        s.writelines(f)

    # currently assuming that serial in parallel here is the best
    print(temp_filename)
#    cmdline = ['mpirun', '-np', '8', opts.shark_binary, temp_filename]
    repo_root = os.path.normpath(os.path.join(_autocalib_dir, '..'))
    cmdline = [opts.shark_binary, temp_filename]
    _exec_shark('Executing MAD-SAGE instance', cmdline, cwd=repo_root)

    # changed to now take the weighted sum of logs of the (summed) reduced chi^2 of each constraint
    #
    # Objective = 10**( sum_i w_i * log10(stat_i) ) = prod_i stat_i**w_i.
    # Previously the weight sat inside the log, log10(stat_i * w_i), which only
    # multiplied the objective by the constant prod_i w_i and so had no effect on
    # which particle was best. With all weights equal to 1 both forms agree.
    # Both statistics are positive (chi2 >= 0; studentT is a negative
    # log-likelihood > 0), so the log is defined; anything else is reported and
    # scored as +inf so the PSO never treats it as a good fit.
    stats = [np.sum(_evaluate(c, statTest, modeldir, subvols)) for c in opts.constraints]
    if all(np.isfinite(v) and v > 0 for v in stats):
        total = 10**sum(c.weight * np.log10(v) for c, v in zip(opts.constraints, stats))
    else:
        logger.warning('Non-positive or non-finite statistic %r for particle %r; scoring as inf',
                       stats, values)
        total = np.inf
    logger.info('Particle %r evaluated to %f (par file: %s)', values, total, os.path.basename(temp_filename))

    # MAD-SAGE galaxy output is large and is always removed; the small temp .par
    # file is kept in the output directory (see above)
    shutil.rmtree(modeldir)
    return total
