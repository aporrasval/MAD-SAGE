#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2019
# Copyright by UWA (in the framework of the ICRAR)
#
# Originally contributed by Mawson Sammons
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
"""
Contains the methods used to quantify the "goodness" of fit
between model and observed parameters
"""

import functools
import math
import sys

import numpy as np
import scipy.special


if sys.version_info[0] == 3:
    b2s = lambda b: b.decode('ascii')
else:
    b2s = lambda b: b

def load_space(fname):
    """Loads the search space from a csv file"""
    space = np.genfromtxt(fname, delimiter=',',
              dtype=[('name', 'object'), ('plot_label', 'object'),
                     ('is_log', 'int'),
                     ('lb', 'float'), ('ub', 'float')])
    space['name'] = [b2s(n) for n in space['name']]
    space['plot_label'] = [b2s(n) for n in space['plot_label']]
    return space

# ---------------------------------------------------------------------------
# Search-space transforms.
#
# The space file gives lb/ub in physical units for every parameter. When
# is_log == 1 the PSO should explore that parameter uniformly in log10, so the
# swarm works in "search coordinates": log10(value) for is_log parameters and
# the value itself otherwise. Everything written to disk (the temporary .par
# file, the track files, the log) is converted back to physical units.
# ---------------------------------------------------------------------------

def search_bounds(space):
    """Lower/upper bounds in search coordinates (log10 for is_log parameters)"""
    is_log = space['is_log'].astype(bool)
    lb = np.array(space['lb'], dtype=float)
    ub = np.array(space['ub'], dtype=float)
    if np.any((lb[is_log] <= 0) | (ub[is_log] <= 0)):
        raise ValueError('Log-sampled parameters must have positive bounds: %r'
                         % list(np.array(space['name'])[is_log]))
    lb[is_log] = np.log10(lb[is_log])
    ub[is_log] = np.log10(ub[is_log])
    return lb, ub

def to_physical(x, space):
    """Convert particle position(s) from search coordinates to physical values.
    Works on a single particle (D,) or a whole swarm (S, D)."""
    x = np.array(x, dtype=float)
    is_log = space['is_log'].astype(bool)
    x[..., is_log] = 10.0 ** x[..., is_log]
    return x

def npsum(f):
    """Return the array-wise sum of the returned array"""
    @functools.wraps(f)
    def wrapper(*args):
        return np.sum(f(*args))
    return wrapper

@npsum
def chi2(obs, mod, err):
    return (mod - obs) ** 2 / (err ** 2) / (len(mod) - 2.0) # hard-coded for 2 parameters being fit here for a reduced chi-squared.  Note that this array actually has to be summed to get the reduced chi^2 of the fit

@npsum
def studentT(obs, mod, err):
    """Negative log-likelihood of the normalised residuals under a Student-t
    distribution (summed over points by @npsum). Smaller is better, so the PSO
    can minimise it directly, and it is always positive (a t pdf never exceeds
    1/sqrt(2*pi)), so it can be combined across constraints in log space.

    nu is set by matching the t-distribution variance, nu/(nu-2), to the mean
    squared normalised residual `var`, i.e. nu = 2*var/(var-1). This only exists
    for var > 1; for var <= 1 the residuals are no broader than Gaussian, so the
    Gaussian (nu -> infinity) limit is used instead.

    Fixes relative to the previous version (which made every PSO particle NaN):
      * returned +log(t), a log-likelihood that is negative and must be
        maximised, while the PSO minimises and execution.py takes log10 of it;
      * used x = (mod-obs)^2/err instead of the squared normalised residual
        (mod-obs)^2/err^2;
      * scipy.special.gamma overflows for large nu; gammaln is used instead.
    """
    sigma = (obs - mod) / err
    var = sum(sigma ** 2 / len(sigma))
    x = sigma ** 2
    if var > 1:
        nu = (2 * var) / (var - 1)
        log_t = (
            scipy.special.gammaln((nu + 1) / 2.0)
            - scipy.special.gammaln(nu / 2.0)
            - 0.5 * np.log(nu * math.pi)
            - (nu + 1) / 2.0 * np.log1p(x / nu)
        )
    else:
        log_t = -0.5 * x - 0.5 * np.log(2 * math.pi)
    return -log_t

stat_tests = {
    'student-t': studentT,
    'chi2': chi2
}
