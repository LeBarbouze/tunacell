#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
script: bivariate-analysis.py

This script follows univariate-analysis-2.py.

In this script we show how to load/compute bivariate analysis, where
cross correlations are computed between couple of observables (of similar
type)
"""

from __future__ import print_function
from builtins import input
import argparse
import time

from tqdm import tqdm

import matplotlib.pyplot as plt
import pandas as pd

from tunacell import Experiment, Observable, FilterSet
from tunacell.base.observable import FunctionalObservable
from tunacell.base.observable import set_observable_list
from tunacell.filters.cells import FilterCellIDparity

from tunacell.stats.api import (compute_univariate, load_univariate,
                            compute_stationary, load_stationary, NoValidTimes,
                            compute_bivariate, load_bivariate,
                            compute_stationary_bivariate, load_stationary_bivariate)
from tunacell.stats.single import UnivariateIOError, StationaryUnivariateIOError
from tunacell.stats.two import BivariateIOError, StationaryBivariateIOError
from tunacell.stats.utils import Regions, CompuParams
from tunacell.plotting.dynamics import plot_stationary

from tunacell.plotting.statics import scatter_plot

# close all open plots
plt.close('all')

# Arguments
argparser = argparse.ArgumentParser()
argparser.add_argument('-e', '--experiment', type=str,
                       help='Path to experiment root folder',
                       default='~/tmptunacell/simutest')
argparser.add_argument('-i', '--interactive',
                       help='Ask user to press Enter between parts',
                       action='store_true')
argparser.add_argument('--time', type=float,
                      help='Time per figure when non-interactive mode is on',
                      default=3)

args = argparser.parse_args()
single_plot_timing = args.time

msg = ('==============tunacell=tutorial==============\n'
       '==                                         ==\n'
       '==            Bivariate analysis           ==\n'
       '==                                         ==\n'
       '== This tutorial shows more details about  ==\n'
       '== the bivariate analysis (statistics of a ==\n'
       '== couple of observables):                 ==\n'
       '==   * import/export of univariate results ==\n'
       '==   * computation of bivariate statistics ==\n'
       '==   * at stationary (cross-correlations)  ==\n'
       '==  (refer to comments in code to get more ==\n'
       '==   details)                              ==\n'
       '==                                         ==\n'
       '==============tunacell=tutorial==============\n')
print(msg)
print()

# =============================================================================
# We start with the same settings as in univariate-analysis-2.py.
# We define the same observables and load their univariate analysis, that
# is needed for further computing
# =============================================================================

# define the Parser instance, no filter applied
path_to_exp = args.experiment
exp = Experiment(path_to_exp)
# define a condition
even = FilterCellIDparity('even')
condition = FilterSet(label='evenID', filtercell=even)

# Reference values
md = exp.metadata
params = md['ornstein_uhlenbeck_params']
ref_mean = params['target']
ref_var = params['noise']/(2 * params['spring'])
ref_decayrate = params['spring']

# TIME-LAPSE OBSERVABLES (time-series per cell)

# exact growth rate (model)
ou = Observable(name='exact-growth-rate', raw='ou')

# local estimate of growth rate by using the differentiation of size measurement
# (the raw column 'exp_ou_int' plays the role of cell size in our simulations)
gr = Observable(name='approx-growth-rate', raw='exp_ou_int',
                differentiate=True, scale='log',
                local_fit=True, time_window=15.)

# dynamic, functional observable: twice the growth rate
ou2 = FunctionalObservable(name='double-growth-rate', f=lambda x : 2 * x, observables=[ou, ])

# time-aligned upon root cell division for size analysis
# fixing tref allows to align timeseries to a common origin; the 'root' option
# means that it will be aligned to each colony root cell division time
size = Observable(name='size', raw='exp_ou_int', tref='root')

continuous_obs = [ou, gr, ou2, size]

# SOME CELL-CYCLE TYPE OBSERVABLES (one value per cell)

# cell-cycle average growth rate
average_gr = Observable(name='average-growth-rate', raw='ou',
                        differentiate=False, scale='linear',
                        local_fit=False, mode='average', timing='g')

# size at cell division
division_size = Observable(name='division-size', raw='exp_ou_int',
                           differentiate=False, scale='log',
                           local_fit=False, mode='division', timing='g')

# increase in cell size timed at division time
increase = Observable(name='added-size', raw='exp_ou_int',
                      mode='net-increase-additive', timing='d')

cycle_obs = [average_gr, division_size, increase]

msg = 'Loading univariate results...'
dashes = len(msg) * '*'
print(msg + '\n' + dashes)

univariates_store = {}
for obs in continuous_obs + cycle_obs:
    print('* {} ...'.format(obs.name))
    try:
        univ = load_univariate(exp, obs, cset=[condition, ])
    except UnivariateIOError:
        univ = compute_univariate(exp, obs, cset=[condition, ])
        univ.export_text()  # save as text files
    # store univariate object in a dic indexed by observable
    univariates_store[obs] = univ

    # some options for plotting functions
    trefs = [40., 80., 150.]
    grefs = [1, 2]
    if obs in [ou, gr]:
        kwargs = {'mean_ref': ref_mean,
                  'var_ref': ref_var}
        kwargs2 = {'show_exp_decay': ref_decayrate,
                   'trefs': trefs}
    elif obs in [ou2, ]:
        kwargs = {'mean_ref': 2 * ref_mean,
                  'var_ref': 4 * ref_var}
        kwargs2 = {'show_exp_decay': ref_decayrate,
                   'trefs': trefs}
    elif obs in [size, increase]:
        kwargs = {}
        kwargs2 = {'trefs': trefs}
    elif obs in [average_gr, ]:
        kwargs = {'mean_ref': ref_mean}
        kwargs2 = {'trefs': grefs}
    else:
        kwargs = {}
        kwargs2 = {'trefs': grefs}
#    print('Ok')

regions = Regions(exp)
# regions.reset()  # eliminate all regions except 'ALL'
steady_region = regions.get('ALL')

# and we need to use some computation options (more on that elsewhere)
# define computation options
options = CompuParams()  # leaving to default is safe

# =============================================================================
# We proceed first to dynamic bivariate analysis, which computes matrices
# of covariances.
# To do so, we need to define couples of observables of the same type, i.e.
# either couple of dynamic observables, either couple of cell-cycle observables
# with both absolute timing or both generation timing.
# =============================================================================

couples = [(ou, gr), (average_gr, division_size)]
# note that ordering with a couple matters

msg = 'Computation of bivariate statistics...'
dashes = len(msg) * '*'
print(msg + '\n' + dashes)

for o1, o2 in couples:
    print('* Couple {} - {} ...'.format(o1.name, o2.name))
    u1 = univariates_store[o1]
    u2 = univariates_store[o2]
    try:
        biv = load_bivariate(u1, u2)
    except BivariateIOError:
        biv = compute_bivariate(u1, u2)
        biv.export_text()
#    print('Ok')
    # export master result as dataframe and look at random rows
    print('Looking at some examples computed for master...')
    df = biv.master.as_dataframe()
    if len(df[df['counts'] > 0]) > 10:
        excerpt = df[df['counts'] > 0].sample(10).sort_index()
    else:
        excerpt = df[df['counts'] > 0]
    print('{}'.format(excerpt))
    print()

# =============================================================================
# Now we move to the more informative cross correlation function at stationarity
# =============================================================================
msg = 'Cross-correlation at stationarity'
dashes = len(msg) * '*'
print(msg + '\n' + dashes)

figs = []
for o1, o2 in couples:
    print('* Couple {} - {} ...'.format(o1.name, o2.name))
    u1 = univariates_store[o1]
    u2 = univariates_store[o2]
    try:
        biv = load_stationary_bivariate(u1, u2, steady_region, options)
    except StationaryBivariateIOError:
        biv = compute_stationary_bivariate(u1, u2, steady_region, options)
        biv.export_text()
#    print('Ok')
    # export master result as dataframe and look at random rows
    print('Looking at some examples computed for master...')
    df = biv.master.as_dataframe()
    if len(df[df['counts'] > 0]) > 10:
        excerpt = df[df['counts'] > 0].sample(10).sort_index()
    else:
        excerpt = df[df['counts'] > 0]
    print('{}'.format(excerpt))
    print()

    if o1 == ou:
        kwargs = {'show_exp_decay': ref_decayrate}
    else:
        kwargs = {}
    fig = plot_stationary(biv, save=True, **kwargs)
    fig.show()
    figs.append(fig)

if args.interactive:
    ans = input('Press Enter to close these figures and proceed')
else:
    for seconds in tqdm(range(10*len(figs)), desc='waiting'):
        time.sleep(single_plot_timing/10)
plt.close('all')

# =============================================================================
# We can also analyse the bivariate analysis at stationarity with
# a scatter plot and associated empirical distributions.
# For instance to study the dependency between division_size and cell cycle
# growth rate:
# =============================================================================
msg = 'Scatter plot of division size vs cell cycle growth rate'
dashes = len(msg) * '*'
print(msg + '\n' + dashes)

biv = load_stationary_bivariate(univariates_store[average_gr],
                                univariates_store[division_size],
                                steady_region, options)
fig, ax0, ax1, ax2, hs = scatter_plot(biv, xsize=6, ysize=6,
                                      use_xname=None,
                                      use_yname=None,
                                      groupby=None,
                                      color_index=2,
                                      xunits=r'min$^{{-1}}$',
                                      yunits=r'$\mu m$')
labels = [h.get_label() for h in hs]
ax1.legend(handles=hs, labels=labels, loc='upper left', bbox_to_anchor=(1, 1))

fig.show()
if args.interactive:
    ans = input('Press Enter to close these figures and terminate script')
else:
    for seconds in tqdm(range(10), desc='waiting'):
        time.sleep(single_plot_timing/10)
plt.close('all')

