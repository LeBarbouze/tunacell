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

import sys

import matplotlib.pyplot as plt
import pandas as pd

from tuna import Experiment, Observable, FilterSet
from tuna.base.observable import FunctionalObservable
from tuna.filters.cells import FilterCellIDparity

from tuna.stats.api import (compute_univariate, load_univariate,
                            compute_stationary, load_stationary, NoValidTimes,
                            compute_bivariate, load_bivariate,
                            compute_stationary_bivariate, load_stationary_bivariate)
from tuna.stats.single import UnivariateIOError, StationaryUnivariateIOError
from tuna.stats.two import BivariateIOError, StationaryBivariateIOError
from tuna.stats.utils import Regions, CompuParams
from tuna.plotting.dynamics import plot_onepoint, plot_twopoints, plot_stationary

from tuna.io import text


# close all open plots
plt.close('all')

# =============================================================================
# We start with the same settings as in univariate-analysis-2.py.
# We define the same observables and load their univariate analysis, that
# is needed for further computing
# =============================================================================

# define the Parser instance, no filter applied
path_to_exp = '~/tmptuna/simutest'
exp = Experiment(path_to_exp)
# define a condition
even = FilterCellIDparity('even')
condition = FilterSet(label='evenID', filtercell=even)

# Reference values
md = exp.metadata.loc[exp.label]
ref_mean = md.target
ref_var = md.noise / (2*md.spring)
ref_decayrate = md.spring
tmin = md.start
tmax = md.stop
period = md.period


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


univariates_store = {}
for obs in continuous_obs + cycle_obs:
    print('{} ...'.format(obs.name))
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
    print('Ok')

regions = Regions(exp)
regions.reset()  # eliminate all regions except 'ALL'
regions.add(name='steady', tmin=20., tmax=160.)
steady_region = regions.get('steady')

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

for o1, o2 in couples:
    print('Couple {} - {} ...'.format(o1.name, o2.name))
    u1 = univariates_store[o1]
    u2 = univariates_store[o2]
    try:
        biv = load_bivariate(u1, u2)
    except BivariateIOError:
        biv = compute_bivariate(u1, u2)
        biv.export_text()
    print('Ok')
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
print('Cross-correlation at stationarity\n'
      '---------------------------------')
figs = []
for o1, o2 in couples:
    print('Couple {} - {} ...'.format(o1.name, o2.name))
    u1 = univariates_store[o1]
    u2 = univariates_store[o2]
    try:
        biv = load_stationary_bivariate(u1, u2, steady_region, options)
    except StationaryBivariateIOError:
        biv = compute_stationary_bivariate(u1, u2, steady_region, options)
        biv.export_text()
    print('Ok')
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
        kwargs = {'ref_decay': ref_decayrate}
    else:
        kwargs = {}
    fig = plot_stationary(biv, save=True, **kwargs)
    figs.append(fig)

# when run from ipython, figure should automatically be plotted
try:
    __IPYTHON__
# otherwise call .plot() and wait for pressing Enter
except NameError:
    for fig in figs:
        fig.show()
        if sys.version_info[0] == 2:
            ans = raw_input('Press Enter to proceed...')
        else:
            ans = input('Press Enter to proceed...')

# =============================================================================
# Note that we can get easily non-dynamic bivariate analysis
# '(mixing cell-cycle observables with different reporting times)
# =============================================================================
plt.figure()
couple = [division_size, increase]
print('Bivariate sampling of {} and {} ...'.format(couple[0].name, couple[1].name))
all_dfs = []
for li in exp.iter_lineages(size=100):  # testing
    dfs = []
    for obs in couple:
        ts = li.get_timeseries(obs, cset=[condition, ])
        dfs.append(ts.to_dataframe(sharp_tleft=steady_region.tmin,
                                   sharp_tright=steady_region.tmax))
    all_dfs.append(pd.merge(*dfs, how='outer'))
df = pd.concat(all_dfs, ignore_index=True)
excerpt = df.sample(10).sort_index()
print('{}'.format(excerpt))

plt.scatter(df[couple[0].name], df[couple[1].name])
plt.xlabel(couple[0].name)
plt.ylabel(couple[1].name)

# when run from ipython, figure should automatically be plotted
try:
    __IPYTHON__
# otherwise call .plot() and wait for pressing Enter
except NameError:
    plt.show()
    if sys.version_info[0] == 2:
        ans = raw_input('Press Enter to proceed...')
    else:
        ans = input('Press Enter to proceed...')
