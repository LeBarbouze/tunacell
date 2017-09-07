#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
script: univariate-analysis-2.py

followinf univariate-analysis.py, used to explore more features of univariate
analysis

In this script we show how to load/compute univariate analysis on other
types of observables, and how to build good-looking autocorrelation functions
by using the stationary hypothesis, where only time differences matter, which
in practice increases the effective sample size.
"""

from __future__ import print_function

import sys

import matplotlib.pyplot as plt

from tuna import Experiment, Observable, FilterSet
from tuna.base.observable import FunctionalObservable
from tuna.filters.cells import FilterCellIDparity

from tuna.stats.api import (compute_univariate, load_univariate,
                            compute_stationary, load_stationary, NoValidTimes)
from tuna.stats.single import UnivariateIOError, StationaryUnivariateIOError
from tuna.stats.utils import Regions, CompuParams
from tuna.plotting.dynamics import plot_onepoint, plot_twopoints, plot_stationary


def press_enter(figs):
    """Convenient function to print figures and press Enter"""
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
    return
    
# close all open plots
plt.close('all')

# =============================================================================
# We start with the same settings as in univariate-analysis.py.
# We first load the univariate analysis that we performed, and exported
# in univariate-anbalysis.py (please run that script before starting this one)
# =============================================================================

# define the exp instance, no filter applied
path_to_exp = '~/tmptuna/simutest'
exp = Experiment(path_to_exp)
# define a condition
even = FilterCellIDparity('even')
condition = FilterSet(label='evenID', filtercell=even)

ou = Observable(name='exact-growth-rate', raw='ou')

# Reference values
md = exp.metadata.loc[exp.label]
ref_mean = md.target
ref_var = md.noise / (2*md.spring)
ref_decayrate = md.spring

# loading univariate analysis for the ou observable
univariate = load_univariate(exp, ou, cset=[condition, ])

# =============================================================================
# The last command would raise an exception of type UnivariateIOError if
# computations had not been exported before. We can use this property
# to try loading results, and if it fails, start the computation.
# Below we do so for a few other observables
# =============================================================================

# TIME-LAPSE OBSERVABLES (time-series per cell)

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
figs = []
print('Computing dynamic univariate statistics...')
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

    fig = plot_onepoint(univ, show_ci=True, save=True, **kwargs)
    figs.append(fig)
    fig2 = plot_twopoints(univ, save=True, **kwargs2)
    figs.append(fig2)

press_enter(figs)


# =============================================================================
# A look at the onepoint functions allows the user to identify regions of time
# where the process looks stationary. There is a function to define such
# regions, and in fact, we already use one in the previous computations,
# defined by default: the region 'ALL' that comprises all time values.
# Here we mention the used region explicitely, and we stick to the 'ALL'
# region since we find the process stationary on the entire time course
# =============================================================================

regions = Regions(exp)
regions.reset()  # eliminate all regions except 'ALL'
steady_region = regions.get('ALL')

# and we need to use some computation options (more on that elsewhere)
# define computation options
options = CompuParams()  # leaving to default is safe

# =============================================================================
# Now we proceed in the same way: try to load, if it fails, compute.
# We call the plotting function accordingly.
# =============================================================================
print('Computing stationary univariate statistics...')
figs = []
for obs in continuous_obs + cycle_obs:
    print('{} ...'.format(obs.name))
    # need the univariate object to compute stationary statistics
    univ = univariates_store[obs]
    try:
        stat = load_stationary(univ, steady_region, options)
    except StationaryUnivariateIOError:
        try:
            stat = compute_stationary(univ, steady_region, options)
            stat.export_text()  # save as text files
        except NoValidTimes:
            stat = None
    print('Ok')
    # plotting features
    if obs in [ou, gr, ou2]:
        kwargs = {'ref_decay': ref_decayrate, 'fitlog': True,
                  'fitting_time_max': 45.}  # must be larger than expected decay time
    else:
        kwargs = {}

    if stat is not None:
        fig = plot_stationary(stat, save=True, **kwargs)
        figs.append(fig)

press_enter(figs)


# =============================================================================
# For the sake of demonstration, we define here another, smaller region
# =============================================================================

regions.add(name='beginning', tmin=0., tmax=100.)
reg = regions.get('beginning')

# we just start the computation for exact-growth-rate
univ = univariates_store[ou]
try:
    stat = load_stationary(univ, reg, options)
except StationaryUnivariateIOError:
    stat = compute_stationary(univ, reg, options)
    stat.export_text()

fig = plot_stationary(stat, save=True, ref_decay=ref_decayrate, fitlog=True)
press_enter([fig, ])