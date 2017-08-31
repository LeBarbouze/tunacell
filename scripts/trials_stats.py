#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
dynamics.py
"""

from __future__ import print_function

import os
import matplotlib.pyplot as plt

from tuna import Parser, Observable, FilterSet
from tuna.observable import FunctionalObservable
from tuna.filters.cells import FilterCellIDparity
from tuna.plotting.dynamics import plot_onepoint, plot_twopoints, plot_stationary
# import api functions
from tuna.stats.api import (compute_univariate_dynamics,
                            compute_stationary_univariate,
                            compute_bivariate,
                            compute_stationary_bivariate,
                            initialize_univariate,
                            initialize_stationary_univariate,
                            initialize_bivariate,
                            initialize_stationary_bivariate)
from tuna.stats.single import UnivariateIOError, StationaryUnivariateIOError
from tuna.stats.two import BivariateIOError, StationaryBivariateIOError
from tuna.stats.utils import Regions, CompuParams
plt.close('all')

# define the Parser instance, no filter applied
path_to_exp = '~/tmptuna/simutest'
parser = Parser(path_to_exp)
# define a condition
even = FilterCellIDparity('even')
condition = FilterSet(label='evenID', filtercell=even)

# define dynamic observables
ou = Observable(name='exact-growth-rate', raw='ou')
gr = Observable(name='approx-growth-rate', raw='exp_ou_int',
                differentiate=True, scale='log',
                local_fit=True, time_window=15.)
# dynamic, functional observable
ou2 = FunctionalObservable(name='double-growth-rate', f=lambda x : 2 * x, observables=[ou, ])
# time-aligned upon root cell division for size analysis
size = Observable(name='size', raw='exp_ou_int', tref='root')

# define cell-cycle observables
average_gr = Observable(name='averate-growth-rate', raw='ou',
                        differentiate=False, scale='linear',
                        local_fit=False, mode='average', timing='g')
division_size = Observable(name='division-size', raw='exp_ou_int',
                           differentiate=False, scale='log',
                           local_fit=False, mode='division', timing='g')

# %% loop over both observables
univs = []
stats = []

# %% Reference values
md = parser.experiment.metadata.loc[parser.experiment.label]
ref_mean = md.target
ref_decayrate = md.spring
tmin = md.start
tmax = md.stop
period = md.period
# %% univ OBJECTS
for obs in [ou, gr, average_gr, division_size, size, ou2]:

    # Statistics: if import fails, compute
    try:
        univ = initialize_univariate(parser, obs, cset=[condition, ])
        univ.import_from_text()
    except UnivariateIOError as uerr:
        print('Impossible to load univariate {}'.format(uerr))
        print('Launching computation')
        univ = compute_univariate_dynamics(parser, obs, cset=[condition, ])
        univ.export_text()
    univs.append(univ)

    # make and save plots
    if obs == ou2:
        fig = plot_onepoint(univ, show_ci=True, mean_ref=2*ref_mean, save=True)
    elif obs == division_size:
        fig = plot_onepoint(univ, show_ci=True, save=True)
    elif obs == size:
        fig = plot_onepoint(univ, show_ci=True, save=True)
    else:
        fig = plot_onepoint(univ, show_ci=True, mean_ref=ref_mean, save=True)
    if obs.timing != 'g':
        fig2 = plot_twopoints(univ, trefs=[40., 80., 120.], save=True)
    else:
        fig2 = plot_twopoints(univ, trefs=[0, 1, 2], save=True)

print('univ objects: OK')

# %% define region(s) for steady state analysis

regs = Regions(parser)
regs.add(label='ALL', tmin=None, tmax=None)

region = regs.get('ALL')

# define computation options
options = CompuParams()
# %% STATIONARY univs
for index, obs in enumerate([ou, gr, average_gr, division_size]):
    univ = univs[index]
    try:
        stat = initialize_stationary_univariate(univ, region, options)
        stat.import_from_text()
    except StationaryUnivariateIOError as suerr:
        print('Impossible to load stationary {}'.format(suerr))
        print('Launching computation')
        stat = compute_stationary_univariate(univ, region, options)
        stat.export_text()
    stats.append(stat)

    fig = plot_stationary(stat, fitlog=True, epsilon=.1,
                          ref_decay=ref_decayrate, save=True)

print('Stationary univ objects: OK')


# %% CROSS CORRELATION
# unpack univariate objects
u_ou, u_gr, u_avou, u_lendiv, u_ou2, u_size = univs

bivs = []
for u1, u2 in [(u_ou, u_gr), (u_avou, u_lendiv)]:
    try:
        biv = initialize_bivariate(u1, u2)
        biv.import_from_text()
    except BivariateIOError:
        biv = compute_bivariate(u1, u2)
        biv.export_text()
    bivs.append(biv)
print('Cross-correlation: OK')

# %% STATIONARY CROSS-CORRELATIONS
sbivs = []
for u1, u2 in [(u_ou, u_gr), (u_avou, u_lendiv)]:
    try:
        sbiv = initialize_stationary_bivariate(u1, u2, region, options)
        sbiv.import_from_text()
    except StationaryBivariateIOError:
        sbiv = compute_stationary_bivariate(u1, u2, region, options)
        sbiv.export_text()
        sbivs.append(sbiv)

    fig = plot_stationary(sbiv, save=True)

print('Stationary cross-correlation: OK')