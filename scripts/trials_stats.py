#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
dynamics.py
"""

from __future__ import print_function

import os
import matplotlib.pyplot as plt

from tuna import Parser, Observable, FilterSet
from tuna.filters.cells import FilterCellIDparity
from tuna.plotting.dynamics import UnivariatePlot, plot_stationary
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
ou = Observable(label='exact-growth-rate', raw='ou')
gr = Observable(label='approx-growth-rate', raw='exp_ou_int',
                differentiate=True, scale='log',
                local_fit=True, time_window=15.)

# define cell-cycle observables
average_gr = Observable(label='averate-growth-rate', raw='ou',
                        differentiate=False, scale='linear',
                        local_fit=False, mode='average', timing='g')
division_size = Observable(label='division-size', raw='exp_ou_int',
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
for obs in [ou, gr, average_gr, division_size]:

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

    uplt = UnivariatePlot(univ)
    uplt.make_onepoint(show_ci=True, mean_ref=ref_mean)
    if obs.timing != 'g':
        uplt.make_twopoints(trefs=[40., 80., 120.])
    else:
        uplt.make_twopoints(trefs=[0, 1, 2])
    uplt.save(extension='.png')

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

    fig = plot_stationary(stat, fitlog=True, epsilon=.1, ref_decay=ref_decayrate)
    folderpath = univ.master._get_obs_path()
    basename = 'stationary_' + stat.region.name
    fname = os.path.join(folderpath, basename + '.png')
    fig.savefig(fname)

print('Stationary univ objects: OK')


# %% CROSS CORRELATION
# unpack univariate objects
u_ou, u_gr, u_avou, u_lendiv = univs

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

    fig = plot_stationary(sbiv)
    folderpath = sbiv.master._get_obs_path()
    basename = 'stationary_bivariate_' + stat.region.name
    fname = os.path.join(folderpath, basename + '.png')
    fig.savefig(fname)
print('Stationary cross-correlation: OK')