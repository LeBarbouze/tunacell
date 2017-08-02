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
ou = Observable(raw='ou')
gr = Observable(raw='exp_ou_int', differentiate=True, scale='log',
                local_fit=True, time_window=15.)

# define cell-cycle observables
average_gr = Observable(raw='ou', differentiate=False, scale='linear', local_fit=False, mode='average', timing='g')
division_size = Observable(raw='exp_ou_int', differentiate=False, scale='log', local_fit=False, mode='division', timing='g')

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

#    sglplt = UnivariatePlot(univ)
#    sglplt.make_onepoint(mean_show_sd=True, mean_ref=ref_mean)
#    sglplt.make_twopoints(trefs=[40., 80.])
#    sglplt.save(extension='.png')

print('univ objects: OK')

# %% define region(s) for steady state analysis

regs = Regions(parser)
regs.add(label='ALL', tmin=None, tmax=None)

region = regs.get('ALL')

## define computation options

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

#    fig = plot_stationary(stat, fitlog=True, epsilon=.1, ref_decay=ref_decayrate)
#    folderpath = univ.master._get_obs_path()
#    basename = 'stationary_' + region.name
#    fname = os.path.join(folderpath, basename + '.png')
#    fig.savefig(fname)

print('Stationary univ objects: OK')


# %% CROSS CORRELATION

#try:
#    two = initialize_bivariate(*univs)
#    two.import_from_text()
#except BivariateIOError:
#    two = compute_bivariate(*univs)
#    two.export_text()
#
#print('Cross-correlation: OK')

# %% STATIONARY CROSS-CORRELATIONS
couple = univs[2:]
s1, s2 = couple
try:
    stwo = initialize_stationary_bivariate(s1, s2, region, options)
    stwo.import_from_text()
except StationaryBivariateIOError:
    stwo = compute_stationary_bivariate(s1, s2, region, options)
    stwo.export_text()

print('Stationary cross-correlation: OK')

plot_stationary(stwo)