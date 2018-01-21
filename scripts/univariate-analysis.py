#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
script: univariate-analysis.py

used to demonstrate how to start analysing time-lapse dynamic data by
computing one-point, and two-point functions for single observables

To do so we will use numerically simulated data (generated with the tunasimu
script) in the simutest folder. No filterset is applied (we consider all data
is valid for statistics). A toy condition is defined by taking the ensemble of
cells with an even identifier. We define a few observables to illustrate tunacell's
capabilities.
"""

from __future__ import print_function

import argparse
from builtins import input
import time

from tqdm import tqdm

import matplotlib.pyplot as plt

from tunacell import Experiment, Observable, FilterSet
from tunacell.filters.cells import FilterCellIDparity
from tunacell.stats.api import compute_univariate
from tunacell.plotting.dynamics import plot_onepoint, plot_twopoints
from tunacell.io import text


# close all open plots
plt.close('all')

# =============================================================================
# First steps is to define which experiment we're analyzing, with which
# filterset and which conditions
# =============================================================================

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
                      default=5)

args = argparser.parse_args()
single_plot_timing = args.time

msg = ('==============tunacell=tutorial==============\n'
       '==                                         ==\n'
       '==        Univariate analysis (1/2)        ==\n'
       '==                                         ==\n'
       '== This tutorial presents an overview of   ==\n'
       '== the univariate analysis (statistics of  ==\n'
       '== single, dynamic observable):            ==\n'
       '==   * computation of moments              ==\n'
       '==   * exploration of Univariate object    ==\n'
       '==   * plotting dynamic moments            ==\n'
       '==  (refer to comments in code to get more ==\n'
       '==   details)                              ==\n'
       '==                                         ==\n'
       '==============tunacell=tutorial==============\n')
print(msg)
print()
# define the Experiment instance without filtering
path_to_exp = args.experiment
exp = Experiment(path_to_exp)
# define a condition
even = FilterCellIDparity('even')
condition = FilterSet(label='evenID', filtercell=even)

# PART I: discovering the analysis on a simple, well defined observable

# =============================================================================
# We start the analysis on a first, simple observable, the 'ou' output from
# numerical simulations, that we use as a model for the cells' exact growth rate
# =============================================================================

ou = Observable(name='exact-growth-rate', raw='ou')

# call the compute_univariate function
msg = 'Launching computation of statistics of the dynamics for {}...'.format(ou.name)
dashes = len(msg) * '*'
print(msg + '\n' + dashes)
univariate = compute_univariate(exp, ou, cset=[condition, ])
# note that we left other input parameters to default...
print('Done')
print()

# =============================================================================
# Exploring the univariate object
# The univariate object stores results for each condition. One item is called
# 'master' as no condition has been applied (use the entire ensemble);
# other items are called as the repr of the FilterSet used for defining the
# condition, here repr(condition).
# Then each item stores a few arrays that corresponds to analysis results
# =============================================================================
msg = 'Inspecting univariate object items..'
dashes = len(msg) * '*'
print(msg + '\n' + dashes)
n_lines = 10
for item_name in ['master', repr(condition)]:
    item = univariate[item_name]
    if item_name == 'master':
        msg = 'Looking at master item (all samples):'
    else:
        msg = 'Looking at {} item (conditioned samples):'.format(item.applied_filter.label)
    dashes = len(msg) * '-'
    print(msg + '\n' + dashes)
    print('One-point function table excerpt:')
    item.display_onepoint(n_lines)
    print()
    print('Two-point function table excerpt:')
    print('(all couple of time-points are present in the entire table)')
    item.display_twopoint(n_lines)
    print()

# =============================================================================
# Now we call plotting functions. For this exact process we know some values
# from theory: we'll import simulation parameters from metadata to show
# the theoretical mean-value of the process, as well as its theoretical
# variance. This way, we can plot exact values together with numerical
# estimates.
# =============================================================================

# Reference values
md = exp.metadata
params = md['ornstein_uhlenbeck_params']
ref_mean = params['target']
ref_var = params['noise']/(2 * params['spring'])
ref_decayrate = params['spring']

# Plotting one-point functions
fig = plot_onepoint(univariate, mean_ref=ref_mean, var_ref=ref_var,
                    show_ci=True, verbose=True)
fig.show()
if args.interactive:
    ans = input('Press Enter to proceed')
else:
    for seconds in tqdm(range(10), desc='waiting'):
        time.sleep(single_plot_timing/10)

# =============================================================================
# Okay the figure fig is divided in 3 scopes:
# * top : number of samples vs time
# * center : average value vs time (here with confidence interval, and ref value)
# * bottom : variance vs time (here with ref value)
# The master curve is shown, as well as each conditioned item (here only the
# evenID condition; observe that number of samples is roughly half of master's
# which is what we expect; we also expect no difference in average value, or
# variance)
# Now let's plot two point functions. There are as many 'sampled' two point
# functions as there are evaluation times (check len(univariate.eval_times));
# we would not see much if we were to plot all of them. We choose to plot
# only three examples, for three time of references.
# =============================================================================

# Plotting two-point functions
fig2 = plot_twopoints(univariate, condition_label='master', trefs=[40., 80., 150.],
                      show_exp_decay=ref_decayrate, verbose=True)
fig2.show()
if args.interactive:
    ans = input('Press Enter to proceed')
else:
    for seconds in tqdm(range(10), desc='waiting'):
        time.sleep(single_plot_timing/10)
# =============================================================================
# Again the figure fig2 is divided in 3 scopes:
# * top : number of samples (i.e. number of independent lineages going from
#   tref to running t; it is between 10^2 to 10^3, quite low to get solid
#   estimates)
# * center : autocorrelation functions, one for each tref (very noisy)
# * bottom : same functions but translated times t-tref
# With the exponential decay guide (we know from theory that autocorrelation
# function decays exponentially for the OU process), it looks plausible
# that our estimates are 'correct', only too low sampling. To gain confidence
# we can turn to the stationary analysis (other script)
# =============================================================================

# =============================================================================
# Saving our computations: if we judge that these results should be stored,
# we can export them as text files. This way, results can be loaded, avoiding
# to perform computations every time.
# =============================================================================

univariate.export_text()  # save results as text files in structured folders

# =============================================================================
# You can go check your ~/tmptunacell/simutest folder. There should be an analysis
# folder with subfolders: filterset folder > observable folder > condition folders
# where results are stored.
# There are also a few sets of functions defined in the tunacell.io.text module
# that allow the user to inspect the various folders in a console,
# and to load objects (at the exception of FunctionalObservable instances)
# =============================================================================

text.print_filtersets(exp)
text.print_observables(exp, exp.fset)
text.print_conditions(exp, exp.fset, ou)

if args.interactive:
    ans = input('Press Enter to proceed')
else:
    time.sleep(single_plot_timing)
