#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
script: plotting-samples.py

This script shows how to plot samples: it first creates a new instance of
simulated experiement, this one shorter than simutest, so that plots do not
extend over too much rows.
"""

from __future__ import print_function
from builtins import input
import sys

import argparse
import time
import subprocess
from tqdm import tqdm

import matplotlib.pyplot as plt

from tunacell import Experiment, Parser, Observable, FilterSet
from tunacell.filters.cells import FilterCellIDparity
from tunacell.plotting.samples import SamplePlot

from tunacell.stats.api import compute_univariate


plt.close('all')


# =============================================================================
# Arguments
# We only set the presentation timing: the delay argument defines the time
# delay between display of two consecutive figures.
# =============================================================================
argparser = argparse.ArgumentParser()

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
       '==             Plotting samples            ==\n'
       '==                                         ==\n'
       '== This tutorial presents an overview of   ==\n'
       '== the basic plotting features of tunacell ==\n'
       '==                                         ==\n'
       '==============tunacell=tutorial==============\n')
print(msg)


# =============================================================================
# We run a numerical simulation of growing and dividing cells with defaults
# parameters except the following
# -f : we force the experiment (to overide previously run simulations)
# -l : we label the experiment 'simushort'
# --stop : we change the time range upper bound to 120 minutes ; this is chosen
# to make colonies with a lower number of leaves (compared to default 180 mins)
# --seed : we chose the seed for the random generator to make comparable
# results
# =============================================================================
name = 'simushort'
if int(sys.version[0]) == 3:
    subprocess.run(['tunasimu', '-f', '-l', name, '--stop', '120.', '--seed', '167389'])
else:
    subprocess.call(['tunasimu', '-f', '-l', name, '--stop', '120.', '--seed', '167389'])

# =============================================================================
# To choose samples from the experiment, we use the Parser class that takes
# an Experiment instance as input. It provides methods to add random samples
# to a buffer of samples that one wants to look at.
# Then we define a filter on the parity of cell identifiers, in order to
# test visually how conditions are handled by tunacell.
# And we define two observables that will be plotted: length and growth rate.
# =============================================================================
path_to_exp = '~/tmptunacell/' + name
exp = Experiment(path_to_exp)
parser = Parser(exp)
# add 10 random samples
parser.add_sample(10)

print('For this example we added 10 random samples')
print(parser)
# define a condition
even = FilterCellIDparity('even')
condition = FilterSet(filtercell=even)

# define observable
length = Observable(name='length', raw='exp_ou_int')
ou = Observable(name='growth-rate', raw='ou')

# =============================================================================
# Now that the parser is defined and that some samples are stored in his buffer
# we will call the plotting functions on particular structures, colonies and
# lineages, from which time-series will be extracted by tunacell.
# We show on following subplots the various options that the plotting module
# provides to help identify patterns of the dynamics.
# =============================================================================
# get one colony
colony = parser.get_colony(0)

figs = []

print()
msg = ('Plotting timeseries from a given colony\n'
       '---------------------------------------')
print(msg)

colplt = SamplePlot([colony, ], parser=parser, conditions=[condition, ])
print('* default settings')
colplt.make_plot(length)
colplt.save(user_bname='colony0', add_obs=False, with_data_text=False,
            extension='.png')
figs.append(colplt.fig)
figs[-1].show()

print('* marking samples that verify or not a condition')
colplt.make_plot(length, report_condition=repr(condition),)
colplt.save(user_bname='colony0-even', add_obs=False, with_data_text=False,
            extension='.png')
figs.append(colplt.fig)
figs[-1].show()

print('* changing color for each cell (colour corresponds to generation index)')
colplt.make_plot(length, report_condition=repr(condition), change_cell_color=True)
colplt.save(user_bname='colony0-even-cell-color', add_obs=False, with_data_text=False,
            extension='.png')
figs.append(colplt.fig)
figs[-1].show()

print('* changing color for each lineage')
colplt.make_plot(length, report_condition=repr(condition), change_lineage_color=True)
colplt.save(user_bname='colony0-even-lineage-color', add_obs=False, with_data_text=False,
            extension='.png')
figs.append(colplt.fig)
figs[-1].show()

num_sup = 3
print('* superimposing {} lineages'.format(num_sup))
colplt.make_plot(length, report_condition=repr(condition), change_lineage_color=True,
                 superimpose=num_sup)
colplt.save(user_bname='colony0-even-super3', add_obs=False,
            with_data_text=False, extension='.png')
figs.append(colplt.fig)
figs[-1].show()

if args.interactive:
    ans = input('Press Enter to close these figures and proceed to next figure batch')
else:
    for seconds in tqdm(range(10*len(figs)), desc='waiting'):
        time.sleep(single_plot_timing/10)
plt.close('all')

figs = []

# We change samples using iterators
plt.close('all')
print()
msg = ('Plotting timeseries obtained by iterating through samples\n'
       '---------------------------------------------------------')
print(msg)
splt = SamplePlot(parser.iter_colonies(mode='samples', size=2),
                   parser=parser, conditions=[condition, ])
print('* Iterating over 2 colonies from our samples, changing colour for each colony')
splt.make_plot(length, report_condition=repr(condition), change_colony_color=True)
splt.save(user_bname='colonies-even', add_obs=False, with_data_text=False,
           extension='.png')
figs.append(splt.fig)
figs[-1].show()

print('* Changing to {}'.format(ou.name))
splt.make_plot(ou, report_condition=repr(condition), change_colony_color=True)
splt.save(user_bname='colonies-ou-even', add_obs=False, with_data_text=False,
           extension='.png')
figs.append(splt.fig)
figs[-1].show()

if args.interactive:
    ans = input('Press Enter to close these figures and proceed to next figure batch')
else:
    for seconds in tqdm(range(10*len(figs)), desc='waiting'):
        time.sleep(single_plot_timing/10)
plt.close('all')
figs = []

print('* First 5 colonies, superimposing 2 lineages per subplot')
splt = SamplePlot(exp.iter_colonies(size=5), parser=parser,
                   conditions=[condition, ])
splt.make_plot(ou, report_condition=repr(condition), change_colony_color=True,
                superimpose=2)
splt.save(user_bname='colonies5-ou-even', add_obs=False, with_data_text=False,
           extension='.png')
figs.append(splt.fig)
figs[-1].show()

print('* Superimposing all timeseries on a single subplot')
#splt5 = SamplePlot(ou, parser.iter_colonies(size=5), parser=parser,
#                   conditions=[condition, ])
splt.make_plot(ou, report_condition=repr(condition), change_colony_color=True,
                superimpose='all', show_markers=False, alpha=.5)
splt.save(user_bname='lineages-from-colonies5', add_obs=False,
           with_data_text=False, extension='.png')
figs.append(splt.fig)
figs[-1].show()


if args.interactive:
    ans = input('Press Enter to close these figures and proceed to next figure batch')
else:
    for seconds in tqdm(range(10*len(figs)), desc='waiting'):
        time.sleep(single_plot_timing/10)
plt.close('all')
figs = []

print('* Iterating over 10 lineages, changing colour per lineage')
splt = SamplePlot(exp.iter_lineages(size=10, shuffle=True),
                   parser=parser, conditions=[condition, ])
splt.make_plot(ou, report_condition=repr(condition), change_lineage_color=True,
                superimpose='all', alpha=.5, show_markers=False)
splt.save(user_bname='lineages10', add_obs=False, with_data_text=False,
           extension='.png')
figs.append(splt.fig)
figs[-1].show()

print('* Adding reference values for average, variance')
# metadata
md = parser.experiment.metadata
# ou expectation values
params = md['ornstein_uhlenbeck_params']
ref_mean = params['target']
ref_var = params['noise']/(2 * params['spring'])

splt.make_plot(ou, report_condition=repr(condition), change_lineage_color=True,
                superimpose='all', alpha=.5, show_markers=False,
                ref_mean=ref_mean, ref_var=ref_var)
splt.save(user_bname='lineages10-with-ref', add_obs=False,
           with_data_text=False, extension='.png')
figs.append(splt.fig)
figs[-1].show()

# we compute univariate statistics so that we can include results in plot
univariate = compute_univariate(exp, ou, cset=[condition, ])
univariate.export_text()
# (this will be explained in univariate-analysis.py)

print('* Adding statistic estimates (when they have been computed)')
splt.make_plot(ou, report_condition=repr(condition), change_lineage_color=True,
                superimpose='all', alpha=.5, show_markers=False,
                data_statistics=True)
splt.save(user_bname='lineages10-with-stats', add_obs=False,
           with_data_text=False, extension='.png')
figs.append(splt.fig)
figs[-1].show()
if args.interactive:
    ans = input('Press enter to close all figures and finish script')
else:
    for seconds in tqdm(range(10*len(figs)), desc='waiting'):
        time.sleep(single_plot_timing/10)
plt.close('all')
