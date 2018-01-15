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

import argparse
import time
import subprocess

import matplotlib.pyplot as plt

from tunacell import Experiment, Parser, Observable, FilterSet
from tunacell.filters.cells import FilterCellIDparity
from tunacell.plotting.samples import SamplePlot


plt.close('all')


# Arguments
argparser = argparse.ArgumentParser()
argparser.add_argument('-t', '--delay', type=float,
                       help='time-delay between two consecutive figures',
                       default=2)
args = argparser.parse_args()


delay = args.delay
long_delay = 4 * delay

# call simurun to create simushort
name = 'simushort'
subprocess.run(['tunasimu', '-f', '-l', name, '--stop', '120.'])

# define the Parser instance, no filter applied
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

## Starting to plot one colony with various options
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
time.sleep(delay)
#plt.close(colplt.fig)


print('* marking samples that verify or not a condition')
colplt.make_plot(length, report_condition=repr(condition),)
colplt.save(user_bname='colony0-even', add_obs=False, with_data_text=False,
            extension='.png')
figs.append(colplt.fig)
figs[-1].show()
time.sleep(delay)


print('* changing color for each cell (colour corresponds to generation index)')
colplt.make_plot(length, report_condition=repr(condition), change_cell_color=True)
colplt.save(user_bname='colony0-even-cell-color', add_obs=False, with_data_text=False,
            extension='.png')
figs.append(colplt.fig)
figs[-1].show()
time.sleep(delay)

print('* changing color for each lineage')
colplt.make_plot(length, report_condition=repr(condition), change_lineage_color=True)
colplt.save(user_bname='colony0-even-lineage-color', add_obs=False, with_data_text=False,
            extension='.png')
figs.append(colplt.fig)
figs[-1].show()
time.sleep(delay)

num_sup = 3
print('* superimposing {} lineages'.format(num_sup))
colplt.make_plot(length, report_condition=repr(condition), change_lineage_color=True,
                 superimpose=num_sup)
colplt.save(user_bname='colony0-even-super3', add_obs=False,
            with_data_text=False, extension='.png')
figs.append(colplt.fig)
figs[-1].show()
#time.sleep(long_delay)

ans = input('Press Enter to close these figures and proceed to next figure batch')
plt.close('all')

figs = []

## We change samples using iterators
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
time.sleep(delay)

print('* Changing to {}'.format(ou.name))
splt.make_plot(ou, report_condition=repr(condition), change_colony_color=True)
splt.save(user_bname='colonies-ou-even', add_obs=False, with_data_text=False,
           extension='.png')
figs.append(splt.fig)
figs[-1].show()
#time.sleep(long_delay)

ans = input('Press Enter to close these figures and proceed to next figure batch')
plt.close('all')


print('* First 5 colonies, superimposing 2 lineages per subplot')
splt = SamplePlot(exp.iter_colonies(size=5), parser=parser,
                   conditions=[condition, ])
splt.make_plot(ou, report_condition=repr(condition), change_colony_color=True,
                superimpose=2)
splt.save(user_bname='colonies5-ou-even', add_obs=False, with_data_text=False,
           extension='.png')
figs.append(splt.fig)
figs[-1].show()
time.sleep(delay)

print('* Superimposing all timeseries on a single subplot')
#splt5 = SamplePlot(ou, parser.iter_colonies(size=5), parser=parser,
#                   conditions=[condition, ])
splt.make_plot(ou, report_condition=repr(condition), change_colony_color=True,
                superimpose='all', show_markers=False, alpha=.5)
splt.save(user_bname='lineages-from-colonies5', add_obs=False,
           with_data_text=False, extension='.png')
figs.append(splt.fig)
figs[-1].show()
#time.sleep(long_delay)

ans = input('Press Enter to close these figures and proceed to next figure batch')
plt.close('all')

print('* Iterating over 10 lineages, changing colour per lineage')
splt = SamplePlot(exp.iter_lineages(size=10, shuffle=True),
                   parser=parser, conditions=[condition, ])
splt.make_plot(ou, report_condition=repr(condition), change_lineage_color=True,
                superimpose='all', alpha=.5, show_markers=False)
splt.save(user_bname='lineages10', add_obs=False, with_data_text=False,
           extension='.png')
figs.append(splt.fig)
figs[-1].show()
time.sleep(delay)

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
time.sleep(delay)


print('* Adding statistic estimates (when they have been computed)')
splt.make_plot(ou, report_condition=repr(condition), change_lineage_color=True,
                superimpose='all', alpha=.5, show_markers=False,
                data_statistics=True)
splt.save(user_bname='lineages10-with-stats', add_obs=False,
           with_data_text=False, extension='.png')
figs.append(splt.fig)
figs[-1].show()
#time.sleep(long_delay)

ans = input('Press enter to close all figures and finish script')
plt.close('all')
#
#print('* A firework summart...')
#for fig in figs:
#    fig.show()
#
#ans = input('Press enter to close all figures and finish script')
