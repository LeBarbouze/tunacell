#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
script: plotting-samples.py

This script shows how to plot samples: it first creates a new instance of
simulated experiement, this one shorter than simutest, so that plots do not
extend over too much rows.
"""

from __future__ import print_function

import subprocess

import matplotlib.pyplot as plt

from tunacell import Experiment, Parser, Observable, FilterSet
from tunacell.filters.cells import FilterCellIDparity
from tunacell.plotting.samples import SamplePlot

from tutorial import press_enter

plt.close('all')

# call simurun to create simushort
subprocess.run(['simurun.py', '-f', '-l', 'simushort', '--stop', '120.'])

# define the Parser instance, no filter applied
path_to_exp = '~/tmptunacell/simushort'
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

# %% Starting to plot one colony with various options
# get one colony
colony = parser.get_colony(0)

print()
msg = ('Plotting timeseries from a given colony\n'
       '---------------------------------------')
print(msg)

colplt = SamplePlot([colony, ], parser=parser, conditions=[condition, ])
print('* default settings')
colplt.make_plot(length)
press_enter(colplt.fig)
colplt.save(user_bname='colony0', add_obs=False, with_data_text=False,
            extension='.png')

print('* marking samples that verify or not a condition')
colplt.make_plot(length, report_condition=repr(condition),)
press_enter(colplt.fig)
colplt.save(user_bname='colony0-even', add_obs=False, with_data_text=False,
            extension='.png')

print('* changing color for each cell (colour corresponds to generation index)')
colplt.make_plot(length, report_condition=repr(condition), change_cell_color=True)
press_enter(colplt.fig)
colplt.save(user_bname='colony0-even-cell-color', add_obs=False, with_data_text=False,
            extension='.png')

print('* changing color for each lineage')
colplt.make_plot(length, report_condition=repr(condition), change_lineage_color=True)
press_enter(colplt.fig)
colplt.save(user_bname='colony0-even-lineage-color', add_obs=False, with_data_text=False,
            extension='.png')

num_sup = 3
print('* superimposing {} lineages'.format(num_sup))
colplt.make_plot(length, report_condition=repr(condition), change_lineage_color=True,
                 superimpose=num_sup)
press_enter(colplt.fig)
colplt.save(user_bname='colony0-even-super3', add_obs=False,
            with_data_text=False, extension='.png')

# %% We change samples using iterators
plt.close('all')
print()
msg = ('Plotting timeseries obtained by iterating through samples\n'
       '---------------------------------------------------------')
print(msg)
splt = SamplePlot(parser.iter_colonies(mode='samples', size=2),
                   parser=parser, conditions=[condition, ])
print('* Iterating over 2 colonies from our samples, changing colour for each colony')
splt.make_plot(length, report_condition=repr(condition), change_colony_color=True)
press_enter(splt.fig)
splt.save(user_bname='colonies-even', add_obs=False, with_data_text=False,
           extension='.png')

print('* Changing to {}'.format(ou.name))
splt.make_plot(ou, report_condition=repr(condition), change_colony_color=True)
press_enter(splt.fig)
splt.save(user_bname='colonies-ou-even', add_obs=False, with_data_text=False,
           extension='.png')

print('* First 5 colonies, superimposing 2 lineages per subplot')
splt = SamplePlot(parser.iter_colonies(size=5), parser=parser,
                   conditions=[condition, ])
splt.make_plot(ou, report_condition=repr(condition), change_colony_color=True,
                superimpose=2)
press_enter(splt.fig)
splt.save(user_bname='colonies5-ou-even', add_obs=False, with_data_text=False,
           extension='.png')

print('* Superimposing all timeseries on a single subplot')
#splt5 = SamplePlot(ou, parser.iter_colonies(size=5), parser=parser,
#                   conditions=[condition, ])
splt.make_plot(ou, report_condition=repr(condition), change_colony_color=True,
                superimpose='all', show_markers=False, alpha=.5)
press_enter(splt.fig)
splt.save(user_bname='lineages-from-colonies5', add_obs=False,
           with_data_text=False, extension='.png')

print('* Iterating over 10 lineages, changing colour per lineage')
splt = SamplePlot(parser.iter_lineages(size=10, shuffle=True),
                   parser=parser, conditions=[condition, ])
splt.make_plot(ou, report_condition=repr(condition), change_lineage_color=True,
                superimpose='all', alpha=.5, show_markers=False)
press_enter(splt.fig)
splt.save(user_bname='lineages10', add_obs=False, with_data_text=False,
           extension='.png')

print('* Adding reference values for average, variance')
# metadata
md = parser.experiment.metadata.loc[parser.experiment.label]
# ou expectation values
ref_mean = md.target
ref_var = md.noise/(2 * md.spring)

splt.make_plot(ou, report_condition=repr(condition), change_lineage_color=True,
                superimpose='all', alpha=.5, show_markers=False,
                ref_mean=ref_mean, ref_var=ref_var)
press_enter(splt.fig)
splt.save(user_bname='lineages10-with-ref', add_obs=False,
           with_data_text=False, extension='.png')


print('* Adding statistic estimates (when they have been computed)')
splt.make_plot(ou, report_condition=repr(condition), change_lineage_color=True,
                superimpose='all', alpha=.5, show_markers=False,
                data_statistics=True)
press_enter(splt.fig)
splt.save(user_bname='lineages10-with-stats', add_obs=False,
           with_data_text=False, extension='.png')
