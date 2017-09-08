#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
script: plotting-samples.py

This script shows how to plot samples.
"""

from __future__ import print_function


import matplotlib.pyplot as plt

from tuna import Experiment, Parser, Observable, FilterSet
from tuna.filters.cells import FilterCellIDparity
from tuna.plotting.samples import SamplePlot

from tutorial import press_enter, args

plt.close('all')
# define the Parser instance, no filter applied
path_to_exp = args.experiment
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

colplt = SamplePlot(length, colony, parser=parser, conditions=[condition, ])
print('* default settings')
colplt.make_plot()
press_enter(colplt.fig)
colplt.save(user_bname='colony0', add_obs=False, with_data_text=False,
            extension='.png')

print('* marking samples that verify or not a condition')
colplt.make_plot(report_condition=repr(condition),)
press_enter(colplt.fig)
colplt.save(user_bname='colony0-even', add_obs=False, with_data_text=False,
            extension='.png')

print('* changing color for each cell (colour corresponds to generation index)')
colplt.make_plot(report_condition=repr(condition), change_cell_color=True)
press_enter(colplt.fig)
colplt.save(user_bname='colony0-even-cell-color', add_obs=False, with_data_text=False,
            extension='.png')

print('* changing color for each lineage')
colplt.make_plot(report_condition=repr(condition), change_lineage_color=True)
press_enter(colplt.fig)
colplt.save(user_bname='colony0-even-lineage-color', add_obs=False, with_data_text=False,
            extension='.png')

num_sup = 3
print('* superimposing {} lineages'.format(num_sup))
colplt.make_plot(report_condition=repr(condition), change_lineage_color=True,
                 superimpose=num_sup)

press_enter(colplt.fig)
colplt.save(user_bname='colony0-even-super3', add_obs=False,
            with_data_text=False, extension='.png')

# %% We change samples using iterators
print()
msg = ('Plotting timeseries obtained by iterating through samples\n'
       '---------------------------------------------------------')
print(msg)
splt2 = SamplePlot(length, parser.iter_colonies(mode='samples', size=2),
                   parser=parser, conditions=[condition, ])
print('* Iterating over 2 colonies from our samples, changing colour for each colony')
splt2.make_plot(report_condition=repr(condition), change_colony_color=True)
press_enter(splt2.fig)
splt2.save(user_bname='colonies-even', add_obs=False, with_data_text=False,
           extension='.png')

splt3 = SamplePlot(ou, parser.iter_colonies(mode='samples', size=2),
                   parser=parser, conditions=[condition, ])
print('* Changing to {}'.format(ou.name))
splt3.make_plot(report_condition=repr(condition), change_colony_color=True)
press_enter(splt3.fig)
splt3.save(user_bname='colonies-ou-even', add_obs=False, with_data_text=False,
           extension='.png')

print('* First 5 colonies, superimposing 2 colonies per subplot')
splt4 = SamplePlot(ou, parser.iter_colonies(size=5), parser=parser,
                   conditions=[condition, ])
splt4.make_plot(report_condition=repr(condition), change_colony_color=True,
                superimpose=2)
press_enter(splt4.fig)
splt4.save(user_bname='colonies5-ou-even', add_obs=False, with_data_text=False,
           extension='.png')

print('* Superimposing all timeseries on a single subplot')
splt5 = SamplePlot(ou, parser.iter_colonies(size=5), parser=parser,
                   conditions=[condition, ])
splt5.make_plot(report_condition=repr(condition), change_colony_color=True,
                superimpose='all', show_markers=False, alpha=.5)
press_enter(splt5.fig)
splt5.save(user_bname='lineages-from-colonies5', add_obs=False,
           with_data_text=False, extension='.png')

print('* Iterating over 10 lineages, changing colour per lineage')
splt6 = SamplePlot(ou, parser.iter_lineages(size=10, shuffle=True),
                   parser=parser, conditions=[condition, ])
splt6.make_plot(report_condition=repr(condition), change_lineage_color=True,
                superimpose='all', alpha=.5, show_markers=False)
press_enter(splt6.fig)
splt6.save(user_bname='lineages10', add_obs=False, with_data_text=False,
           extension='.png')

print('* Adding reference values for average, variance')
# metadata
md = parser.experiment.metadata.loc[parser.experiment.label]
# ou expectation values
ref_mean = md.target
ref_var = md.noise/(2 * md.spring)

splt6.make_plot(report_condition=repr(condition), change_lineage_color=True,
                superimpose='all', alpha=.5, show_markers=False,
                ref_mean=ref_mean, ref_var=ref_var)
press_enter(splt6.fig)
splt6.save(user_bname='lineages10-with-ref', add_obs=False,
           with_data_text=False, extension='.png')


print('* Adding statistic estimates')
splt6.make_plot(report_condition=repr(condition), change_lineage_color=True,
                superimpose='all', alpha=.5, show_markers=False,
                data_statistics=True)
press_enter(splt6.fig)
splt6.save(user_bname='lineages10-with-stats', add_obs=False,
           with_data_text=False, extension='.png')