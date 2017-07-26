#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
plotting-samples.py
"""

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

from tuna import Parser, Observable, FilterSet
from tuna.filters.cells import FilterCellIDparity
from tuna.plotting.samples import SamplePlot

plt.close('all')
# define the Parser instance, no filter applied
path_to_exp = '~/tmptuna/simushort'
parser = Parser(path_to_exp)
# add 10 random samples
parser.add_sample(10)
# define a condition
even = FilterCellIDparity('even')
condition = FilterSet(filtercell=even)

# define observable
length = Observable(raw='exp_ou_int')
ou = Observable(raw='ou')

# get one colony
colony = parser.get_colony(0)

# %% plot  colony with default settings
colplt = SamplePlot(length, colony, parser=parser, conditions=[condition, ])
colplt.make_plot()
colplt.save(user_bname='colony0', add_obs=False, with_data_text=False,
            extension='.png')

# %% plot colony with condition
colplt.make_plot(report_condition=repr(condition),)
colplt.save(user_bname='colony0-even', add_obs=False, with_data_text=False,
            extension='.png')

# %% change cell color
colplt.make_plot(report_condition=repr(condition), change_cell_color=True)
colplt.save(user_bname='colony0-even-cell-color', add_obs=False, with_data_text=False,
            extension='.png')

# %% change lineage color
colplt.make_plot(report_condition=repr(condition), change_lineage_color=True)
colplt.save(user_bname='colony0-even-lineage-color', add_obs=False, with_data_text=False,
            extension='.png')

# %% superimpose a number of lineages
colplt.make_plot(report_condition=repr(condition), change_lineage_color=True,
                 superimpose=3)
colplt.save(user_bname='colony0-even-super3', add_obs=False,
            with_data_text=False, extension='.png')

# %% plot two contiguous colonies from samples with condition
splt2 = SamplePlot(length, parser.iter_colonies(mode='samples', size=2),
                   parser=parser, conditions=[condition, ])
splt2.make_plot(report_condition=repr(condition), change_colony_color=True)
splt2.save(user_bname='colonies-even', add_obs=False, with_data_text=False,
           extension='.png')

# %% plot two contiguous colonies from samples with condition
splt3 = SamplePlot(ou, parser.iter_colonies(mode='samples', size=2),
                   parser=parser, conditions=[condition, ])
splt3.make_plot(report_condition=repr(condition), change_colony_color=True)
splt3.save(user_bname='colonies-ou-even', add_obs=False, with_data_text=False,
           extension='.png')

# %% plot 5 first colonies
splt4 = SamplePlot(ou, parser.iter_colonies(size=5), parser=parser,
                   conditions=[condition, ])
splt4.make_plot(report_condition=repr(condition), change_colony_color=True,
                superimpose=2)
splt4.save(user_bname='colonies5-ou-even', add_obs=False, with_data_text=False,
           extension='.png')

# %% parse identical lineages, but collapse in one axe
splt5 = SamplePlot(ou, parser.iter_colonies(size=5), parser=parser,
                   conditions=[condition, ])
splt5.make_plot(report_condition=repr(condition),change_colony_color=True,
                superimpose='all', show_markers=False, alpha=.5)
splt5.save(user_bname='lineages-from-colonies5', add_obs=False,
           with_data_text=False, extension='.png')

# %% iterating over lineages directly
splt6 = SamplePlot(ou, parser.iter_lineages(size=10, shuffle=True),
                   parser=parser, conditions=[condition, ])
splt6.make_plot(report_condition=repr(condition), change_lineage_color=True,
                superimpose='all', alpha=.5, show_markers=False)
splt6.save(user_bname='lineages10', add_obs=False, with_data_text=False,
           extension='.png')

# %% adding theoretical values
# metadata
md = parser.experiment.metadata.loc[parser.experiment.label]
# ou expectation values
ref_mean = md.target
ref_var = md.noise/(2 * md.spring)

splt6.make_plot(report_condition=repr(condition), change_lineage_color=True,
                superimpose='all', alpha=.5, show_markers=False,
                ref_mean=ref_mean, ref_var=ref_var)
splt6.save(user_bname='lineages10-with-ref', add_obs=False,
           with_data_text=False, extension='.png')

# %% adding data statistics -- works only when they have been computed/exported

splt6.make_plot(report_condition=repr(condition), change_lineage_color=True,
                superimpose='all', alpha=.5, show_markers=False,
                data_statistics=True)
splt6.save(user_bname='lineages10-with-stats', add_obs=False,
           with_data_text=False, extension='.png')