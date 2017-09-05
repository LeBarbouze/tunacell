#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:52:55 2017

@author: joachim
"""
from __future__ import print_function

import argparse
import os

# Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--path', type=str,
                    help='Parent directory in which simulation is stored',
                    default='~/tmptuna')
parser.add_argument('-l', '--label', type=str,
                    help='Label of the experiment/simulation',
                    default='simutest')
args = parser.parse_args()

# %% Loads data from simulations

from tuna import Parser

path_to_exp = os.path.join(os.path.abspath(os.path.expanduser(args.path)),
                           args.label)
parser = Parser(path_to_exp)
# this is equivalent to Parser('simutest') if you cd in appropriate folder

print(parser)
print('**')

# %% Add a random sample

parser.add_sample(1)
print(parser)
print('**')

# %% get cell corresponding to sample index 0

cell = parser.get_cell(0)
print(cell)
print('**')
print(cell.data)
print(cell.data.dtype.names)
print('**')

# %% get colony corresponding to sample index 0

colony = parser.get_colony(0)
colony.show()
print('**')

# %% define observable from raw data

from tuna import Observable
obs = Observable(name='size', raw='exp_ou_int', tref='root')
print(obs)
print(obs.as_string_table())
print('**')

# %% plotting timeseries of our colony

from tuna.plotting.samples import SamplePlot
myplot = SamplePlot(obs, colony, parser=parser)
myplot.make_plot()
myplot.save(user_bname='tutorial_sample', add_obs=False, extension='.png')

# %% Statistics of the dynamics

# looking at ou observable
ou = Observable(name='growth-rate', raw='ou')

from tuna.stats.api import compute_univariate_dynamics
univariate = compute_univariate(parser, ou)
univariate.export_text()

# %% Plotting the statistics

from tuna.plotting.dynamics import plot_onepoint, plot_twopoints
fig = plot_onepoint(univariate, show_ci=True, save=True)
fig2 = plot_twopoints(univariate, save=True)