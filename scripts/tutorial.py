# -*- coding: utf-8 -*-
"""
script: tutorial.py

This script is a recap of commands used in the 10 minute tutorial of tunacell
documentation.
"""
from __future__ import print_function
from builtins import input  # future package 
import argparse
import time


# Arguments
argparser = argparse.ArgumentParser()
argparser.add_argument('-e', '--experiment', type=str,
                       help='Path to experiment root folder',
                       default='~/tmptunacell/simutest')

args = argparser.parse_args()


if __name__ == '__main__':
    # %% Loads data from text files (exported from numerical simulations)
    from tunacell import Experiment, Parser

    path_to_exp = args.experiment
    exp = Experiment(path_to_exp)  # this defines the Experiment object
    
    print(exp)
    print('**')
    
    # %% To collect small samples, we use the Parser object
    parser = Parser(path_to_exp)
    
    # Add a known sample
    parser.add_sample(('container_079', 4))  # this one works out on default settings
    parser.add_sample(1)  # adds 1 random sample default settings have not been used
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
    print(colony)  # Python 3 compatibility --- treelib
    print('**')
    
    # %% define observable from raw data
    
    from tunacell import Observable
    obs = Observable(name='size', raw='exp_ou_int')
    print(obs)
    print(obs.as_string_table())
    print('**')
    
    # %% plotting timeseries of our colony
    
    from tunacell.plotting.samples import SamplePlot
    myplot = SamplePlot([colony, ], parser=parser)
    myplot.make_plot(obs)
    myplot.save(user_bname='tutorial_sample', add_obs=False, extension='.png')
    print('Plotting timeseries')
    myplot.fig.show()
    time.sleep(2)
    print('**')
    
    # %% Statistics of the dynamics
    
    # looking at ou observable
    ou = Observable(name='growth-rate', raw='ou')
    
    from tunacell.stats.api import compute_univariate
    univariate = compute_univariate(exp, ou)
    univariate.export_text()
    
    # %% Plotting the statistics
    
    from tunacell.plotting.dynamics import plot_onepoint, plot_twopoints
    fig = plot_onepoint(univariate, show_ci=True, save=True)
    print('Plotting one-point functions')
    fig.show()
    print('**')
    print('Plotting two-point functions')
    fig2 = plot_twopoints(univariate, save=True)
    fig2.show()
    print('**')
    ans = input('Press Enter to close all files')