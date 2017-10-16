# -*- coding: utf-8 -*-
"""
script: tutorial.py

This script is a recap of commands used in the 10 minute tutorial of tunacell
documentation.
"""
from __future__ import print_function
from builtins import input  # future package 
import argparse


def press_enter(*figs):
    """Convenient function to print figures and press Enter"""
    # when run from ipython, figure should automatically be plotted
    try:
        __IPYTHON__
    # otherwise call .plot() and wait for pressing Enter
    except NameError:
        for fig in figs:
            fig.show()
        ans = input('Press Enter to proceed...')
    return

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
    
    # Add a random sample
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
    myplot = SamplePlot(obs, colony, parser=parser)
    myplot.make_plot()
    press_enter(myplot.fig)
    myplot.save(user_bname='tutorial_sample', add_obs=False, extension='.png')
    
    # %% Statistics of the dynamics
    
    # looking at ou observable
    ou = Observable(name='growth-rate', raw='ou')
    
    from tunacell.stats.api import compute_univariate
    univariate = compute_univariate(exp, ou)
    univariate.export_text()
    
    # %% Plotting the statistics
    
    from tunacell.plotting.dynamics import plot_onepoint, plot_twopoints
    fig = plot_onepoint(univariate, show_ci=True, save=True)
    fig2 = plot_twopoints(univariate, save=True)
    press_enter(fig, fig2)