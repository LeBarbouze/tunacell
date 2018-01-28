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
argparser.add_argument('-i', '--interactive',
                       help='Ask user to press Enter between parts',
                       action='store_true')
argparser.add_argument('--time', type=float,
                      help='Time per figure when non-interactive mode is on',
                      default=3)
args = argparser.parse_args()

single_plot_timing = args.time

args = argparser.parse_args()


if __name__ == '__main__':
    msg = ('==============tunacell=tutorial=============\n'
           '==                                         ==\n'
           '==             Welcome tutorial            ==\n'
           '==                                         ==\n'
           '== This tutorial presents an overview of   ==\n'
           '== the basic features of tunacell:         ==\n'
           '==   * load an experiment,                 ==\n'
           '==   * parse samples manually,             ==\n'
           '==   * define observables,                 ==\n'
           '==   * plot timeseries from samples,       ==\n'
           '==   * compute dynamic moments,            ==\n'
           '==   * plot the dynamic moments            ==\n'
           '==                                         ==\n'
           '==============tunacell=tutorial==============\n')
    print(msg)
    print()
    # Loads data from text files (exported from numerical simulations)
    print('Loading data, printing Experiment instance')
    print('******************************************')
    from tunacell import Experiment, Parser

    path_to_exp = args.experiment
    exp = Experiment(path_to_exp)  # this defines the Experiment object
    
    print(exp)
    
    if args.interactive:
        ans = input('Press Enter to proceed')
    print()
    
    # To collect small samples, we use the Parser object
    print('Parsing samples with the Parser class')
    print('*************************************')
    parser = Parser(path_to_exp)
    
    # Add a known sample
    parser.add_sample(('container_079', 4))  # this one works out on default settings
    parser.add_sample(1)  # adds 1 random sample default settings have not been used
    print(parser)
    print()
    if args.interactive:
        ans = input('Press Enter to proceed')
    print()
    
    # get cell corresponding to sample index 0
    print('Collecting cell structure for the first sample')
    print('**********************************************')
    cell = parser.get_cell(0)
    print('* Cell representation')
    print(cell)
    print()
    print('* Cell data')
    print(cell.data)
    print('* Cell data column labels')
    print(cell.data.dtype.names)
    print()
    if args.interactive:
        ans = input('Press Enter to proceed')
    print()
    
    # get colony corresponding to sample index 0
    print('Collecting colony structure for the first sample')
    print('************************************************')
    colony = parser.get_colony(0)
    print(colony)  # Python 3 compatibility --- treelib
    if args.interactive:
        ans = input('Press Enter to proceed')
    print()
    
    # define observable from raw data
    print('Define the size Observable')
    print('**************************')
    from tunacell import Observable
    obs = Observable(name='size', raw='exp_ou_int')
    print(obs)
    print(obs.as_string_table())
    if args.interactive:
        ans = input('Press Enter to proceed')
    print()
    
    # plotting timeseries of our colony
    print('Plotting timeseries of cell size for cells in our chosen colony')
    print('***************************************************************')
    from tunacell.plotting.samples import SamplePlot
    colony.decompose(seed=357)  # so that decomposition is always identical
    myplot = SamplePlot([colony, ], parser=parser)
    myplot.make_plot(obs)
    myplot.save(user_bname='tutorial_sample', add_obs=False, extension='.png')
    print('Plotting timeseries')
    myplot.fig.show()
    print()
    if args.interactive:
        ans = input('Press Enter to proceed')
    else:
        time.sleep(single_plot_timing)
    print()
    
    # Statistics of the dynamics
    print('Computing the dynamic moments of the growth-rate observable')
    print('***********************************************************')
    # looking at ou observable
    ou = Observable(name='growth-rate', raw='ou')

    from tunacell.stats.api import compute_univariate
    univariate = compute_univariate(exp, ou)
    univariate.export_text()
    print()
    if args.interactive:
        ans = input('Press Enter to proceed')
    print()
    
    # Plotting the statistics
    print('Plotting the dynamic moments of growth rate')
    print('*******************************************')
    from tunacell.plotting.dynamics import plot_onepoint, plot_twopoints
    fig = plot_onepoint(univariate, show_ci=True, save=True)
    print('* Plotting one-point functions')
    fig.show()
    print('* Plotting two-point functions')
    fig2 = plot_twopoints(univariate, save=True)
    fig2.show()
    print()
    if args.interactive:
        ans = input('Press Enter to finish tutorial (and close plots)')
    else:
        time.sleep(2*single_plot_timing)
