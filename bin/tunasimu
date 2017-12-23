#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
script: simurun.py

used to run numerical simulations from terminal. Command is:

    python simurun.py [[-p|--path <path>]
                       [-l|--label <exp-name>]
                       [-f|--force]
                       [-s|--samples <number-of-samples>]
                       [--alpha <growth-rate-target>]
                       [--alpha_sd <standard-deviation-for-alpha]
                       [--alpha_autocorr <autocorrelation time for alpha>]
                       [--div_lambda <lambda division parameter>]
                       [--div_size_cutoff <target size for division>]
                       [--div_use_alpha {'parameter', 'birth'}]
                       [--div_fluctuation_mode {'gamma', 'none'}]
                       [--div_relative_fluctuations <ratio sd/mean>]
                       [--initial_birth_size_mode {'fixed', 'lognormal'}]
                       [--initial_birth_size_sigma <sigma for initial birth size>]
                       [--start <starting-time-value>]
                       [--stop <stoping-time-value>]
                       [--period <time-interval-between-acquisitions>]
                      ]

All parameters have default values (see code).
"""
from __future__ import print_function
from builtins import input  # future package

import logging
import argparse
import os
import shutil
import numpy as np

from tunacell.simu.ou import OUParams, OUSimulation, SimuParams, DivisionParams, SampleInitialSize


logging.basicConfig(level=logging.INFO)



def setup_args():
    """Configure argparse parser to get all parameters for numerical simulation
    """
    # Arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # files
    parser.add_argument('--dir', type=str,
                        help='Parent directory in which simulation is stored',
                        default='~/tmptunacell')
    parser.add_argument('-l', '--label', type=str,
                        help='Label of the experiment/simulation',
                        default='simutest')
    parser.add_argument('-r', '--random-seed', type=int,
                        help='Random generator seed (use a specific number for reproducibility)',
                        default=None)
    parser.add_argument('-f', '--force',
                        help='Force writting new data (erase old simulations)',
                        action='store_true')
    parser.add_argument('-s', '--samples', type=int,
                        help='Number of simulated container samples',
                        default=100)
    # OU parameters
    parser.add_argument('--alpha', type=float,
                        help='Growth rate target (in [/min])',
                        default=np.log(2.)/60.)  # corresponding to size doubling time = 1 hour, in /min
    parser.add_argument('--alpha_sd', type=float,
                        help='Equilibrium standard deviation for alpha (in [/min])',
                        default=np.log(2.)/(10.*60.))  #1/10 of default target, in /min
    parser.add_argument('--alpha_autocorr', type=float,
                        help='Autocorrelation time for alpha (in [min]',
                        default=30.)  # in mins
    # division parameters
    parser.add_argument('--div_lambda', type=float,
                        help='Lambda parameter for interdivision timing (between 0 and 1)',
                        default=0)  # timer
    parser.add_argument('--div_size_cutoff', type=float,
                        help='Division size cutoff (in arbitrary units)',
                        default=2.)
    parser.add_argument('--div_use_alpha', type=str,
                        help='Use parameter, or birth value growth rate',
                        default='parameter')
    parser.add_argument('--div_fluctuation_mode', type=str,
                        help='[gamma, none] whether to add or not fluctuations to interdivision time',
                        default='gamma')
    parser.add_argument('--div_relative_fluctuations', type=float,
                        help='Ratio sd/predicted linear interdivision time',
                        default=.1)
    
    # initial birth size sampling
    parser.add_argument('--initial_birth_size_mode', type=str,
                        help='[fixed, lognormal]', default='fixed')
    parser.add_argument('--initial_birth_size_sigma', type=float,
                        help='sigma param for the lognormal distribution',
                        default=2.*np.log(2.))
    
    # simulation parameters
    parser.add_argument('--start', type=float,
                        help='Time at which simulation starts',
                        default=0.)
    parser.add_argument('--stop', type=float,
                        help='Time at which simulation stops',
                        default=180.)
    parser.add_argument('--period', type=float,
                        help=('Time period between two consecutive time-lapse'
                              ' acquisitions'),
                        default=5.)
    args = parser.parse_args()
#    print(vars(args).keys())
    return args


def main(args):
    """Setup numerical simulation parameters from command line arguments
    
    Parameters
    ----------
    args : :class:`argparse.Namespace` instance
        any Python object is valid as it contains following attributes:
        ['dir', 'label', 'random_seed', 'force', 'samples',
        'alpha', 'alpha_sd', 'alpha_autocorr', 'div_lambda', 'div_size_cutoff',
        'div_use_alpha', 'div_fluctuation_mode', 'div_relative_fluctuations',
        'initial_birth_size_mode', 'initial_birth_size_sigma', 'start', 'stop',
        'period']
    """


    path = os.path.abspath(os.path.expanduser(args.dir))
    if not os.path.exists(path):
        os.makedirs(path)
    print('Path: {}'.format(path))
    print('Label: {}'.format(args.label))
    # %% SET PARAMETERS FOR NUMERICAL SIMULATION
    
    # general simulation parameters
    tstart = args.start  # initial time of simulation
    tstop = args.stop  # final time of simulation (data is not recorded beyond this)
    dt = args.period  # time interval for data recording (sampling of process)
    nsamples = args.samples  # number of containers that will be created by simulation
    
    tree_per_cont = 2  # number of trees per container
    
    simuParams = SimuParams(nbr_container=nsamples,
                            nbr_colony_per_container=tree_per_cont,
                            start=tstart, stop=tstop, interval=dt)
    
    # OU parameters
    target_value = args.alpha  # average value of the random process
    spring = 1./args.alpha_autocorr  # spring constant (brings back random walker to target value)
    noise_intensity = 2. * spring * (args.alpha_sd)**2   # noise intensity for random walk
    
    ouParams = OUParams(target=target_value, spring=spring, noise=noise_intensity)
    
    # DIVISION TIMING
    divParams = DivisionParams(lambda_=args.div_lambda,
                               size_target=args.div_size_cutoff,
                               mode=args.div_fluctuation_mode,
                               relative_fluctuations=args.div_relative_fluctuations)
    
    # INITIAL BIRTH SIZE
    birthsizeParams = SampleInitialSize(size_cutoff=args.div_size_cutoff,
                                        mode=args.initial_birth_size_mode,
                                        sigma=args.initial_birth_size_sigma)
    
    # %% SIMULATION: initiate, define, run/write
    np.random.seed(42)  # fixed seed to reproduce with default results
    
    exp = OUSimulation(label=args.label,
                       simuParams=simuParams, divisionParams=divParams,
                       ouParams=ouParams, birthsizeParams=birthsizeParams)
    
    # check that experiment has been saved before
    ans = 'go'
    current_name = args.label
    exp_path = os.path.join(path, args.label)
    if args.force:
        shutil.rmtree(exp_path)
    while os.path.exists(exp_path) and ans != 'a':
        print('Experiment {} already exists.'.format(current_name))
        ans = input('Override [o], Change experiment name [c], Abort [a]: ')
        if ans == 'c':
            new_name = current_name
            while new_name == current_name:
                new_name = input('NEW NAME: ')
            exp.label = new_name
            exp_path = os.path.join(path, new_name)
        # when overriding, need to erase everything first
        elif ans == 'o':
            shutil.rmtree(exp_path)
    
    # export except if process has been aborted
    if ans != 'a':
        exp.raw_text_export(path=path)
    
    
    
if __name__ == '__main__':
    args = setup_args()
    main(args)
    
